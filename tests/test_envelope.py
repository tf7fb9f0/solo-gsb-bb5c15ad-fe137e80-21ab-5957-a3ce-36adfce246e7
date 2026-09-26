"""端到端测试：覆盖数值正确性、完整卷积与 422 校验。"""

from decimal import Decimal

from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

CARBON = [
    {"mass_number": 12, "exact_mass": "12.000000", "abundance": 0.9893},
    {"mass_number": 13, "exact_mass": "13.003355", "abundance": 0.0107},
]
OXYGEN = [
    {"mass_number": 16, "exact_mass": "15.994915", "abundance": 0.99757},
    {"mass_number": 17, "exact_mass": "16.999132", "abundance": 0.00038},
    {"mass_number": 18, "exact_mass": "17.999161", "abundance": 0.00205},
]
# 质量数跨度恰为 4：允许
WIDE = [
    {"mass_number": 10, "exact_mass": "10.0", "abundance": 0.5},
    {"mass_number": 14, "exact_mass": "14.0", "abundance": 0.5},
]


def _envelope(elements):
    return client.post("/envelope", json={"elements": elements})


def test_single_element_c1_methane_like_binomial_probabilities():
    resp = _envelope([{"atom_count": 2, "isotopes": CARBON}])
    assert resp.status_code == 200
    data = resp.json()
    assert [p["offset"] for p in data["peaks"]] == [0, 1, 2]

    p0, p1, p2 = (p["probability"] for p in data["peaks"])
    assert abs(p0 - 0.9893**2) < 1e-12
    assert abs(p1 - 2 * 0.9893 * 0.0107) < 1e-12
    assert abs(p2 - 0.0107**2) < 1e-12

    assert abs(data["total_probability"] - 1.0) < 1e-12
    assert data["baseline_mass"] == 24.0


def test_centroid_masses_are_probability_weighted():
    resp = _envelope([{"atom_count": 2, "isotopes": CARBON}])
    peaks = {p["offset"]: p for p in resp.json()["peaks"]}

    p, q = Decimal("0.9893"), Decimal("0.0107")
    m12, m13 = Decimal("12.000000"), Decimal("13.003355")
    expected = {
        0: 2 * m12,
        1: m12 + m13,
        2: 2 * m13,
    }
    for offset, mass in expected.items():
        assert abs(Decimal(str(peaks[offset]["centroid_mass"])) - mass) < Decimal("1e-9")


def test_multi_element_convolution_co2():
    resp = _envelope(
        [
            {"symbol": "C", "atom_count": 1, "isotopes": CARBON},
            {"symbol": "O", "atom_count": 2, "isotopes": OXYGEN},
        ]
    )
    assert resp.status_code == 200
    peaks = {p["offset"]: p for p in resp.json()["peaks"]}

    pc, po = 0.9893, 0.99757
    # offset 0：C12 且两个 O 都是 O16
    assert abs(peaks[0]["probability"] - pc * po**2) < 1e-12
    # offset 1：C13(O16)2 或 C12·O16·O17（两个位置）
    assert abs(peaks[1]["probability"] - (0.0107 * po**2 + pc * 2 * 0.99757 * 0.00038)) < 1e-12

    assert abs(resp.json()["total_probability"] - 1.0) < 1e-12


def test_peaks_sorted_by_offset_and_baseline():
    resp = _envelope(
        [
            {"symbol": "C", "atom_count": 3, "isotopes": CARBON},
            {"symbol": "O", "atom_count": 1, "isotopes": OXYGEN},
        ]
    )
    data = resp.json()
    offsets = [p["offset"] for p in data["peaks"]]
    assert offsets == sorted(offsets)
    assert offsets[0] == 0
    # 3 个 C 最大偏移 3，1 个 O 最大偏移 2
    assert offsets[-1] == 5
    assert abs(data["baseline_mass"] - (3 * 12.0 + 15.994915)) < 1e-12


def test_complete_convolution_keeps_all_nonzero_terms():
    """20 个稀同位素丰度 1e-6 的原子：最大偏移 20，且每个偏移都必须保留。"""
    rare = [
        {"mass_number": 1, "exact_mass": "1.0", "abundance": "0.999999"},
        {"mass_number": 2, "exact_mass": "2.0", "abundance": "0.000001"},
    ]
    resp = _envelope([{"symbol": "X", "atom_count": 20, "isotopes": rare}])
    assert resp.status_code == 200
    data = resp.json()
    offsets = [p["offset"] for p in data["peaks"]]
    # 完整卷积：0～20 共 21 个非零峰，没有因为概率低而提前丢弃
    assert offsets == list(range(21))
    # offset 20 的概率为 1e-120，在 Decimal 内部保留；输出转 float 后下溢为 0，
    # 因此这里改用内部引擎验证非零项的保留。
    from app.envelope import compute_envelope
    from app.schemas import EnvelopeRequest

    request = EnvelopeRequest.model_validate(
        {"elements": [{"symbol": "X", "atom_count": 20, "isotopes": rare}]}
    )
    result = compute_envelope(request)
    assert len(result["peaks"]) == 21
    assert result["peaks"][-1]["offset"] == 20
    assert result["peaks"][-1]["probability"] == Decimal("1e-120")
    assert result["total_probability"] == Decimal(1)


def test_span_four_is_allowed_but_five_is_rejected():
    ok = _envelope([{"atom_count": 1, "isotopes": WIDE}])
    assert ok.status_code == 200
    assert ok.json()["peaks"][-1]["offset"] == 4

    too_wide = [
        {"mass_number": 10, "exact_mass": "10.0", "abundance": 0.5},
        {"mass_number": 15, "exact_mass": "15.0", "abundance": 0.5},
    ]
    bad = _envelope([{"atom_count": 1, "isotopes": too_wide}])
    assert bad.status_code == 422
    loc = [" -> ".join(str(x) for x in e["loc"]) for e in bad.json()["detail"]]
    assert any("elements" in x and "isotopes" in x for x in loc)


def test_abundance_must_be_nonnegative_and_sum_to_one():
    negative = [
        {"mass_number": 12, "exact_mass": "12.0", "abundance": 1.2},
        {"mass_number": 13, "exact_mass": "13.0", "abundance": -0.2},
    ]
    resp = _envelope([{"atom_count": 1, "isotopes": negative}])
    assert resp.status_code == 422
    locs = [".".join(str(x) for x in e["loc"]) for e in resp.json()["detail"]]
    assert any("abundance" in loc for loc in locs)

    bad_sum = [
        {"mass_number": 12, "exact_mass": "12.0", "abundance": 0.9},
        {"mass_number": 13, "exact_mass": "13.0", "abundance": 0.05},
    ]
    resp = _envelope([{"atom_count": 1, "isotopes": bad_sum}])
    assert resp.status_code == 422
    locs = [".".join(str(x) for x in e["loc"]) for e in resp.json()["detail"]]
    assert any("isotopes" in loc for loc in locs)


def test_total_atom_count_bounds():
    resp_zero = _envelope([{"atom_count": 0, "isotopes": CARBON}])
    assert resp_zero.status_code == 422

    resp_big = _envelope([{"atom_count": 61, "isotopes": CARBON}])
    assert resp_big.status_code == 422

    resp_ok = _envelope(
        [
            {"atom_count": 30, "isotopes": CARBON},
            {"atom_count": 30, "isotopes": OXYGEN},
        ]
    )
    assert resp_ok.status_code == 200
    assert len(resp_ok.json()["peaks"]) == 91  # offset 0～90


def test_zero_abundance_isotopes_do_not_create_peaks():
    isotopes = [
        {"mass_number": 12, "exact_mass": "12.0", "abundance": 1.0},
        {"mass_number": 13, "exact_mass": "13.003355", "abundance": 0.0},
    ]
    resp = _envelope([{"atom_count": 3, "isotopes": isotopes}])
    assert resp.status_code == 200
    data = resp.json()
    assert [p["offset"] for p in data["peaks"]] == [0]
    assert data["peaks"][0]["probability"] == 1.0
    assert data["peaks"][0]["centroid_mass"] == 36.0


def test_basic_type_and_missing_field_errors_have_locations():
    resp = client.post("/envelope", json={"elements": [{"atom_count": -1, "isotopes": CARBON}]})
    assert resp.status_code == 422
    locs = [".".join(str(x) for x in e["loc"]) for e in resp.json()["detail"]]
    assert any("atom_count" in loc for loc in locs)

    resp = client.post("/envelope", json={"elements": []})
    assert resp.status_code == 422

    resp = client.post("/envelope", json={})
    assert resp.status_code == 422
    locs = [".".join(str(x) for x in e["loc"]) for e in resp.json()["detail"]]
    assert any("elements" in loc for loc in locs)
