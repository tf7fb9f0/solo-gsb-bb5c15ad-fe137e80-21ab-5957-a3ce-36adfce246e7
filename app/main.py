"""FastAPI 应用入口。"""

from decimal import Decimal

from fastapi import FastAPI

from .envelope import compute_envelope
from .schemas import EnvelopeRequest, EnvelopeResponse, PeakOut

app = FastAPI(
    title="同位素峰包络计算 API",
    version="1.0.0",
    description=(
        "根据各元素的原子数及同位素质量数、精确质量与丰度，"
        "通过独立抽样完整卷积计算分子的质量偏移分布与精确质量质心。"
    ),
)


@app.get("/health", tags=["health"], summary="健康检查")
def health() -> dict:
    return {"status": "ok"}


@app.post(
    "/envelope",
    response_model=EnvelopeResponse,
    tags=["envelope"],
    summary="计算同位素峰包络",
)
def envelope(request: EnvelopeRequest) -> EnvelopeResponse:
    result = compute_envelope(request)

    peaks = [
        PeakOut(
            offset=peak["offset"],
            # JSON number 为 IEEE-754 双精度；卷积全程使用 Decimal，仅输出时转换
            probability=float(peak["probability"]),
            centroid_mass=float(peak["centroid_mass"]),
        )
        for peak in result["peaks"]
    ]
    return EnvelopeResponse(
        baseline_mass=float(result["baseline_mass"]),
        total_probability=float(result["total_probability"]),
        peaks=peaks,
    )
