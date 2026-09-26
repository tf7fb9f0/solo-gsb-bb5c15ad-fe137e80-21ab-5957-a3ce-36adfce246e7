"""同位素质量偏移分布的卷积引擎。

分布表示：``{offset: {"p": Decimal, "mass": Decimal}}``

* ``p``    —— 该整数质量偏移下所有原子配置的概率之和；
* ``mass`` —— 同一偏移下概率加权的精确质量累加值（Σ p_config · m_config），
  质心 = ``mass / p``。

始终进行精确的完整卷积，相同偏移合并概率与质量加权值，绝不在中途丢弃
“低概率”项；所有计算使用 Decimal，避免浮点误差。
"""

from decimal import Decimal, localcontext
from typing import Dict, List, Tuple

from .schemas import ABUNDANCE_SUM_TOLERANCE, EnvelopeRequest, ElementIn

# 卷积内部的 Decimal 有效数字精度。最多 60 个原子、每项 241 个，普通输入
# 远不需要这么多位数，取 80 位足以保证任何非零中间项都不会被舍入为 0。
DECIMAL_PRECISION = 80

# 分布条目的键名，便于保持类型清晰
_PROB = "p"
_MASS = "mass"
Entry = Dict[str, Decimal]
Distribution = Dict[int, Entry]


def _single_atom_distribution(element: ElementIn) -> Tuple[Distribution, Decimal]:
    """构造单原子分布，并归一化丰度、确定该元素的精确质量基准。

    返回 ``(分布, baseline_exact_mass)``；同一质量数的多个同位素先合并。
    ``baseline_exact_mass`` 为最小质量数下各同位素的概率加权平均精确质量，
    即该元素单原子在偏移 0 处的质心。
    """
    base_mass_number = min(iso.mass_number for iso in element.isotopes)

    abundance_sum = sum((iso.abundance for iso in element.isotopes), Decimal(0))
    # 校验阶段已保证合计与 1 的偏差不超过 ABUNDANCE_SUM_TOLERANCE，这里除以
    # 实际合计做归一化，使最终概率总和严格为 1，且不改变各峰的相对比例。
    normalized: Dict[int, Tuple[Decimal, Decimal]] = {}
    base_group_mass = Decimal(0)
    base_group_prob = Decimal(0)
    for iso in element.isotopes:
        offset = iso.mass_number - base_mass_number
        probability = iso.abundance / abundance_sum
        if offset == 0:
            base_group_mass += probability * iso.exact_mass
            base_group_prob += probability
        if offset in normalized:
            prev_prob, prev_mass = normalized[offset]
            prob = prev_prob + probability
            # 同质量数（同 offset）同位素先合并：质量按概率加权
            mass = prev_mass + probability * iso.exact_mass
            normalized[offset] = (prob, mass)
        else:
            normalized[offset] = (probability, probability * iso.exact_mass)

    distribution: Distribution = {}
    for offset, (probability, weighted_mass) in normalized.items():
        distribution[offset] = {_PROB: probability, _MASS: weighted_mass}
    # 基准质量 = 最小质量数分组内按丰度条件加权的平均精确质量，即单原子
    # 偏移 0 峰的质心（保证 baseline_mass 与第一个峰的 centroid_mass 一致）。
    baseline_exact_mass = base_group_mass / base_group_prob
    return distribution, baseline_exact_mass


def _convolve(left: Distribution, right: Distribution) -> Distribution:
    """两个独立分布做卷积，相同偏移合并概率与概率加权质量。

    合并规则（偏移相同的两个原子配置，配置质量为绝对精确质量）：
        p_new = p₁p₂
        W_new = p₁p₂·(m₁ + m₂) = W₁·p₂ + p₁·W₂，其中 W = p·m
    质心在最终输出时由 W / p 求出。
    """
    result: Distribution = {}
    for off_left, entry_left in left.items():
        p_left = entry_left[_PROB]
        w_left = entry_left[_MASS]
        for off_right, entry_right in right.items():
            p_right = entry_right[_PROB]
            w_right = entry_right[_MASS]
            probability = p_left * p_right
            weighted_mass = w_left * p_right + p_left * w_right
            offset = off_left + off_right
            if offset in result:
                result[offset][_PROB] += probability
                result[offset][_MASS] += weighted_mass
            else:
                result[offset] = {_PROB: probability, _MASS: weighted_mass}
    return result


def _power_distribution(single: Distribution, count: int) -> Distribution:
    """通过平方倍增求同一种元素 ``count`` 个原子的独立抽样分布（完整卷积）。"""
    if count == 0:
        return {0: {_PROB: Decimal(1), _MASS: Decimal(0)}}

    result: Distribution = {0: {_PROB: Decimal(1), _MASS: Decimal(0)}}
    base = single
    remaining = count
    while remaining > 0:
        if remaining & 1:
            result = _convolve(result, base)
        remaining >>= 1
        if remaining > 0:
            base = _convolve(base, base)
    return result


def _assert_normalized_input(elements: List[ElementIn]) -> None:
    """防御性校验：丰度合计与 1 的偏差需在容差内（schemas 层已校验）。"""
    for element in elements:
        abundance_sum = sum((iso.abundance for iso in element.isotopes), Decimal(0))
        if abs(abundance_sum - Decimal(1)) > ABUNDANCE_SUM_TOLERANCE:
            raise ValueError("同位素丰度合计必须为 1")


def compute_envelope(request: EnvelopeRequest) -> dict:
    """计算分子的同位素峰包络。

    以各元素最小质量数为整数质量偏移基准；每种元素的 ``atom_count`` 个原子
    独立抽样后做完整卷积，元素之间再卷积；相同偏移合并概率并计算精确质量
    质心。返回可直接序列化为 :class:`EnvelopeResponse` 的字典。
    """
    with localcontext() as context:
        context.prec = DECIMAL_PRECISION

        _assert_normalized_input(request.elements)

        molecule: Distribution = {0: {_PROB: Decimal(1), _MASS: Decimal(0)}}
        baseline_mass = Decimal(0)

        for element in request.elements:
            single, element_baseline = _single_atom_distribution(element)
            baseline_mass += element_baseline * element.atom_count
            element_distribution = _power_distribution(single, element.atom_count)
            molecule = _convolve(molecule, element_distribution)

        peaks = []
        total_probability = Decimal(0)
        for offset in sorted(molecule.keys()):
            probability = molecule[offset][_PROB]
            # 不做任何概率阈值截断；只有数学上真正为 0 的项才不输出。
            if probability == 0:
                continue
            total_probability += probability
            # W = Σ p_config·m_config 已用绝对精确质量累加，故 W/p 即该峰
            # 的概率加权精确质量质心；baseline_mass 仅作为基准单独返回。
            centroid = molecule[offset][_MASS] / probability
            peaks.append(
                {
                    "offset": offset,
                    "probability": probability,
                    "centroid_mass": centroid,
                }
            )

        return {
            "baseline_mass": baseline_mass,
            "total_probability": total_probability,
            "peaks": peaks,
        }
