"""请求 / 响应模型与输入校验。"""

from decimal import Decimal
from typing import List, Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

# 丰度合计与 1 的最大允许偏差。普通 JSON 浮点（如 0.9893 + 0.0107）经
# Decimal(float) 解析后会带有 ~1e-16 量级的误差，因此给一个极小的容差；
# 进入卷积前会对丰度做归一化，保证输出概率总和严格为 1。
ABUNDANCE_SUM_TOLERANCE = Decimal("1e-9")

MIN_TOTAL_ATOMS = 1
MAX_TOTAL_ATOMS = 60
MAX_MASS_NUMBER_SPAN = 4


class IsotopeIn(BaseModel):
    """单个同位素：整数质量数、十进制精确质量、丰度。"""

    model_config = ConfigDict(extra="forbid")

    mass_number: int = Field(..., ge=1, description="整数质量数，必须为正整数")
    exact_mass: Decimal = Field(..., description="十进制精确质量（u），必须为有限正数")
    abundance: Decimal = Field(..., description="丰度，范围 [0, 1]，同一元素内合计为 1")

    @field_validator("exact_mass")
    @classmethod
    def _check_exact_mass(cls, value: Decimal) -> Decimal:
        # 先判有限再比较大小：Decimal NaN 参与排序会抛 InvalidOperation。
        if not value.is_finite() or value <= 0:
            raise ValueError("exact_mass 必须为有限正数")
        return value

    @field_validator("abundance")
    @classmethod
    def _check_abundance(cls, value: Decimal) -> Decimal:
        if not value.is_finite() or value < 0 or value > 1:
            raise ValueError("abundance 必须为 0 到 1 之间的有限数")
        return value


class ElementIn(BaseModel):
    """一种元素：非负整数原子数及其同位素列表。"""

    model_config = ConfigDict(extra="forbid")

    symbol: Optional[str] = Field(default=None, max_length=16, description="可选元素符号，仅用于阅读")
    atom_count: int = Field(..., ge=0, description="该元素的原子数，非负整数")
    isotopes: List[IsotopeIn] = Field(..., min_length=1, description="该元素的同位素列表，至少 1 个")

    @field_validator("isotopes")
    @classmethod
    def _check_isotopes(cls, value: List[IsotopeIn]) -> List[IsotopeIn]:
        abundance_sum = sum((iso.abundance for iso in value), Decimal(0))
        if abs(abundance_sum - Decimal(1)) > ABUNDANCE_SUM_TOLERANCE:
            raise ValueError(
                f"同一元素各同位素丰度合计必须为 1，实际合计为 {abundance_sum}"
            )

        min_mass = min(iso.mass_number for iso in value)
        max_mass = max(iso.mass_number for iso in value)
        if max_mass - min_mass > MAX_MASS_NUMBER_SPAN:
            raise ValueError(
                f"同一元素同位素质量数跨度不能超过 {MAX_MASS_NUMBER_SPAN}，"
                f"实际跨度为 {max_mass - min_mass}（{min_mass}～{max_mass}）"
            )
        return value


class EnvelopeRequest(BaseModel):
    """POST /envelope 请求体。"""

    model_config = ConfigDict(extra="forbid")

    elements: List[ElementIn] = Field(..., min_length=1, description="分子所含元素列表，至少 1 种")

    @model_validator(mode="after")
    def _check_total_atom_count(self) -> "EnvelopeRequest":
        total = sum(element.atom_count for element in self.elements)
        if not MIN_TOTAL_ATOMS <= total <= MAX_TOTAL_ATOMS:
            raise ValueError(
                f"所有元素原子数合计必须在 {MIN_TOTAL_ATOMS}～{MAX_TOTAL_ATOMS} 之间，"
                f"实际合计为 {total}"
            )
        return self


class PeakOut(BaseModel):
    """单个质量偏移峰。"""

    offset: int = Field(..., description="相对最小质量数基准的整数质量偏移")
    probability: float = Field(..., description="该偏移的合并概率")
    centroid_mass: float = Field(..., description="概率加权的精确质量质心（u）")


class EnvelopeResponse(BaseModel):
    """POST /envelope 响应体。"""

    baseline_mass: float = Field(..., description="所有原子均取最小质量数同位素时的精确质量（u）")
    total_probability: float = Field(..., description="全部峰的概率总和，归一化后应为 1")
    peaks: List[PeakOut] = Field(..., description="按 offset 升序排列的非零概率峰列表")
