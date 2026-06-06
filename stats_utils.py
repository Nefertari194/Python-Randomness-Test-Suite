from __future__ import print_function

import math

from gamma_functions import gammaincc


try:
    # 优先使用科学计算库的实现：精度/稳定性通常更好
    from scipy import stats as _scipy_stats
except Exception:
    _scipy_stats = None


try:
    # 用于特殊函数的高精度实现
    from scipy import special as _scipy_special
except Exception:
    _scipy_special = None


def erfc(x):
    # 互补误差函数：很多测试的正态分布尾概率会用到
    if _scipy_special is not None:
        return float(_scipy_special.erfc(x))
    return math.erfc(x)


def chi2_sf(chi2, df):
    # 卡方分布右尾概率（自由度由参数指定）
    if _scipy_stats is not None:
        return float(_scipy_stats.chi2.sf(chi2, df))
    # 没有科学计算库时，用不完全伽马函数的互补形式计算
    return float(gammaincc(float(df) / 2.0, float(chi2) / 2.0))


def normal_two_sided_p_from_z(z):
    # 标准正态双侧显著性概率
    return float(erfc(abs(float(z)) / math.sqrt(2.0)))

