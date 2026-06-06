from __future__ import print_function

import math

from bit_utils import count_ones_zeroes
from stats_utils import normal_two_sided_p_from_z


def runs(bits, alpha=0.01, verbose=False):
    # 游程检验：统计 0/1 切换次数是否符合随机序列的期望
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    _, ones = count_ones_zeroes(bits)
    prop = float(ones) / float(n)
    if verbose:
        print("  prop ", prop)

    tau = 2.0 / math.sqrt(float(n))
    if verbose:
        print("  tau ", tau)

    if abs(prop - 0.5) > tau:
        # 标准里要求 0/1 比例不能偏太多，否则统计量假设不成立
        return False, None, {"reason": "0/1 比例偏得有点多，这个测试按标准就不做了"}

    v_obs = 1.0
    for i in range(n - 1):
        if bits[i] != bits[i + 1]:
            v_obs += 1.0

    if verbose:
        print("  vobs ", v_obs)

    denom = 2.0 * math.sqrt(2.0 * float(n)) * prop * (1.0 - prop)
    if denom == 0.0:
        return False, None, {"reason": "分母为 0，无法计算"}
    z = abs(v_obs - (2.0 * float(n) * prop * (1.0 - prop))) / denom
    p = normal_two_sided_p_from_z(z)
    return p >= alpha, p, None

