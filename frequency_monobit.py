from __future__ import print_function

import math

from bit_utils import count_ones_zeroes
from stats_utils import normal_two_sided_p_from_z


def frequency_monobit(bits, alpha=0.01, verbose=False):
    # 单比特频数检验：看 0/1 个数是否显著偏离 1/2
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    zeroes, ones = count_ones_zeroes(bits)
    s_obs = abs(ones - zeroes)
    if verbose:
        print("  1 的个数 = %d" % ones)
        print("  0 的个数 = %d" % zeroes)

    z = float(s_obs) / (math.sqrt(float(n)) * math.sqrt(2.0))
    # 显著性概率使用标准正态双侧尾概率
    p = normal_two_sided_p_from_z(z)
    return p >= alpha, p, None

