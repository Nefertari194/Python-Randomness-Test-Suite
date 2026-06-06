from __future__ import print_function

import math

from bit_utils import bits_to_int
from stats_utils import chi2_sf


def approximate_entropy(bits, alpha=0.01, verbose=False):
    # 近似熵检验：比较两种相邻模式长度的分布熵差，构造卡方统计量
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    m = int(math.floor(math.log(n, 2))) - 6
    if m < 2:
        m = 2
    if m > 3:
        m = 3

    if verbose:
        print("  n         = ", n)
        print("  m         = ", m)

    phi_m = []
    for iterm in range(m, m + 2):
        # 为了统计环状模式，把前(模式长度-1)个比特拼到末尾
        padded_bits = bits + bits[0 : iterm - 1]

        counts = []
        for i in range(2**iterm):
            # 统计每个模式（当前模式长度）出现次数
            count = 0
            for j in range(n):
                if bits_to_int(padded_bits[j : j + iterm]) == i:
                    count += 1
            counts.append(count)
            if verbose:
                print("  Pattern %d of %d, count = %d" % (i + 1, 2**iterm, count))

        c_i = [float(x) / float(n) for x in counts]

        s = 0.0
        for p in c_i:
            if p > 0.0:
                s += p * math.log(p)
        phi_m.append(s)
        if verbose:
            print("  phi(%d)    = %f" % (iterm, s))

    app_en = phi_m[0] - phi_m[1]
    if verbose:
        print("  AppEn(%d)  = %f" % (m, app_en))
    chi_sq = 2.0 * float(n) * (math.log(2.0) - app_en)
    if verbose:
        print("  ChiSquare = ", chi_sq)

    df = 2**m
    # 显著性概率：卡方分布右尾概率
    p_value = chi2_sf(chi_sq, df)
    return p_value >= alpha, p_value, None

