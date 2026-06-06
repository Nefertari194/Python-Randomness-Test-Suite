from __future__ import print_function

import math

from stats_utils import chi2_sf


def _prob(m, i):
    m8 = [0.2148, 0.3672, 0.2305, 0.1875]
    m128 = [0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124]
    m512 = [0.1170, 0.2460, 0.2523, 0.1755, 0.1027, 0.1124]
    m1000 = [0.1307, 0.2437, 0.2452, 0.1714, 0.1002, 0.1088]
    m10000 = [0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727]
    if m == 8:
        return m8[i]
    if m == 128:
        return m128[i]
    if m == 512:
        return m512[i]
    if m == 1000:
        return m1000[i]
    return m10000[i]


def longest_run_ones(bits, alpha=0.01, verbose=False):
    # 最长连续 1 检验：把序列分块，统计每块中最长连续 1 的长度落在各区间的频数，再做卡方检验
    n = len(bits)
    if n < 128:
        return False, None, {"reason": "长度不够：这个测试至少要 128 位"}
    if n < 6272:
        m = 8
    elif n < 750000:
        m = 128
    else:
        m = 10000

    if m == 8:
        k = 3
        min_n = 16
    elif m == 128:
        k = 5
        min_n = 49
    else:
        k = 6
        min_n = 75

    n_blocks = int(math.floor(float(n) / float(m)))
    if n_blocks < min_n:
        return False, None, {"reason": "有效块数太少（N=%d < %d）" % (n_blocks, min_n)}

    v = [0, 0, 0, 0, 0, 0, 0]
    for i in range(n_blocks):
        block = bits[i * m : (i + 1) * m]
        run = 0
        longest = 0
        for b in block:
            if b == 1:
                run += 1
                if run > longest:
                    longest = run
            else:
                run = 0

        if m == 8:
            if longest <= 1:
                v[0] += 1
            elif longest == 2:
                v[1] += 1
            elif longest == 3:
                v[2] += 1
            else:
                v[3] += 1
        elif m == 128:
            if longest <= 4:
                v[0] += 1
            elif longest == 5:
                v[1] += 1
            elif longest == 6:
                v[2] += 1
            elif longest == 7:
                v[3] += 1
            elif longest == 8:
                v[4] += 1
            else:
                v[5] += 1
        else:
            if longest <= 10:
                v[0] += 1
            elif longest == 11:
                v[1] += 1
            elif longest == 12:
                v[2] += 1
            elif longest == 13:
                v[3] += 1
            elif longest == 14:
                v[4] += 1
            elif longest == 15:
                v[5] += 1
            else:
                v[6] += 1

    chi_sq = 0.0
    for i in range(k + 1):
        p_i = _prob(m, i)
        # 期望频数 = 块数 × 对应区间的理论概率
        chi_sq += ((v[i] - (n_blocks * p_i)) ** 2) / (n_blocks * p_i)

    if verbose:
        print("  n = " + str(n))
        print("  K = " + str(k))
        print("  M = " + str(m))
        print("  N = " + str(n_blocks))
        print("  chi_sq = " + str(chi_sq))

    p = chi2_sf(chi_sq, k)
    return p >= alpha, p, None

