from __future__ import print_function

import math

from stats_utils import chi2_sf


def berelekamp_massey(bits):
    # 伯勒坎普-梅西算法：求二进制序列的最小线性反馈移位寄存器长度（线性复杂度）
    n = len(bits)
    b = [0 for _ in bits]
    c = [0 for _ in bits]
    b[0] = 1
    c[0] = 1

    l = 0
    m = -1
    nn = 0
    while nn < n:
        d = bits[nn]
        for i in range(1, l + 1):
            d = d ^ (c[i] & bits[nn - i])
        if d != 0:
            t = c[:]
            for i in range(0, n - nn + m):
                c[nn - m + i] = c[nn - m + i] ^ b[i]
            if l <= (nn / 2):
                l = nn + 1 - l
                m = nn
                b = t
        nn = nn + 1

    return l, c[0:l]


def linear_complexity(bits, patternlen=None, alpha=0.01, verbose=False):
    # 线性复杂度检验：把序列分块，计算每块线性复杂度，落桶后做卡方检验
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    if patternlen is not None:
        m = int(patternlen)
    else:
        if n < 1000000:
            return False, None, {"reason": "长度不够：这个测试默认要至少 10^6 位"}
        m = 512

    k = 6
    n_blocks = int(math.floor(float(n) / float(m)))
    if n_blocks <= 0:
        return False, None, {"reason": "分块后块数为 0"}

    if verbose:
        print("  M = ", m)
        print("  N = ", n_blocks)
        print("  K = ", k)

    lc = []
    for i in range(n_blocks):
        x = bits[(i * m) : ((i + 1) * m)]
        lc.append(berelekamp_massey(x)[0])

    a = float(m) / 2.0
    b = ((((-1) ** (m + 1)) + 9.0)) / 36.0
    c = ((m / 3.0) + (2.0 / 9.0)) / (2**m)
    mu = a + b - c

    t_list = []
    for i in range(n_blocks):
        x = ((-1.0) ** m) * (lc[i] - mu) + (2.0 / 9.0)
        t_list.append(x)

    v = [0, 0, 0, 0, 0, 0, 0]
    # 按标准区间把统计量分到 7 个桶
    for t in t_list:
        if t <= -2.5:
            v[0] += 1
        elif t <= -1.5:
            v[1] += 1
        elif t <= -0.5:
            v[2] += 1
        elif t <= 0.5:
            v[3] += 1
        elif t <= 1.5:
            v[4] += 1
        elif t <= 2.5:
            v[5] += 1
        else:
            v[6] += 1

    pi = [0.010417, 0.03125, 0.125, 0.5, 0.25, 0.0625, 0.020833]
    chi_sq = 0.0
    for i in range(k + 1):
        chi_sq += ((v[i] - (n_blocks * pi[i])) ** 2.0) / (n_blocks * pi[i])

    if verbose:
        print("  chisq = ", chi_sq)

    p = chi2_sf(chi_sq, k)
    if verbose:
        print("  P = ", p)
    return p >= alpha, p, None

