from __future__ import print_function

import math

import gf2matrix

from stats_utils import chi2_sf


def binary_matrix_rank(bits, m=32, q=32, alpha=0.01, verbose=False):
    # 二进制矩阵秩检验：把序列切成“行×列”的矩阵块，统计满秩/秩-1/其它 的分布并做卡方检验
    n = len(bits)
    m = int(m)
    q = int(q)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}
    if m <= 0 or q <= 0:
        return False, None, {"reason": "矩阵尺寸必须是正数"}

    n_blocks = int(math.floor(float(n) / float(m * q)))
    if verbose:
        print("  块数 %d" % n_blocks)
        print("  用掉的 bit 数: %d" % (n_blocks * m * q))
        print("  丢掉的 bit 数: %d" % (n - (n_blocks * m * q)))

    if n_blocks < 38:
        return False, None, {"reason": "块数太少：至少要 38 个块才能做这个测试"}

    # 按标准公式计算理论概率：满秩、秩-1、其它
    r = m
    product = 1.0
    for i in range(r):
        upper1 = 1.0 - (2.0 ** (i - q))
        upper2 = 1.0 - (2.0 ** (i - m))
        lower = 1.0 - (2.0 ** (i - r))
        product *= (upper1 * upper2) / lower
    fr_prob = product * (2.0 ** ((r * (q + m - r)) - (m * q)))

    r = m - 1
    product = 1.0
    for i in range(r):
        upper1 = 1.0 - (2.0 ** (i - q))
        upper2 = 1.0 - (2.0 ** (i - m))
        lower = 1.0 - (2.0 ** (i - r))
        product *= (upper1 * upper2) / lower
    frm1_prob = product * (2.0 ** ((r * (q + m - r)) - (m * q)))

    lr_prob = 1.0 - (fr_prob + frm1_prob)

    fm = 0
    fmm = 0
    remainder = 0
    for blknum in range(n_blocks):
        block = bits[blknum * (m * q) : (blknum + 1) * (m * q)]
        # 二元伽罗瓦域(2)上的矩阵，行变换求秩
        matrix = gf2matrix.matrix_from_bits(m, q, block, blknum)
        rank = gf2matrix.rank(m, q, matrix, blknum)

        if rank == m:
            fm += 1
        elif rank == (m - 1):
            fmm += 1
        else:
            remainder += 1

    chi_sq = ((fm - (fr_prob * n_blocks)) ** 2) / (fr_prob * n_blocks)
    chi_sq += ((fmm - (frm1_prob * n_blocks)) ** 2) / (frm1_prob * n_blocks)
    chi_sq += ((remainder - (lr_prob * n_blocks)) ** 2) / (lr_prob * n_blocks)
    # 该卡方统计量自由度为 2
    p = chi2_sf(chi_sq, 2)

    if verbose:
        print("  满秩矩阵个数  = ", fm)
        print("  秩少 1 的个数 = ", fmm)
        print("  其它情况个数  = ", remainder)
        print("  卡方值 = ", chi_sq)

    return p >= alpha, p, None

