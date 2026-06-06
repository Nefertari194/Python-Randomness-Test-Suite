from __future__ import print_function

import math

from bit_utils import count_ones_zeroes
from stats_utils import chi2_sf


def frequency_block(bits, alpha=0.01, verbose=False, block_size=None, num_blocks=None):
    # 分块频数检验：把序列切成若干块，每块看 1 的比例是否偏离 1/2，并做卡方检验
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}
    if n < 100:
        return False, None, {"reason": "长度不够：这个测试至少要 100 位"}

    if block_size is not None:
        m = int(block_size)
        if m <= 0:
            return False, None, {"reason": "块大小得是正数"}
        n_blocks = int(math.floor(float(n) / float(m)))
    elif num_blocks is not None:
        n_blocks = int(num_blocks)
        if n_blocks <= 0:
            return False, None, {"reason": "块数得是正数"}
        m = int(math.floor(float(n) / float(n_blocks)))
    else:
        m = 20
        n_blocks = int(math.floor(float(n) / float(m)))
        if n_blocks > 99:
            n_blocks = 99
            m = int(math.floor(float(n) / float(n_blocks)))

    if n_blocks <= 0 or m <= 0:
        return False, None, {"reason": "分块参数不合适"}
    if n_blocks * m > n:
        return False, None, {"reason": "块数和块大小加起来超过了序列长度"}

    if verbose:
        print("  n = %d" % n)
        print("  N = %d" % n_blocks)
        print("  M = %d" % m)

    chi_sq = 0.0
    for i in range(n_blocks):
        block = bits[i * m : (i + 1) * m]
        _, ones = count_ones_zeroes(block)
        pi = float(ones) / float(m)
        # 每块贡献：4×块长度×(比例-1/2)^2
        chi_sq += 4.0 * float(m) * ((pi - 0.5) ** 2)

    # 显著性概率：卡方分布右尾概率（自由度等于块数）
    p = chi2_sf(chi_sq, n_blocks)
    return p >= alpha, p, None

