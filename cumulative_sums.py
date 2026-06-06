from __future__ import print_function

import math

from stats_utils import erfc


def _norm_cdf(x):
    # 标准正态分布的累积分布函数
    return 0.5 * erfc(-float(x) * math.sqrt(0.5))


def _p_value(n, z):
    # 按标准公式计算累加和检验的显著性概率
    n = float(n)
    z = float(z)

    sum_a = 0.0
    startk = int(math.floor((((-n / z) + 1.0) / 4.0)))
    endk = int(math.floor((((n / z) - 1.0) / 4.0)))
    for k in range(startk, endk + 1):
        c1 = (((4.0 * k) + 1.0) * z) / math.sqrt(n)
        c2 = (((4.0 * k) - 1.0) * z) / math.sqrt(n)
        sum_a += _norm_cdf(c1) - _norm_cdf(c2)

    sum_b = 0.0
    startk = int(math.floor((((-n / z) - 3.0) / 4.0)))
    endk = int(math.floor((((n / z) - 1.0) / 4.0)))
    for k in range(startk, endk + 1):
        c1 = (((4.0 * k) + 3.0) * z) / math.sqrt(n)
        c2 = (((4.0 * k) + 1.0) * z) / math.sqrt(n)
        sum_b += _norm_cdf(c1) - _norm_cdf(c2)

    return 1.0 - sum_a + sum_b


def cumulative_sums(bits, alpha=0.01, verbose=False):
    # 累加和检验：看随机游走的最大偏离是否显著
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    # 把 0/1 映射为 -1/+1，形成随机游走
    x = [(b * 2) - 1 for b in bits]

    pos = 0
    forward_max = 0
    for e in x:
        pos += e
        if abs(pos) > forward_max:
            forward_max = abs(pos)

    pos = 0
    backward_max = 0
    for e in reversed(x):
        pos += e
        if abs(pos) > backward_max:
            backward_max = abs(pos)

    p_forward = _p_value(n, forward_max)
    p_backward = _p_value(n, backward_max)

    if verbose:
        print("  p_forward  = ", p_forward)
        print("  p_backward = ", p_backward)

    plist = [p_forward, p_backward]
    return (p_forward >= alpha and p_backward >= alpha), None, plist

