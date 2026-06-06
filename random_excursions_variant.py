from __future__ import print_function

import math

from stats_utils import erfc


def random_excursions_variant(bits, alpha=0.01, verbose=False, min_J=500):
    # 随机游走变体：统计各状态访问次数与期望的偏差，转换为正态近似的显著性概率
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    # 0/1 -> -1/+1
    x = [(bit * 2) - 1 for bit in bits]

    pos = 0
    s = []
    for e in x:
        pos += e
        s.append(pos)
    sprime = [0] + s + [0]

    j = 0
    for value in sprime[1:]:
        if value == 0:
            j += 1
    if verbose:
        print("J=", j)
    if j <= 0:
        return False, None, {"reason": "循环数 J 太小（%d）" % j}
    if j < int(min_J):
        return False, None, {"reason": "循环数 J 太小（%d < %d），这个测试结果不靠谱" % (j, int(min_J))}

    count = [0 for _ in range(19)]
    # 只统计状态 -9..9 的访问次数
    for value in sprime:
        if abs(value) < 10:
            count[value + 9] += 1

    success = True
    plist = []
    for state in range(-9, 10):
        if state == 0:
            continue
        top = abs(count[state + 9] - j)
        bottom = math.sqrt(2.0 * float(j) * ((4.0 * abs(state)) - 2.0))
        if bottom == 0.0:
            return False, None, {"reason": "分母为 0，无法计算"}
        p = float(erfc(float(top) / float(bottom)))
        plist.append(p)
        if p < alpha:
            success = False
        if verbose:
            err = " Not Random" if p < alpha else ""
            print("x = %1.0f\t count=%d\tp = %f %s" % (state, count[state + 9], p, err))

    return success, None, plist

