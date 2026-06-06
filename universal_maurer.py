from __future__ import print_function

import math

from stats_utils import erfc


def _pattern_to_int(pattern):
    n = 0
    for bit in pattern:
        n = (n << 1) + int(bit)
    return n


def universal_maurer(bits, patternlen=None, initblocks=None, alpha=0.01, verbose=False):
    # 通用统计检验：用“模式上次出现位置”的平均对数距离来衡量随机性
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    if patternlen is not None:
        l = int(patternlen)
    else:
        ns = [
            904960,
            2068480,
            4654080,
            10342400,
            22753280,
            49643520,
            107560960,
            231669760,
            496435200,
            1059061760,
        ]
        l = 6
        if n < 387840:
            return False, None, {"reason": "长度不够：这个测试至少要 387840 位（现在 %d）" % n}
        for threshold in ns:
            if n >= threshold:
                l += 1

    if l <= 0:
        return False, None, {"reason": "L 必须是正数"}

    nblocks = int(math.floor(float(n) / float(l)))
    if initblocks is not None:
        q = int(initblocks)
    else:
        q = int(10 * (2**l))
    k = nblocks - q
    if k <= 0:
        return False, None, {"reason": "数据太短，分出来的 K 不够用（K=%d）" % k}

    nsymbols = int(2**l)
    # 记录表：保存每个模式最近一次出现的块下标（从 1 开始）
    t = [0 for _ in range(nsymbols)]
    for i in range(q):
        pattern = bits[i * l : (i + 1) * l]
        idx = _pattern_to_int(pattern)
        t[idx] = i + 1

    s = 0.0
    for i in range(q, nblocks):
        pattern = bits[i * l : (i + 1) * l]
        j = _pattern_to_int(pattern)
        dist = (i + 1) - t[j]
        t[j] = i + 1
        s += math.log(dist, 2)

    if verbose:
        print("  sum =", s)

    fn = s / float(k)
    if verbose:
        print("  fn =", fn)

    ev_table = [
        0,
        0.73264948,
        1.5374383,
        2.40160681,
        3.31122472,
        4.25342659,
        5.2177052,
        6.1962507,
        7.1836656,
        8.1764248,
        9.1723243,
        10.170032,
        11.168765,
        12.168070,
        13.167693,
        14.167488,
        15.167379,
    ]
    var_table = [
        0,
        0.690,
        1.338,
        1.901,
        2.358,
        2.705,
        2.954,
        3.125,
        3.238,
        3.311,
        3.356,
        3.384,
        3.401,
        3.410,
        3.416,
        3.419,
        3.421,
    ]

    if l >= len(ev_table):
        return False, None, {"reason": "L 太大，超出参数表（L=%d）" % l}

    denom = (0.7 - 0.8 / float(l) + (4.0 + 32.0 / float(l)) * (pow(float(k), -3.0 / float(l))) / 15.0) * (
        math.sqrt(var_table[l] / float(k))
    ) * math.sqrt(2.0)
    if denom == 0.0:
        return False, None, {"reason": "分母为 0，无法计算"}

    # 用参数表做正态近似，得到统计量，再用互补误差函数得到显著性概率
    mag = abs((fn - ev_table[l]) / denom)
    p = erfc(mag)
    return p >= alpha, p, None

