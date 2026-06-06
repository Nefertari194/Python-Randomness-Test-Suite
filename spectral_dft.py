from __future__ import print_function

import math
import sys

import numpy

from stats_utils import normal_two_sided_p_from_z


def spectral_dft(bits, alpha=0.01, verbose=False):
    # 频谱（离散傅里叶）检验：把 0/1 映射到 -1/+1，做快速傅里叶变换，看幅值超过阈值的数量是否符合期望
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    if (n % 2) == 1:
        bits = bits[:-1]
        n = len(bits)
        if n <= 0:
            return False, None, {"reason": "数据是空的"}

    ts = [(b * 2) - 1 for b in bits]
    # 快速傅里叶变换得到频域系数
    fs = numpy.fft.fft(numpy.array(ts))

    if sys.version_info > (3, 0):
        mags = abs(fs)[: n // 2]
    else:
        mags = abs(fs)[: n / 2]

    t = math.sqrt(math.log(1.0 / 0.05) * float(n))
    # 理论期望：前半段频谱中有 95% 的点低于阈值
    n0 = 0.95 * float(n) / 2.0
    n1 = float(sum(1.0 for mag in mags if mag < t))
    if verbose:
        print("  N0 = %f" % n0)
        print("  N1 = %f" % n1)

    denom = math.sqrt((float(n) * 0.95 * 0.05) / 4.0)
    if denom == 0.0:
        return False, None, {"reason": "分母为 0，无法计算"}
    d = (n1 - n0) / denom
    # 用正态近似把统计量转成双侧显著性概率
    p = normal_two_sided_p_from_z(d)
    return p >= alpha, p, None

