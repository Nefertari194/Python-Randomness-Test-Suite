from __future__ import print_function

import math

from stats_utils import chi2_sf


def _pr(u, eta):
    if u == 0:
        return math.exp(-eta)

    s = 0.0
    for l in range(1, u + 1):
        term = -eta - (u * math.log(2.0)) + (l * math.log(eta)) - math.lgamma(l + 1) + math.lgamma(u) - math.lgamma(l) - math.lgamma(u - l + 1)
        s += math.exp(term)
    return s


def template_overlapping(bits, blen=6, alpha=0.01, verbose=False, block_size=None, num_blocks=None, allow_degraded=False):
    # 重叠模板匹配：统计每块里（允许重叠的）模板命中次数的分布，再做卡方检验
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    m = int(blen)
    if m < 2:
        m = 2
    b = [1 for _ in range(m)]

    k = 5
    m_std = 1062
    n_std = 968

    if allow_degraded:
        if block_size is not None:
            m_blk = int(block_size)
            if m_blk <= 0:
                return False, None, {"reason": "块大小得是正数"}
            n_blk = int(math.floor(float(n) / float(m_blk)))
        elif num_blocks is not None:
            n_blk = int(num_blocks)
            if n_blk <= 0:
                return False, None, {"reason": "块数得是正数"}
            m_blk = int(math.floor(float(n) / float(n_blk)))
        else:
            n_blk = max(1, int(math.floor(float(n) / 128.0)))
            m_blk = int(math.floor(float(n) / float(n_blk)))

        if n_blk <= 0 or m_blk <= 0:
            return False, None, {"reason": "分块参数不合适"}
        if m_blk <= m:
            return False, None, {"reason": "块太短了，模板放不下"}
        if n_blk < 2:
            return False, None, {"reason": "块数太少，结果不太靠谱"}
        m_use = m_blk
        n_use = n_blk
    else:
        if n < (m_std * n_std):
            return False, None, {"reason": "长度不够：至少要 %d 位（现在 %d）" % (m_std * n_std, n)}
        m_use = m_std
        n_use = n_std

    blocks = []
    for i in range(n_use):
        blocks.append(bits[i * m_use : (i + 1) * m_use])

    v = [0 for _ in range(k + 1)]
    for block in blocks:
        count = 0
        for position in range(m_use - m + 1):
            if block[position : position + m] == b:
                count += 1

        if count >= k:
            v[k] += 1
        else:
            v[count] += 1

    lambd = (float(m_use - m + 1)) / float(2**m)
    eta = lambd / 2.0

    # 命中次数的理论分布（按标准的近似公式计算）
    pi = [0.0 for _ in range(k + 1)]
    s = 0.0
    for i in range(k):
        pi[i] = _pr(i, eta)
        s += pi[i]
    pi[k] = 1.0 - s

    chi_sq = 0.0
    for i in range(k + 1):
        chi_sq += ((v[i] - (n_use * pi[i])) ** 2) / (n_use * pi[i])

    p = chi2_sf(chi_sq, 5)
    if verbose:
        print("  B = ", b)
        print("  m = ", m)
        print("  M = ", m_use)
        print("  N = ", n_use)
        print("  K = ", k)
        print("  v[j] =  ", v)
        print("  chisq = ", chi_sq)

    return p >= alpha, p, None

