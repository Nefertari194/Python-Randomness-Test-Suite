from __future__ import print_function

import math
import os

from bitio import iter_bit_sequences
from stats_utils import chi2_sf


def _int_to_pattern(n, m):
    patt = []
    for i in range(m):
        patt.append((n >> i) & 1)
    return patt


def _count_pattern(patt, bits, n):
    thecount = 0
    m = len(patt)
    for i in range(n):
        match = True
        for j in range(m):
            if patt[j] != bits[i + j]:
                match = False
        if match:
            thecount += 1
    return thecount


def _psi_sq(m, n, padded_bits):
    # 计算 ψ² 统计量（模式计数的平方和）
    psi_sq_m = 0.0
    for i in range(2**m):
        pattern = _int_to_pattern(i, m)
        count = _count_pattern(pattern, padded_bits, n)
        psi_sq_m += float(count**2)
    psi_sq_m = psi_sq_m * float(2**m) / float(n)
    return psi_sq_m - float(n)


def _pick_default_m(n, max_m=16):
    if n <= 0:
        return None
    m = int(math.floor(math.log(n, 2))) - 2
    if m < 2:
        return None
    return min(int(max_m), int(m))


def _check_prereq(n, m, max_m=16):
    if n <= 0:
        return False, "数据是空的"
    if m is None:
        return False, "算不出 m（数据太短）"
    if m < 2:
        return False, "m 不能小于 2"
    if m > int(max_m):
        return False, "m 太大了（默认最多 %d），这个算起来会很慢" % int(max_m)
    if n < (2 ** (m + 2)):
        return False, "长度不够：需要 n >= 2^(m+2) 才比较靠谱（现在 n=%d, m=%d）" % (n, m)
    return True, None


def serial_patterns(bits, patternlen=None, alpha=0.01, verbose=False, max_m=16):
    # 序列检验：统计三种模式长度下的 ψ²，并构造两个卡方统计量
    n = len(bits)
    if patternlen is not None:
        m = int(patternlen)
    else:
        m = _pick_default_m(n, max_m=max_m)

    ok, reason = _check_prereq(n, m, max_m=max_m)
    if not ok:
        if verbose:
            print("  条件不足：%s" % reason)
        return False, None, {"reason": reason, "p_values": None, "m": m, "n": n}

    # 循环统计：把前(模式长度-1)个比特拼到后面
    padded_bits = bits + bits[0 : m - 1]

    psi_sq_m = _psi_sq(m, n, padded_bits)
    psi_sq_mm1 = _psi_sq(m - 1, n, padded_bits)
    psi_sq_mm2 = _psi_sq(m - 2, n, padded_bits)

    delta1 = psi_sq_m - psi_sq_mm1
    delta2 = psi_sq_m - (2.0 * psi_sq_mm1) + psi_sq_mm2

    # 两个统计量分别对应不同自由度的卡方分布
    p1 = chi2_sf(delta1, 2 ** (m - 1))
    p2 = chi2_sf(delta2, 2 ** (m - 2))

    if verbose:
        print("  psi_sq_m   = ", psi_sq_m)
        print("  psi_sq_mm1 = ", psi_sq_mm1)
        print("  psi_sq_mm2 = ", psi_sq_mm2)
        print("  delta1     = ", delta1)
        print("  delta2     = ", delta2)
        print("  P1         = ", p1)
        print("  P2         = ", p2)

    success = (p1 >= alpha) and (p2 >= alpha)
    return success, None, {"p_values": [p1, p2], "m": m, "n": n, "alpha": alpha}


def run_serial_patterns(config):
    key_path = config.get("key_path")
    if not key_path:
        raise ValueError("config 里得有 key_path")

    m = config.get("m", None)
    alpha = float(config.get("alpha", 0.01))
    input_format = config.get("input_format", "auto")
    bigendian = bool(config.get("bigendian", True))
    recursive = bool(config.get("recursive", False))
    max_m = int(config.get("max_m", 16))
    verbose = bool(config.get("verbose", False))
    print_results = bool(config.get("print_results", True))

    results = []
    for item in iter_bit_sequences(key_path, input_format=input_format, bigendian=bigendian, recursive=recursive):
        bits = item["bits"]
        source = item["source"]
        index = item["index"]

        success, _, detail = serial_patterns(bits, patternlen=m, alpha=alpha, verbose=verbose, max_m=max_m)
        p_values = detail.get("p_values") if isinstance(detail, dict) else None
        reason = detail.get("reason") if isinstance(detail, dict) else None

        if p_values is None:
            status = "skip"
        else:
            status = "pass" if success else "fail"

        row = {
            "source": source,
            "index": index,
            "status": status,
            "p_values": p_values,
            "reason": reason,
            "n": detail.get("n") if isinstance(detail, dict) else len(bits),
            "m": detail.get("m") if isinstance(detail, dict) else m,
        }
        results.append(row)

        if print_results:
            tag = "%s#%d" % (os.path.basename(source), index)
            if status == "skip":
                print("%s  skip  %s" % (tag, reason))
            else:
                p1, p2 = p_values
                print("%s  %s  P1=%.6g  P2=%.6g" % (tag, status, p1, p2))

    return results

