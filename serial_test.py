

from __future__ import print_function

import math
from gamma_functions import gammaincc

import os

from bitio import iter_bit_sequences

def int2patt(n,m):
    pattern = list()
    for i in range(m):
        pattern.append((n >> i) & 1)
    return pattern
    
def countpattern(patt,bits,n):
    thecount = 0
    for i in range(n):
        match = True
        for j in range(len(patt)):
            if patt[j] != bits[i+j]:
                match = False
        if match:
            thecount += 1
    return thecount

def psi_sq_mv1(m, n, padded_bits):
    psi_sq_m = 0.0
    for i in range(2**m):
        pattern = int2patt(i,m)
        count = countpattern(pattern, padded_bits, n)
        psi_sq_m += (count**2)
    psi_sq_m = psi_sq_m * (2**m)/n 
    psi_sq_m -= n
    return psi_sq_m            
         
def _pick_default_m(n, max_m=16):
    if n <= 0:
        return None
    m = int(math.floor(math.log(n, 2))) - 2
    if m < 2:
        return None
    return min(max_m, m)


def _check_serial_prereq(n, m, max_m=16):
    if n <= 0:
        return False, "数据是空的"
    if m is None:
        return False, "算不出 m（数据太短）"
    if m < 2:
        return False, "m 不能小于 2"
    if m > max_m:
        return False, "m 太大了（默认最多 %d），这个算起来会很慢" % max_m
    if n < (2 ** (m + 2)):
        return False, "长度不够：需要 n >= 2^(m+2) 才比较靠谱（现在 n=%d, m=%d）" % (n, m)
    return True, None


def serial_test(bits, patternlen=None, alpha=0.01, verbose=False, max_m=16):
    n = len(bits)
    if patternlen is not None:
        m = int(patternlen)
    else:
        m = _pick_default_m(n, max_m=max_m)

    ok, reason = _check_serial_prereq(n, m, max_m=max_m)
    if not ok:
        if verbose:
            print("  条件不足：%s" % reason)
        return False, None, {"reason": reason, "p_values": None, "m": m, "n": n}

    # 先把前面 m-1 个比特拼到后面，方便循环统计
    padded_bits = bits + bits[0:m-1]

    psi_sq_m   = psi_sq_mv1(m, n, padded_bits)
    psi_sq_mm1 = psi_sq_mv1(m-1, n, padded_bits)
    psi_sq_mm2 = psi_sq_mv1(m-2, n, padded_bits)    
    
    delta1 = psi_sq_m - psi_sq_mm1
    delta2 = psi_sq_m - (2*psi_sq_mm1) + psi_sq_mm2
    
    P1 = gammaincc(2**(m-2),delta1/2.0)
    P2 = gammaincc(2**(m-3),delta2/2.0)

    if verbose:
        print("  psi_sq_m   = ", psi_sq_m)
        print("  psi_sq_mm1 = ", psi_sq_mm1)
        print("  psi_sq_mm2 = ", psi_sq_mm2)
        print("  delta1     = ", delta1)
        print("  delta2     = ", delta2)
        print("  P1         = ", P1)
        print("  P2         = ", P2)

    success = (P1 >= alpha) and (P2 >= alpha)
    return success, None, {"p_values": [P1, P2], "m": m, "n": n, "alpha": alpha}


def run_serial_test(config):
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

        success, _, detail = serial_test(bits, patternlen=m, alpha=alpha, verbose=verbose, max_m=max_m)
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

if __name__ == "__main__":
    # 这里就简单跑个小例子，主要看一下函数能不能正常跑通
    bits = [0, 0, 1, 1, 0, 1, 1, 1] * 16
    success, _, detail = serial_test(bits, patternlen=3, alpha=0.01, verbose=True)
    if detail["p_values"] is None:
        print("结果：skip  %s" % detail.get("reason"))
    else:
        p1, p2 = detail["p_values"]
        print("结果：%s  P1=%.6g  P2=%.6g" % ("pass" if success else "fail", p1, p2))
    
