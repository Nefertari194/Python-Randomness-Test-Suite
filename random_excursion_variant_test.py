

from __future__ import print_function

import math

# 随机游走变体测试
def random_excursion_variant_test(bits, alpha=0.01, verbose=False, min_J=500):
    n = len(bits)

    x = list()             # 把 0/1 变成 -1/+1
    for bit in bits:
        x.append((bit * 2)-1)

    # 算累加和
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)    
    sprime = [0]+s+[0] # 两头补 0

    # 统计循环数
    J = 0
    for value in sprime[1:]:
        if value == 0:
            J += 1
    if verbose:
        print("J=",J)
    if J <= 0:
        return False, None, {"reason": "循环数 J 太小（%d）" % J}
    # 统计每个状态出现次数
    count = [0 for _ in range(19)]
    for value in sprime:
        if (abs(value) < 10):
            count[value + 9] += 1

    # 算 p 值
    success = True
    plist = list()
    for x in range(-9,10):
        if x != 0:
            top = abs(count[x + 9]-J)
            bottom = math.sqrt(2.0 * J *((4.0*abs(x))-2.0))
            p = math.erfc(top/bottom)
            plist.append(p)
            if p < alpha:
                err = " Not Random"
                success = False
            else:
                err = ""
            if verbose:
                print("x = %1.0f\t count=%d\tp = %f %s"  % (x, count[x + 9], p, err))
            
    if (J < min_J):
        return False, None, {"reason": "循环数 J 太小（%d < %d），这个测试结果不靠谱" % (J, min_J)}
    return (success,None,plist)
