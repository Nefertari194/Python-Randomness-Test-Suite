

from __future__ import print_function

import math
from gamma_functions import *

# 随机游走测试
def random_excursion_test(bits, alpha=0.01, verbose=False, min_J=500):
    n = len(bits)

    x = list()             # 把 0/1 变成 -1/+1
    for bit in bits:
        x.append((bit*2)-1)

    # 算累加和
    pos = 0
    s = list()
    for e in x:
        pos = pos+e
        s.append(pos)    
    sprime = [0]+s+[0] # 两头补 0
    
    # 把序列拆成一段段的循环段
    pos = 1
    cycles = list()
    while (pos < len(sprime)):
        cycle = list()
        cycle.append(0)
        while sprime[pos]!=0:
            cycle.append(sprime[pos])
            pos += 1
        cycle.append(0)
        cycles.append(cycle)
        pos = pos + 1
    
    J = len(cycles)
    if verbose:
        print("J="+str(J))    
    if J <= 0:
        return False, None, {"reason": "循环数 J 太小（%d）" % J}
    if (J < min_J):
        return False, None, {"reason": "循环数 J 太小（%d < %d），这个测试结果不靠谱" % (J, min_J)}
    
    vxk = [['a','b','c','d','e','f'] for y in [-4,-3,-2,-1,1,2,3,4] ]

    # 统计每种状态出现次数
    for k in range(6):
        for index in range(8):
            mapping = [-4,-3,-2,-1,1,2,3,4]
            x = mapping[index]
            cyclecount = 0
            # 统计有多少个循环段里，x 出现了 k 次
            for cycle in cycles:
                oc = 0
                # 数一下这个循环段里 x 出现了几次
                for pos in cycle:
                    if (pos == x):
                        oc += 1
                # 满足条件就加 1
                if (k < 5):
                    if oc == k:
                        cyclecount += 1
                else:
                    if k == 5:
                        if oc >=5:
                            cyclecount += 1
            vxk[index][k] = cyclecount
    
    # 参考概率表
    pixk=[[0.5     ,0.25   ,0.125  ,0.0625  ,0.0312 ,0.0312],
          [0.75    ,0.0625 ,0.0469 ,0.0352  ,0.0264 ,0.0791],
          [0.8333  ,0.0278 ,0.0231 ,0.0193  ,0.0161 ,0.0804],
          [0.875   ,0.0156 ,0.0137 ,0.012   ,0.0105 ,0.0733],
          [0.9     ,0.01   ,0.009  ,0.0081  ,0.0073 ,0.0656],
          [0.9167  ,0.0069 ,0.0064 ,0.0058  ,0.0053 ,0.0588],
          [0.9286  ,0.0051 ,0.0047 ,0.0044  ,0.0041 ,0.0531]]
    
    success = True
    plist = list()
    for index in range(8):
        mapping = [-4,-3,-2,-1,1,2,3,4]
        x = mapping[index]
        chisq = 0.0
        for k in range(6):
            top = float(vxk[index][k]) - (float(J) * (pixk[abs(x)-1][k]))
            top = top*top
            bottom = J * pixk[abs(x)-1][k]
            chisq += top/bottom
        p = gammaincc(5.0/2.0,chisq/2.0)
        plist.append(p)
        if p < alpha:
            err = " Not Random"
            success = False
        else:
            err = ""
        if verbose:
            print("x = %1.0f\tchisq = %f\tp = %f %s"  % (x,chisq,p,err))

    return (success, None, plist)

if __name__ == "__main__":
    bits = [0,1,1,0,1,1,0,1,0,1]
    success, _, plist = random_excursion_test(bits)
    
    print("success =",success)
    print("plist = ",plist)
