
 
from __future__ import print_function

import math
from gamma_functions import *

def berelekamp_massey(bits):
    n = len(bits)
    b = [0 for x in bits]  # 初始化两个数组
    c = [0 for x in bits]
    b[0] = 1
    c[0] = 1
    
    L = 0
    m = -1
    N = 0
    while (N < n):
        # 算偏差
        d = bits[N]
        for i in range(1,L+1):
            d = d ^ (c[i] & bits[N-i])
        if (d != 0):  # d 不为 0 就调整多项式
            t = c[:]
            for i in range(0,n-N+m):
                c[N-m+i] = c[N-m+i] ^ b[i] 
            if (L <= (N/2)):
                L = N + 1 - L
                m = N
                b = t 
        N = N +1
    # 返回复杂度 L 和多项式
    return L , c[0:L]
    
def linear_complexity_test(bits, patternlen=None, alpha=0.01, verbose=False):
    n = len(bits)
    # 第一步：选块大小
    if patternlen is not None:
        M = patternlen  
    else: 
        if n < 1000000:
            return False, None, {"reason": "长度不够：这个测试默认要至少 10^6 位"}
        M = 512
    K = 6 
    N = int(math.floor(n/M))
    if verbose:
        print("  M = ", M)
        print("  N = ", N)
        print("  K = ", K)    
    
    # 第二步：算每个块的线性复杂度
    LC = list()
    for i in range(N):
        x = bits[(i*M):((i+1)*M)]
        LC.append(berelekamp_massey(x)[0])
    
    # 第三步：算均值
    a = float(M)/2.0
    b = ((((-1)**(M+1))+9.0))/36.0
    c = ((M/3.0) + (2.0/9.0))/(2**M)
    mu =  a+b-c
    
    T = list()
    for i in range(N):
        x = ((-1.0)**M) * (LC[i] - mu) + (2.0/9.0)
        T.append(x)
        
    # 第四步：统计落在哪个区间
    v = [0,0,0,0,0,0,0]
    for t in T:
        if t <= -2.5:
            v[0] += 1
        elif t <= -1.5:
            v[1] += 1
        elif t <= -0.5:
            v[2] += 1
        elif t <= 0.5:
            v[3] += 1
        elif t <= 1.5:
            v[4] += 1
        elif t <= 2.5:
            v[5] += 1            
        else:
            v[6] += 1

    # 第五步：算卡方
    pi = [0.010417,0.03125,0.125,0.5,0.25,0.0625,0.020833]
    chisq = 0.0
    for i in range(K+1):
        chisq += ((v[i] - (N*pi[i]))**2.0)/(N*pi[i])
    if verbose:
        print("  chisq = ", chisq)
    # 第六步：算 p 值
    P = gammaincc((K/2.0),(chisq/2.0))
    if verbose:
        print("  P = ", P)
    success = (P >= alpha)
    return (success, P, None)
    
if __name__ == "__main__":
    bits = [1,1,0,1,0,1,1,1,1,0,0,0,1]
    L,poly = berelekamp_massey(bits)

    bits = [1,1,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,
            0,1,1,1,0,1,0,1,1,1,1,0,0,0,1,1,1,0,1,0,1,1,1,1,
            0,0,0,1,1,1,0,1,0,1,1,1,1,0,0,0,1]
    success,p,_ = linear_complexity_test(bits,patternlen=7)
    
    print("L =",L)
    print("p = ",p)
       
