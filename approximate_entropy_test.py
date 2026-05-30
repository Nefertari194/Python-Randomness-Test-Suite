from __future__ import print_function

import math
from gamma_functions import *

def bits_to_int(bits):
    theint = 0
    for i in range(len(bits)):
        theint = (theint << 1) + bits[i]
    return theint
        
def approximate_entropy_test(bits, alpha=0.01, verbose=False):
    n = len(bits)
    
    m = int(math.floor(math.log(n,2)))-6
    if m < 2:
        m = 2
    if m >3 :
        m = 3
        
    if verbose:
        print("  n         = ", n)
        print("  m         = ", m)
    
    Cmi = list()
    phi_m = list()
    for iterm in range(m,m+2):
        # 第一步
        padded_bits=bits+bits[0:iterm-1]
    
        # 第二步
        counts = list()
        for i in range(2**iterm):
            count = 0
            for j in range(n):
                if bits_to_int(padded_bits[j:j+iterm]) == i:
                    count += 1
            counts.append(count)
            if verbose:
                print("  Pattern %d of %d, count = %d" % (i+1,2**iterm, count))
    
        # 第三步
        Ci = list()
        for i in range(2**iterm):
            Ci.append(float(counts[i])/float(n))
        
        Cmi.append(Ci)
    
        # 第四步
        sum = 0.0
        for i in range(2**iterm):
            if (Ci[i] > 0.0):
                sum += Ci[i]*math.log(Ci[i])
        phi_m.append(sum)
        if verbose:
            print("  phi(%d)    = %f" % (iterm, sum))
        
    # 第五步：上面的循环做完就行
    
    # 第六步
    appen_m = phi_m[0] - phi_m[1]
    if verbose:
        print("  AppEn(%d)  = %f" % (m, appen_m))
    chisq = 2*n*(math.log(2) - appen_m)
    if verbose:
        print("  ChiSquare = ", chisq)
    # 第七步
    p = gammaincc(2**(m-1),(chisq/2.0))
    
    success = (p >= alpha)
    return (success, p, None)

if __name__ == "__main__":
    bits = [1,1,0,0,1,0,0,1,0,0,0,0,1,1,1,1,
            1,1,0,1,1,0,1,0,1,0,1,0,0,0,1,0,
            0,0,1,0,0,0,0,1,0,1,1,0,1,0,0,0,
            1,1,0,0,0,0,1,0,0,0,1,1,0,1,0,0,
            1,1,0,0,0,1,0,0,1,1,0,0,0,1,1,0,
            0,1,1,0,0,0,1,0,1,0,0,0,1,0,1,1,
            1,0,0,0]
    success, p, _ = approximate_entropy_test(bits)
    
    print("success =",success)
    print("p = ",p)
    
