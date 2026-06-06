from __future__ import print_function

from stats_utils import chi2_sf


def random_excursions(bits, alpha=0.01, verbose=False, min_J=500):
    # 随机游走检验：把序列映射成随机游走，按“循环段”统计状态访问次数分布，再做卡方检验
    n = len(bits)
    if n <= 0:
        return False, None, {"reason": "数据是空的"}

    # 0/1 -> -1/+1，形成游走步长
    x = [(bit * 2) - 1 for bit in bits]

    pos = 0
    s = []
    for e in x:
        pos += e
        s.append(pos)
    sprime = [0] + s + [0]

    # 按回到 0 的位置把游走拆成多个循环段
    pos = 1
    cycles = []
    while pos < len(sprime):
        cycle = [0]
        while sprime[pos] != 0:
            cycle.append(sprime[pos])
            pos += 1
        cycle.append(0)
        cycles.append(cycle)
        pos = pos + 1

    j = len(cycles)
    if verbose:
        print("J=" + str(j))
    if j <= 0:
        return False, None, {"reason": "循环数 J 太小（%d）" % j}
    if j < int(min_J):
        return False, None, {"reason": "循环数 J 太小（%d < %d），这个测试结果不靠谱" % (j, int(min_J))}

    vxk = [[0, 0, 0, 0, 0, 0] for _ in [-4, -3, -2, -1, 1, 2, 3, 4]]
    mapping = [-4, -3, -2, -1, 1, 2, 3, 4]

    for k in range(6):
        for index in range(8):
            state = mapping[index]
            cyclecount = 0
            for cycle in cycles:
                oc = 0
                for p in cycle:
                    if p == state:
                        oc += 1
                if k < 5:
                    if oc == k:
                        cyclecount += 1
                else:
                    if oc >= 5:
                        cyclecount += 1
            vxk[index][k] = cyclecount

    pixk = [
        [0.5, 0.25, 0.125, 0.0625, 0.0312, 0.0312],
        [0.75, 0.0625, 0.0469, 0.0352, 0.0264, 0.0791],
        [0.8333, 0.0278, 0.0231, 0.0193, 0.0161, 0.0804],
        [0.875, 0.0156, 0.0137, 0.012, 0.0105, 0.0733],
        [0.9, 0.01, 0.009, 0.0081, 0.0073, 0.0656],
        [0.9167, 0.0069, 0.0064, 0.0058, 0.0053, 0.0588],
        [0.9286, 0.0051, 0.0047, 0.0044, 0.0041, 0.0531],
    ]

    success = True
    plist = []
    for index in range(8):
        state = mapping[index]
        chi_sq = 0.0
        for k in range(6):
            expected = float(j) * float(pixk[abs(state) - 1][k])
            top = float(vxk[index][k]) - expected
            chi_sq += (top * top) / expected
        p = chi2_sf(chi_sq, 5)
        plist.append(p)
        if p < alpha:
            success = False
        if verbose:
            err = " Not Random" if p < alpha else ""
            print("x = %1.0f\tchisq = %f\tp = %f %s" % (state, chi_sq, p, err))

    return success, None, plist

