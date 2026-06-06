from __future__ import print_function


def count_ones_zeroes(bits):
    # 统计 0/1 个数：多个测试都会用到
    ones = 0
    zeroes = 0
    for bit in bits:
        if bit == 1:
            ones += 1
        else:
            zeroes += 1
    return zeroes, ones


def bits_to_int(bits):
    # 把一个比特序列（高位在前）转成整数，用于模式计数类测试
    x = 0
    for b in bits:
        x = (x << 1) + int(b)
    return x

