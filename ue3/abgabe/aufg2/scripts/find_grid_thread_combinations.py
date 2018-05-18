#!/usr/bin/env python3
"""
Find grid and thread combinations for different grid sizes for the
laplacian thingy...

"""
import numpy as np

def factors(n):
    """
    Stolen from SX.

    """
    n = int(n)
    while n > 1:
        for i in range(2, int(n) + 1):
            if n % i == 0:
                n /= i
                yield int(i)
                break


def factorization(n):
    """
    Return a list with prime factors.

    """
    res = []
    for factor in factors(n):
        res.append(factor)
    return res


def block_thread(factor_list):
    """
    Generate a block/thread configuration from the factorization.

    """
    factor_list = np.asarray(factor_list)
    factor_list = np.sort(factor_list)

    blockX = 1
    blockY = 1
    threadX = 1
    threadY = 1

    if len(factor_list) % 4 == 0:
        threadX = factor_list[2]
        threadY = factor_list[3]
        blockX = factor_list[0]
        blockY = factor_list[1]
        return blockX, blockY, threadX, threadY

    if len(factor_list) == 2:
        threadX = factor_list[0]
        threadY = factor_list[1]
        return blockX, blockY, threadX, threadY

    if len(factor_list) == 3:
        threadX = factor_list[0]
        threadY = factor_list[1]
        blockY = factor_list[2]
        return blockX, blockY, threadX, threadY

    # if all else fails multiply the last two elements of the list
    new_list = []

    for it, number in enumerate(factor_list):
        if it == 0:
            continue
        elif it == 1:
            new_list.append(factor_list[0]*factor_list[1])
            continue
        else:
            new_list.append(number)
            continue

    # call the function recursively
    return block_thread(new_list)

if __name__ == '__main__':

    with open('factorizations', 'w') as fact_file:
        for n in range(1, 64+1):
            N = n + 2
            Nsq = (N)*(N)
            factor_list = factorization(Nsq)
            blockX, blockY, threadX, threadY = block_thread(factor_list)
            fact_file.write('{}, {}, {}, {}, {},\n'.format(n, blockX, blockY, threadX, threadY))
