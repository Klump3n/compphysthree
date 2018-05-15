#!/usr/bin/env python3
"""
Factor large numbers (a power of two).

"""
if __name__ == '__main__':
    number = 1024*1024          # 2^(20)

    exp = 20
    exp_max = 10

    gridX = []
    gridY = []
    threadX = []
    threadY = []

    gridSize = []
    threadSize = []

    for i in range(20):
        for j in range(20):
            for k in range(20):
                for l in range(20):
                    if (
                            (i + j + k + l == 20) and
                            (k <= 10) and
                            (l <= 10) and
                            (i != j) and
                            (k != l) and
                            (i+j not in gridSize) and
                            (k+l not in threadSize)
                    ):
                        gridSize.append(i+j)
                        threadSize.append(k+l)
                        gridX.append(i)
                        gridY.append(j)
                        threadX.append(k)
                        threadY.append(l)

    print(gridX)
    print(gridY)
    print(threadX)
    print(threadY)
