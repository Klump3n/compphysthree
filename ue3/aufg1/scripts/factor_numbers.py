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

    # brute force loop
    for i in range(20):
        for j in range(20):
            for k in range(20):
                for l in range(20):
                    if (
                            (i + j + k + l == 20)
                            and (k+l <= 10)
                            and (i != j)
                            and (k != l)
                            # and (i+j == 10)
                    ):
                        gridSize.append(i+j)
                        threadSize.append(k+l)
                        gridX.append(i)
                        gridY.append(j)
                        threadX.append(k)
                        threadY.append(l)

    # print(gridX)
    # print(gridY)
    # print(threadX)
    # print(threadY)

    # for it, _ in enumerate(gridX):
    #     print("gridX: {}, gridY: {}, threadX: {}, threadY: {}".format(2**gridX[it], 2**gridY[it], 2**threadX[it], 2**threadY[it]))
    print("gridsize: {}".format(len(gridX)))
    with open('grid_parameters', 'w') as f:
        for it, _ in enumerate(gridX):
            f.write('{}, {}, {}, {},\n'.format(2**gridX[it], 2**gridY[it], 2**threadX[it], 2**threadY[it]))
