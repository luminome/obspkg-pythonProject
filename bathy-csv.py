#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt


if __name__ == '__main__':
    dat = np.loadtxt(open("./data/bathy_med.csv", "rb"), delimiter=";", skiprows=1)
    dy = np.loadtxt(open("./data/lat_vector.csv", "rb"), delimiter=";", skiprows=1)

    h, w = dat.shape
    print(w, h)
    dx = np.linspace(-7.0, 42.5, w)
    print(dy[0], dy[h - 1])
    print(dx[0], dx[w - 1])

    minx = dx[0]
    miny = dy[h - 1]
    maxx = dx[w - 1]
    maxy = dy[0]
    #extent(minx,maxx,miny,maxy)
    # extent floats (left, right, bottom, top)
    # g2 = np.mgrid[dx[0]:dx[w - 1]:w, dy[0]:dy[h - 1]:h]
    # print(g2.shape)

    gdx, gdy = np.meshgrid(dx, dy, indexing='ij')
    coordinate_grid = np.array([gdx, gdy])
    # print(coordinate_grid.shape)
    # print(gdx)
    # print(gdy)
    plt.imshow(dat, extent=(minx, maxx, miny, maxy), interpolation="none")
    plt.scatter(gdx, gdy, s=0.01)
    plt.show()


    exit()


    print(dat.shape)
    print(dy.shape)

    dx = np.linspace(0, 10, dat.shape[1])
    #print(dx)
    print(dx.shape, dy.shape)

    gdx, gdy = np.meshgrid(dx, dy, indexing='ij')
    coordinate_grid = np.array([gdx, gdy])
    print(coordinate_grid.shape)
    print(gdx)
    print(gdy)

    # lex = zip(dx, dy)
    # print(lex)
    # g = np.arange(80, 122, 2)

    # pos = np.array([dx, dy], dtype=object)
    #
    # print(pos.shape)
    # bathy_med.csv
    # lat_vector.csv
    # lat_vector.csv

    def f(x, y):
        return (1 - x / 2 + x ** 5 + y ** 3) * np.exp(-x ** 2 - y ** 2)


    n = 10
    x = np.linspace(-3, 3, 4 * n)
    y = np.linspace(-3, 3, 3 * n)
    X, Y = np.meshgrid(dx, dy)
    plt.imshow(dat, interpolation="none")  #;//f(X, Y))

    plt.xticks(dx, 10)

    plt.show()
    #
    # x = np.linspace(0, dat.shape[1], 1)
    # y = 2 * x + 5
    # plt.title("Matplotlib demo")
    # plt.xlabel("x axis caption")
    # plt.ylabel("y axis caption")
    # plt.plot(dat[:, 0], dat[:, 1])  #x, y)
    # plt.show()
