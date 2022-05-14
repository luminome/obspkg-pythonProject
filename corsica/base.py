#!/usr/bin/env python3

#from numpy import genfromtxt
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from shapely import ops
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString, LinearRing

config = {"simplification": 0.005}


def do_refine(shape):
    x, y = shape.coords.xy
    tck, u = interpolate.splprep([x, y], s=0, k=2)
    ost = 0.000001  # 1/shape.length
    u_new = np.arange(0, 1.0 + ost, ost)
    out = interpolate.splev(u_new, tck)
    return out


if __name__ == '__main__':
    my_data = np.genfromtxt('contour20m_Corsica.csv', delimiter=',')
    print(my_data.shape)

    t_ring = LinearRing(my_data)

    sx, sy = my_data[:, 0], my_data[:, 1]
    plt.plot(sx, sy, 'ro', ms=2.0, color=(0.0, 0.0, 0.0, 0.35))
    plt.plot(sx, sy, linewidth=1, color=(0.0, 0.0, 0.0, 0.35))

    # ring_A = t_ring.simplify(config['simplification'])
    # plt.plot(*ring_A.coords.xy, linewidth=1, color=(1.0, 0.0, 1.0))

    ring_B = t_ring.simplify(config['simplification']/10)
    plt.plot(*ring_B.coords.xy, linewidth=1, color=(1.0, 0.0, 0.0))

    p_out = do_refine(ring_B)
    p_coords = zip(p_out[0], p_out[1])
    plt.plot(p_out[0], p_out[1], linewidth=1, color=(1.0, 0.5, 0.0))

    # p_last = Polygon(p_coords).simplify(config['simplification']/10)
    # sx, sy = p_last.exterior.xy
    # plt.plot(sx, sy, 'ro', ms=3.0, color=(0.0, 0.0, 0.0, 0.8))
    # plt.plot(*p_last.exterior.xy, linewidth=1, color=(0.0, 0.0, 0.0))

    # offsetted_right = ring_A.parallel_offset(0.005, 'right', resolution=4, join_style=2, mitre_limit=0.1)
    # if offsetted_right.type == 'MultiLineString':
    #     for li in offsetted_right.geoms:
    #         offsetted = li.parallel_offset(0.005, 'right', resolution=4, join_style=2, mitre_limit=0.1)
    #         if offsetted.type == 'MultiLineString':
    #             for lp in offsetted.geoms:
    #                 plt.plot(*lp.coords.xy, linewidth=1, color=(0.0, 0.0, 0.0))
    #         else:
    #             plt.plot(*offsetted.coords.xy, linewidth=1, color=(0.0, 0.0, 0.0))
    # else:
    #     print('dolo')
    #     plt.plot(*offsetted_right.coords.xy, linewidth=1, color=(0.0, 0.0, 0.0))




    plt.axis('scaled')
    plt.show()
