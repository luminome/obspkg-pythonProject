from skimage import measure
from skimage.draw import bezier_curve
from scipy import interpolate
from scipy.interpolate.fitpack import insert
from scipy.ndimage.filters import gaussian_filter

from scipy.interpolate import interp1d
import scipy.optimize as spopt

import matplotlib.pyplot as plt
from shapely import ops
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString
from shapely.validation import make_valid
from shapely.geometry import CAP_STYLE, JOIN_STYLE
import json
import math
import pathlib
import numpy as np


# #//STEP ONE: CONVERT CONTOUR TO ARRAY OF BEZIERS.
# #//STEP TWO: FIND INTERSECTION OF SECTOR AND PARAMETRIZED(BOUNDED) CURVE.
# #//step THREE: SLICE CURVE ON SECTOR BOUNDS USING (X,Y) TO FIND POINT AND T TO DO THE REST.




# https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splprep.html
# https://jike.in/?qa=424783/python-approximating-data-with-a-multi-segment-cubic-bezier-curve-and-a-distance-as-well-as-a-curvature-contraint
# https://stackoverflow.com/questions/29934831/matplotlib-draw-spline-from-3-points
# https://pomax.github.io/bezierinfo/#tracing  // AND  previoustable of contentsnextFinding Y, given X

#
# function controlPoints(p) {
#   // given the points array p calculate the control points
#   var pc = [];
#   for (var i = 1; i < p.length - 1; i++) {
#     var dx = p[i - 1].x - p[i + 1].x; // difference x
#     var dy = p[i - 1].y - p[i + 1].y; // difference y
#     // the first control point
#     var x1 = p[i].x - dx * t;
#     var y1 = p[i].y - dy * t;
#     var o1 = {
#       x: x1,
#       y: y1
#     };
#
#     // the second control point
#     var x2 = p[i].x + dx * t;
#     var y2 = p[i].y + dy * t;
#     var o2 = {
#       x: x2,
#       y: y2
#     };
#
#     // building the control points array
#     pc[i] = [];
#     pc[i].push(o1);
#     pc[i].push(o2);
#   }
#   return pc;
# }
#


def ctrl_points(co):
    # https://codepen.io/enxaneta/pen/PqLNLv?editors=0010
    t = 1 / 5

    for i, pt in enumerate(co):
        try:
            dx = co[i - 1][0] - co[i + 1][0]  # #// difference x
            dy = co[i - 1][1] - co[i + 1][1]  # #// difference y

            # #// the first control point
            x1 = co[i][0] - dx * t
            y1 = co[i][1] - dy * t
            o1 = {
              "x": x1,
              "y": y1
            }

            # #// the second control point
            x2 = co[i][0] + dx * t
            y2 = co[i][1] + dy * t
            o2 = {
              "x": x2,
              "y": y2
            }

        except IndexError:
            print(f'no index at{i}')
            pass

        print(pt)


def b_spline_to_bezier_series(tck, per=False):
    """Convert a parametric b-spline into a sequence of Bezier curves of the same degree.

    Inputs:
    tck : (t,c,k) tuple of b-spline knots, coefficients, and degree returned by splprep.
    per : if tck was created as a periodic spline, per *must* be true, else per *must* be false.

    Output:
    A list of Bezier curves of degree k that is equivalent to the input spline.
    Each Bezier curve is an array of shape (k+1,d) where d is the dimension of the
    space; thus the curve includes the starting point, the k-1 internal control
    points, and the endpoint, where each point is of d dimensions.
    """
    # from fitpack import insert
    # from numpy import asarray, unique, split, sum
    t, c, k = tck
    t = np.asarray(t)
    try:
        c[0][0]
    except IndexError:
        # I can't figure out a simple way to convert nonparametric splines to
        # parametric splines. Oh well.
        raise TypeError("Only parametric b-splines are supported.")

    new_tck = tck
    if per:
        # ignore the leading and trailing k knots that exist to enforce periodicity
        knots_to_consider = np.unique(t[k:-k])
    else:
        # the first and last k+1 knots are identical in the non-periodic case, so
        # no need to consider them when increasing the knot multiplicities below
        knots_to_consider = np.unique(t[k+1:-k-1])
        # For each unique knot, bring it's multiplicity up to the next multiple of k+1
        # This removes all continuity constraints between each of the original knots,
        # creating a set of independent Bezier curves.

    desired_multiplicity = k+1
    for x in knots_to_consider:
        current_multiplicity = np.sum(t == x)
        remainder = current_multiplicity % desired_multiplicity
        if remainder != 0:
            # add enough knots to bring the current multiplicity up to the desired multiplicity
            number_to_insert = desired_multiplicity - remainder
            new_tck = insert(x, new_tck, number_to_insert, per)

    tt,cc,kk = new_tck
    # strip off the last k+1 knots, as they are redundant after knot insertion
    bezier_points = np.transpose(cc)[:-desired_multiplicity]
    if per:
        # again, ignore the leading and trailing k knots
        bezier_points = bezier_points[k:-k]
    # group the points into the desired bezier curves

    return np.split(bezier_points, len(bezier_points) / desired_multiplicity, axis = 0)
    # return bezier_points

if __name__ == '__main__':

    src_depths = {
        "data": np.loadtxt(open("../data/bathy_med_tt.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "lons": np.loadtxt(open("../data/lon_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "lats": np.loadtxt(open("../data/lat_vector_t.csv", "rb"), delimiter=";", encoding=None, skiprows=0),
        "reso": [None, 10, 6, 2, 1],
        "degs": 60
    }

    dat = src_depths['data']
    # x = src_depths['lons']
    #xt = np.arange(x.min(), x.max(), 1.0)
    xt = np.arange(-6, 36.5, 1.0)
    yt = np.arange(30, 46, 1.0)
    # x = np.array([0, 1, 2, 3])
    # y = np.array([0.650, 0.660, 0.675, 0.685])
    # my_xticks = ['a', 'b', 'c', 'd']
    # plt.xticks(x, my_xticks)
    # plt.yticks(np.arange(y.min(), y.max(), 0.005))

    #, extent = (-6, 36.5, 30, 46),
    # plt.xticks(xt)
    # plt.yticks(yt)
    #
    # plt.xlim(-6, 36.5)
    # plt.ylim(30, 46)

    plt.imshow(dat, interpolation="none")
    #
    # m = 60
    # xx = 8
    # yy = 4
    #
    # sk = box(xx*m, yy*m, (xx+1)*m, (yy+1)*m)
    # tk = sk.buffer(m*0.2, join_style=2, mitre_limit=5.0)

    def graze(boxtuple, dat=False):
        m = 60
        xx, yy, w = boxtuple
        inside = box(xx * m, yy * m, (xx + w) * m, (yy + w) * m)
        outside = inside.buffer(m * 0.5, join_style=2, mitre_limit=5.0)
        return [inside, outside]

    sectors = [sector for sector in map(graze, [(8, 2.5, 2)])]  #, (10, 2.5, 2)])]  #, (12, 2.5, 2)])]

    def make_contours(data, limit):
        contours = measure.find_contours(data, limit)
        def make_poly(contour):
            #poly = np.flip(contour)  #Polygon(np.flip(contour))
            #poly = poly.simplify(4.0)
            return contour

        return map(make_poly, contours)


    m = 60
    x_vals = [m*10, m*10]
    y_vals = [m*2.4, m*4.5]

    plt.plot(x_vals, y_vals, linewidth=4, color=(0.0, 0.0, 0.0))

    # spline = interp1d(xy[0], xy[1])  # define function based on spline data points
    # line = interp1d(x_vals, y_vals)  # define function based on line data points
    #
    #
    #
    # f = lambda x: spline(x) - line(x)  # difference function, its zero marks the intersection
    # r = spopt.bisect(f, a=max(xy[0][0], x_vals[0]), b=min(xy[0][-1], x_vals[-1]))  # find root via bisection
    #
    # plt.scatter(r, spline(r))
    # print(r, spline(r))
    # plt.show()

    def do_refine(shape):
        x, y = shape.coords.xy
        tck, u = interpolate.splprep([x, y], s=0., k=2)  #, k=3)  #, k=0)
        ost = 0.001  #1/shape.length
        u_new = np.arange(0, 1.0+ost, ost)
        out = interpolate.splev(u_new, tck)

        # spline = interp1d(out[0], out[1], fill_value="extrapolate")
        # line = interp1d(x_vals, y_vals, fill_value="extrapolate")
        #
        # f = lambda gx: spline(gx) - line(gx)  # difference function, its zero marks the intersection
        # r = spopt.bisect(f, a=max(out[0][0], x_vals[0]), b=min(out[0][-1], x_vals[-1]))  # find root via bisection
        # plt.scatter(r, spline(r))

        plt.plot(out[0], out[1], linewidth=2, color=(1.0, 1.0, 1.0))
        # plt.plot(x, y, out[0], out[1], linewidth=1, color=(1.0, 1.0, 1.0))

        tout = b_spline_to_bezier_series(tck, False)
        print(tout)
        xy = [spl1[1:] for spl1 in tout]
        pp = np.concatenate([tout[0][:1]] + xy)
        print(pp)


        #
        # step_n = len(tout)
        # steps = np.zeros(step_n)
        # for n in range(step_n - 1):
        #     step = np.random.choice([-1, 0, 1], size=(1, 2))
        #     print(steps)
        #     steps = np.append(steps, step, axis=0)
        #     # something will be checked after each n
        #
        #
        #
        #
        #
        # xy = [spl1[:1] for spl1 in tout]
        # pp = np.concatenate([tout[0][:1]] + xy)

        # xy = [spl1[0] for spl1 in tout]


        # sx, sy = xy[:, 0], xy[:, 1]
        # plt.scatter(sx, sy, 'ro')
        # # sx, sy = tout[:, 0], tout[:, 1]
        # # plt.plot(sx, sy, 'ro', color=(1.0, 0.0, 1.0), lw=0.25)
        # print(tout)




    def draw_contours_for_level(data, gauss, simplification, interval):
        dlim = np.arange(1, 4000, interval)
        gb = gaussian_filter(data, sigma=.5)

        contours = measure.find_contours(data, 1)
        pass


    gb = gaussian_filter(dat, sigma=.5)
    conty = [p for p in make_contours(gb, -100)]

    for s in sectors:
        s_inside, s_outside = s
        plt.plot(*s_inside.exterior.xy, linewidth=2, color=(1.0, 0.2, 0.2))
        plt.plot(*s_outside.exterior.xy, linewidth=1, color=(0.2, 1.0, 0.2))

        for b in conty:
            ring_green = LineString(np.flip(b))
            ring_green = ring_green.simplify(0.5)

            if s_inside.intersects(ring_green):
                shell = s_outside.intersection(ring_green)
                if shell.type == 'MultiLineString':
                    for sub_shell in shell.geoms:

                        if len(sub_shell.coords) > 4:
                            do_refine(sub_shell)

                        plt.plot(*sub_shell.coords.xy, linewidth=1, color=(0.0, 0.0, 0.0))
                else:

                    if len(shell.coords) > 4:
                        do_refine(shell)

                    plt.plot(*shell.coords.xy, linewidth=1, color=(0.0, 0.0, 0.0))



    #
    # #dlim = [-20, -200, -2000]  #, -50, -100, -200, -500]
    # dlim = np.arange(1, 4000, 200)
    # gb = gaussian_filter(dat, sigma=.5)
    # #dlim = [-20]
    # for d in dlim:
    #     polygons = [p for p in make_contours(gb, -d)]
    #     for b in polygons:
    #         ring_green = LineString(np.flip(b))
    #         #ring_green = ring_green.simplify(1.0)
    #
    #         for s in sectors:
    #             s_inside, s_outside = s
    #             plt.plot(*s_inside.exterior.xy, linewidth=2, color=(1.0, 0.2, 0.2))
    #
    #
    #             if s_inside.intersects(ring_green):
    #                 shell = s_inside.intersection(ring_green)
    #                 if shell.type == 'MultiLineString':
    #                     for sub_shell in shell.geoms:
    #                         plt.plot(*sub_shell.coords.xy, linewidth=1, color=(1.0, 1.0, 0.0))
    #                 else:
    #                     plt.plot(*shell.coords.xy, linewidth=1, color=(1.0, 1.0, 0.0))

            #ringgreen = ringgreen.simplify(3.0)

            #
            # if ringgreen.length > 40.0:
            #     print(ringgreen.length)
            #
            #     plt.plot(*ringgreen.coords.xy, linewidth=1, color=(1.0, 1.0, 0.0))

                #do_refine(ringgreen)
                # print(ringgreen.length)
                #
                # x, y = b[:, 1], b[:, 0]
                #
                # plt.plot(x, y, linewidth=1, color=(1.0, 1.0, 1.0))
            # if len(b.exterior.coords) > 5:
            # #     ringgreen = LineString(list(b.exterior.coords))
            # #     do_refine(ringgreen)
            #
            #     plt.plot(*b.exterior.xy, linewidth=1, color=(1.0, 1.0, 1.0))

    # for b in sectors:
    #     s_inside, s_outside = b
    #     plt.plot(*s_inside.exterior.xy, linewidth=2, color=(1.0, 0.2, 0.2))
    #     plt.plot(*s_outside.exterior.xy, linewidth=1, color=(0.2, 1.0, 0.2))
    #
    #     for poly in polygons:
    #         if poly.area > 1.0:
    #             if s_inside.intersects(poly):
    #                 shell = s_outside.intersection(poly.exterior)
    #                 if shell.type == 'MultiLineString':
    #                     for sub_shell in shell.geoms:
    #
    #                         # # LinearRing of the polygons
    #                         # ringgreen = LineString(list(sub_shell.coords))
    #                         # ringblue = LineString(list(s_inside.exterior.coords))
    #                         # # Union of the rings
    #                         # sub_shell = ringgreen.union(ringblue)
    #                         # #print(sub_shell)
    #                         # sub_shell = ops.linemerge(sub_shell)
    #                         # for ss in sub_shell.geoms:
    #                         plt.plot(*sub_shell.coords.xy, linewidth=4, color=(0.0, 0.0, 0.0))
    #                         if len(sub_shell.coords) > 3:
    #                             do_refine(sub_shell)
    #                 else:
    #                     plt.plot(*shell.coords.xy, linewidth=4, color=(0.0, 0.0, 0.0))
    #                     do_refine(shell)
    #
    # sk, tk = graze((8, 4))
    # plt.plot(*sk.exterior.xy, linewidth=2, color=(1.0, 0.2, 0.2))
    # plt.plot(*tk.exterior.xy, linewidth=1, color=(0.2, 1.0, 0.2))

    plt.show()
    exit()

    #tk = box(5.8 * m, 5.8 * m, 7.2 * m, 7.2 * m)
    # Find contours at a constant value of 0.8
    # contours = measure.find_contours(dat, -100.0)
    #
    #
    #
    # for contour in contours:
    #
    #
    #
    #     polygon = Polygon(np.flip(contour))
    #     polygon = polygon.simplify(1.0)
    #
    #     if polygon.area > 1.0:
    #         if sk.intersects(polygon):
    #
    #             x, y = contour[:, 1], contour[:, 0]
    #
    #             pze = tk.intersection(polygon.exterior)
    #             if pze.type == 'MultiLineString':
    #                 for li in pze.geoms:
    #                     plt.plot(*li.coords.xy, linewidth=4, color=(0.0, 0.0, 0.0))
    #             else:
    #
    #                 x, y = pze.coords.xy
    #                 tck, u = interpolate.splprep([x,y], s=2)
    #                 unew = np.arange(0, 1.01, 0.01)
    #                 out = interpolate.splev(unew, tck)
    #                 plt.plot(x, y, out[0], out[1], linewidth=3, color=(0.0, 0.0, 0.0))
    #                 tout = b_spline_to_bezier_series(tck, False)
    #                 x, y = tout[:, 0], tout[:, 1]
    #                 plt.plot(x, y, 'ro')
    #
    #
    #
    #                 print(out)
    #                 #plt.plot(*pze.coords.xy, linewidth=4, color=(0.0, 0.0, 0.0))
    #                 #ctrl_points(pze.coords)
    #             # print(pze)
    #             # plt.plot(*pze.xy, linewidth=1, color=(0.0, 0.0, 0.0))
    #
    #
    #         plt.plot(*polygon.exterior.xy, linewidth=1, color=(1.0, 1.0, 1.0))

        # print(polygon)
        # break
        #
        # plt.plot(contour[:, 1], contour[:, 0], linewidth=1, color=(1.0, 1.0, 1.0))

    xt = np.arange(-6, 36.5, 1.0)
    yt = np.arange(30, 46, 1.0)
    # plt.xticks(xt)
    # plt.yticks(yt)
    # plt.colorbar()
    plt.show()
