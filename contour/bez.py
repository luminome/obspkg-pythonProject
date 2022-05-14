import numpy as np
from scipy.interpolate import splprep


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


      from http://mail.scipy.org/pipermail/scipy-dev/2007-February/006651.html

    """
    from scipy.interpolate.fitpack import insert
    from numpy import asarray, unique, split, sum, transpose
    t, c, k = tck
    t = asarray(t)
    try:
        c[0][0]
    except:
        # I can't figure out a simple way to convert nonparametric splines to
        # parametric splines. Oh well.
        raise TypeError("Only parametric b-splines are supported.")
    new_tck = tck
    if per:
        # ignore the leading and trailing k knots that exist to enforce periodicity
        knots_to_consider = unique(t[k:-k])
    else:
        # the first and last k+1 knots are identical in the non-periodic case, so
        # no need to consider them when increasing the knot multiplicities below
        knots_to_consider = unique(t[k + 1:-k - 1])
    # For each unique knot, bring it's multiplicity up to the next multiple of k+1
    # This removes all continuity constraints between each of the original knots,
    # creating a set of independent Bezier curves.
    desired_multiplicity = k + 1
    for x in knots_to_consider:
        current_multiplicity = sum(t == x)
        remainder = current_multiplicity % desired_multiplicity
        if remainder != 0:
            # add enough knots to bring the current multiplicity up to the desired multiplicity
            number_to_insert = desired_multiplicity - remainder
            new_tck = insert(x, new_tck, number_to_insert, per)
    tt, cc, kk = new_tck
    # strip off the last k+1 knots, as they are redundant after knot insertion
    bezier_points = transpose(cc)[:-desired_multiplicity]
    if per:
        # again, ignore the leading and trailing k knots
        bezier_points = bezier_points[k:-k]
    # group the points into the desired bezier curves
    return split(bezier_points, len(bezier_points) / desired_multiplicity, axis=0)


def get_bezier_path(x, y, per=False):
    tck, uout = splprep([x, y], s=0., k=3, per=per)
    spl = b_spline_to_bezier_series(tck, per=per)

    # to mpl path
    from matplotlib.path import Path
    xy = [spl1[1:] for spl1 in spl]
    pp = np.concatenate([spl[0][:1]] + xy)
    codes = [Path.MOVETO] + [Path.CURVE4] * (3 * len(xy))

    p = Path(pp, codes)
    return p


if __name__ == '__main__':
    # from mpl_toolkits.axes_grid.axis_artist import BezierPath
    from matplotlib.patches import PathPatch

    theta = np.linspace(0, 2 * np.pi, 13)
    x = np.cos(theta)
    y = np.sin(theta)


    def myplot(ax, x, y, per=False):
        l1, = ax.plot(x, y, "rs", ms=10)
        p = get_bezier_path(x, y, per=per)
        l2, = ax.plot(p.vertices[:, 0], p.vertices[:, 1], "bo")
        bp = PathPatch(p, ec="b", fc="none")
        ax.add_patch(bp)

        ax.legend([l1, l2], ["input", "control points"])


    import matplotlib.pyplot as plt

    plt.clf()
    ax = plt.subplot(121)
    myplot(ax, theta, y)
    ax = plt.subplot(122)
    myplot(ax, x, y, per=True)
    plt.show()


