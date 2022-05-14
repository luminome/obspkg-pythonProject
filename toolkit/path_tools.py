#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap

from scipy import interpolate

from shapely import ops
from shapely.affinity import translate, scale
from shapely.geometry import Polygon, MultiPolygon, Point, box, LineString, LinearRing
from skimage import measure
from skimage import transform
from scipy.ndimage.filters import gaussian_filter


def get_contour(data, sigma, sca, level, density, pos_origin):
    # g_data = gaussian_filter(data, sigma=sigma)

    epm = transform.pyramid_expand(data, upscale=sca, sigma=0.5, preserve_range=True)

    f_contours = measure.find_contours(data, level)
    contours = []
    for ep in f_contours:
        ep = LineString(np.flip(ep))
        ep = scale(ep, xfact=1 / density, yfact=-1 / density, origin=(0, 0))
        ep = translate(ep, xoff=pos_origin[0], yoff=pos_origin[1])
        contours.append(ep)

    return contours


def create_poly(x, y):
    p_shap = Polygon(x, y)
    return p_shap


def do_refine(x, y, s=0, k=2):
    # x, y = shape.coords.xy
    tck, u = interpolate.splprep([x, y], s=s, k=k)  # , k=3)  #, k=0)
    ost = 0.00001  # 1/shape.length
    u_new = np.arange(0, 1.0 + ost, ost)
    out = interpolate.splev(u_new, tck)
    # spline = interp1d(out[0], out[1], fill_value="extrapolate")
    # line = interp1d(x_vals, y_vals, fill_value="extrapolate")
    # f = lambda gx: spline(gx) - line(gx)  # difference function, its zero marks the intersection
    # r = spopt.bisect(f, a=max(out[0][0], x_vals[0]), b=min(out[0][-1], x_vals[-1]))  # find root via bisection
    # plt.scatter(r, spline(r))
    return out


def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection

