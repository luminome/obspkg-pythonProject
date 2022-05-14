# -*- coding: utf-8 -*-
"""
Example demonstrating a variety of scatter plot features.
"""



## Add path to library (just for examples; you do not need this)
## import initExample

from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid


# from collections import namedtuple
# from itertools import chain

app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.resize(800,800)
view = pg.GraphicsLayoutWidget()  ## GraphicsView with GraphicsLayout inserted by default
mw.setCentralWidget(view)
mw.show()
mw.setWindowTitle('pyqtgraph example: ScatterPlot')

## create four areas to add plots
w1 = view.addPlot()
# w2 = view.addViewBox()
# w2.setAspectLocked(True)
# view.nextRow()
# w3 = view.addPlot()
# w4 = view.addPlot()
print("Generating data, this takes a few seconds...")

## Make all plots clickable
clickedPen = pg.mkPen('b', width=2)
lastClicked = []
def clicked(plot, points):
    global lastClicked
    for p in lastClicked:
        p.resetPen()
    print("clicked points", points)
    for p in points:
        p.setPen(clickedPen)
    lastClicked = points

med_marine_regions = gpd.read_file("data/goas_med/club_med.shp")

print(med_marine_regions.head(10))
#//med_marine_regions = med_marine_regions.simplify(1, True)

med_poly = med_marine_regions["geometry"][1]

#fig, ax = plt.subplots(figsize=(10, 6))

print('has', med_poly.geom_type, med_poly.area, med_poly.exterior.geom_type, med_poly.exterior.length)

# for n in med_poly.interiors:
#     #print(n.geom_type, n.length)
#     ns = n.simplify(0.008, True)
#     ns = ns.buffer(0.005)
#     nf = ns.exterior.simplify(0.002, True)
#     pg.plot(*nf.xy, color="blue")

#mext = med_poly.exterior.simplify(1, True)
exterior_coords = med_poly.exterior.coords[:]
print(len(exterior_coords))

    #pg.plot(*med_poly.exterior.xy, pen="r")

## There are a few different ways we can draw scatter plots; each is optimized for different types of data:
## 1) All spots identical and transform-invariant (top-left plot).
## In this case we can get a huge performance boost by pre-rendering the spot
## image and just drawing that image repeatedly.

n = 300
s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(None), brush=pg.mkBrush(255, 255, 255, 120))
pos = med_poly.exterior.xy  #np.random.normal(size=(2,n), scale=1e-5)
#spots = [{'pos': pos[:,i], 'data': 1} for i in range(n)] + [{'pos': [0,0], 'data': 1}]
s1.addPoints(pos)
w1.addItem(s1)
s1.sigClicked.connect(clicked)


## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
