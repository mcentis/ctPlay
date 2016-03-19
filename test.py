from phantom import *
from ctScan import *
from backProjection import *

import matplotlib.pyplot as plt
import numpy as np
from math import *

lis = [circle(np.array([1, 5]), 2, 0.5), circle(np.array([-1, -3]), 1, 0.2)]
#lis = [circle(np.array([-1, -3]), 1, 2)]
pa = phantom(lis)

plt.figure(1)
plt.xlim(-15, 15)
plt.ylim(-15, 15)
pa.draw(plt)

ct = ctScan(pa, radians(0), radians(180), radians(1), -10, 10, 0.15, 15)
ct.scan()
h = ct.sinogram[0]
yEdges = ct.sinogram[1]
xEdges = ct.sinogram[2]

fig = plt.figure(2)
im = plt.imshow(h, interpolation='none', origin='low', extent=[xEdges[0], xEdges[-1], yEdges[0], yEdges[-1]], aspect='auto')
plt.colorbar(im)

rec = backProjection(ct.sinogram, 15, 0.15)
rec.reconstruct()

h = rec.image[0]
yEdges = rec.image[1]
xEdges = rec.image[2]

fig = plt.figure(3)
im = plt.imshow(h, interpolation='none', origin='low', extent=[xEdges[0], xEdges[-1], yEdges[0], yEdges[-1]], aspect='auto')
plt.colorbar(im)


plt.show()
