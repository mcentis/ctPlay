from phantom import *
from ctScan import *
from backProjection import *

import matplotlib.pyplot as plt
import numpy as np
from math import *

def showStuff(inp):
    h = inp[0]
    yEdges = inp[1]
    xEdges = inp[2]
    im = plt.imshow(h, interpolation='none', origin='low', extent=[xEdges[0], xEdges[-1], yEdges[0], yEdges[-1]], aspect='auto')#, cmap='Greys_r')
    plt.colorbar(im)

lis = [circle(np.array([1, 5]), 2, 0.5), circle(np.array([-1, -3]), 1, 0.2), circle(np.array([-3, 0]), 0.1, 5)]
#lis = [circle(np.array([-1, -3]), 1, 2)]
pa = phantom(lis)

plt.figure(1)
plt.xlim(-15, 15)
plt.ylim(-15, 15)
pa.draw(plt)

ct = ctScan(pa, radians(0), radians(180), radians(1), -10, 10, 0.15, 15)

ct.scan()
fig = plt.figure(2)
showStuff(ct.sinogram)

rec = backProjection(ct.sinogram, 15, 0.15)

rec.oneMinusSino()
fig = plt.figure(3)
showStuff(rec.sinogram)

rec.reconstruct()
fig = plt.figure(4)
showStuff(rec.image)
# h = rec.sinogram[0]
# ed = rec.sinogram[2]
# plt.plot(ed[:-1], h[0])

rec.filterSino()
fig = plt.figure(5)
showStuff(rec.sinogram)
# h = rec.sinogram[0]
# ed = rec.sinogram[2]
# plt.plot(ed[:-1], h[0])

rec.reconstruct()
fig = plt.figure(6)
showStuff(rec.image)

plt.show()
