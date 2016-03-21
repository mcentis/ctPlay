import numpy as np
from math import radians
from scipy.fftpack import fft, ifft, fftfreq

class backProjection:
    def __init__(self, sino, dist, voxSize):
        self.sinogram = sino
        self.vSize = float(voxSize)
        self.dist = abs(float(dist))
        self.nVox = int(2 * self.dist / float(voxSize)) # reconstruction grid dimension defined by the rotation center-detector distance
        # outer dimension of the matrix in each dim: -dist -> -dist + nVox * vSize
        self.vMatrix = []
        for j in range(self.nVox):
            self.vMatrix.append([])
            for i in range(self.nVox):
                self.vMatrix[-1].append(0.)

    def rotM(self, ang):
        return np.matrix([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])

    def reconstruct(self):
        linPix = 10
        h = self.sinogram[0]
        yEdges = self.sinogram[1]
        sensEd = self.sinogram[2]
        for iAng in range(len(yEdges) -1): # loop on angles
            ang = yEdges[iAng]
            print 'Projecting angle ', ang
            for iPix in range(len(sensEd) - 1): # loop on pixels
                span = np.arange(sensEd[iPix], sensEd[iPix + 1], linPix)
                weight = h[iAng][iPix] / len(span)
                for l in span:
                    p1 = np.array([self.dist, l]).dot(self.rotM(-ang)).getA1() # select the points as opposite
                    p2 = np.array([self.dist, -l]).dot(self.rotM(-ang + radians(180))).getA1()
                    self.traceRay(p1, p2, weight)
        self.makeImage()

    def filterSino(self):
        h = self.sinogram[0]
        yEdges = self.sinogram[1]
        sensEd = self.sinogram[2]
        fre = fftfreq(len(h[0]), sensEd[1] - sensEd[0])
        pos = sensEd[:-1]
        ang = yEdges[:-1]
        y = []
        x = []
        w = []
        for row, a in zip(h, ang):
            rt = fft(row)
            rtFilt = []
            for r, f in zip(rt, fre):
                rtFilt.append(r * abs(f)) # ramp filter
                #rtFilt.append(r * min(1.0, abs(f))) # other filter
            rit = ifft(rtFilt)
            for r, p in zip(rit, pos):
                y.append(a)
                x.append(p)
                w.append(np.real(r)) # real part, the imaginary one should be really small
        self.sinogram = np.histogram2d(y, x, [len(ang), len(pos)], [[yEdges[0], yEdges[-1]], [sensEd[0], sensEd[-1]]], False, w)

    def oneMinusSino(self):
        h = self.sinogram[0]
        yEdges = self.sinogram[1]
        sensEd = self.sinogram[2]
        fre = fftfreq(len(h[0]), sensEd[1] - sensEd[0])
        pos = sensEd[:-1]
        ang = yEdges[:-1]
        y = []
        x = []
        w = []
        for row, a in zip(h, ang):
            for r, p in zip(row, pos):
                y.append(a)
                x.append(p)
                w.append(1.0-r)
        self.sinogram = np.histogram2d(y, x, [len(ang), len(pos)], [[yEdges[0], yEdges[-1]], [sensEd[0], sensEd[-1]]], False, w)

    def makeImage(self):
        x = []
        y = []
        c = []
        for i in range(self.nVox):
            for j in range(self.nVox):
                c.append(self.vMatrix[j][i])
                x.append(-self.dist + i * self.vSize)
                y.append(-self.dist + (self.nVox -1 -j) * self.vSize)
        self.image = np.histogram2d(y, x, [self.nVox, self.nVox], [[-self.dist, -self.dist + self.nVox * self.vSize], [-self.dist, -self.dist + self.nVox * self.vSize]], False, c)

    # ray tracing from: Siddon 1984 Fast calculation of the exact radiological path for a three dimensional ct array
    # the rays are expressed in each dimension like x = x1 + t(x2 -x1) with 0 <= t <= 1 in each dimension

    def traceRay(self, p1, p2, weight): # follows ray into matrix, change matrix elements accordingly
        # min and max coordinate of planes
        pmin = -self.dist
        pmax = -self.dist + self.nVox * self.vSize

        tmin = []
        tmax = []
        for i in range(len(p1)):
            tmin.append(np.nan) # undefined
            tmax.append(np.nan) # undefined

        for i in range(len(p1)): # determine min and max parameters when possible
            if p2[i] - p1[i] != 0:
                tmin[i] = (pmin - p1[i]) / (p2[i] - p1[i])
                tmax[i] = (pmax - p1[i]) / (p2[i] - p1[i])

        tminAll = np.nanmax([0, np.nanmin([tmin[0], tmax[0]]), np.nanmin([tmin[1], tmax[1]])])
        tmaxAll = np.nanmin([1, np.nanmax([tmin[0], tmax[0]]), np.nanmax([tmin[1], tmax[1]])])
        
        if tminAll >= tmaxAll: # no intersection with matrix (holds only for cases where ray is not parallel to the axis or both points are outside the volume)
            return
        
        # the ray is parallel to one of the axis
        if np.isnan(tmin[1]): # ray parallel to y
            if p1[1] < pmin or p1[1] > pmax:
                return
        if np.isnan(tmin[0]): # ray parallel to x
            if p1[0] < pmin or p1[0] > pmax:
                return

        intersections = [] # list of the t where a plane is encountered
        for i in range(len(p1)):
            if p2[i] - p1[i] == 0:
                continue
            # determine the min and max indexes for each direction
            if p2[i] - p1[i] > 0:
                imin = self.nVox - int((pmax - tminAll * (p2[i] - p1[i]) - p1[i]) / self.vSize)
                imax = int((p1[i] + tmaxAll * (p2[i] - p1[i]) - pmin) / self.vSize)
            else:
                imin = self.nVox - int((pmax - tmaxAll * (p2[i] - p1[i]) - p1[i]) / self.vSize)
                imax = int((p1[i] + tminAll * (p2[i] - p1[i]) - pmin) / self.vSize)

            for j in range(imin, imax + 1): # determine the t for each intersection
                intersections.append((-self.dist + j * self.vSize - p1[i])/(p2[i] - p1[i]))
        
        intersections.append(tminAll)
        intersections.append(tmaxAll)
        intersections = sorted(set(intersections)) # remove duplicates and reorder

        # determine pathlength and indexes for each pair of intersections
        indexLen = []
        np1 = np.array(p1) # to be sure that the the objects are np arrays
        np2 = np.array(p2)
        diff = np2 - np1
        dist = np.sqrt(diff.dot(diff))
        for i in range(1, len(intersections)):
            fmid = (intersections[i] + intersections[i - 1]) / 2
            ile = []
            for j in range(len(p1)): # indexes
                index = int((p1[j] + fmid * (p2[j] - p1[j]) - pmin) / self.vSize)
                ile.append(index)
            ile.append(dist * (intersections[i] - intersections[i - 1]))
            indexLen.append(ile)

        #operation on the matrix 
        for il in indexLen:
            if il[0] < len(self.vMatrix) and il[1] < len(self.vMatrix): # the range of planes is [0, nVox], while the matrix is [0, nVox)
                self.vMatrix[len(self.vMatrix) -1 -il[1]][len(self.vMatrix) -1 -il[0]] += il[2] * weight / dist
