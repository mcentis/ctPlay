from phantom import *
import numpy as np

class ctScan:
    def __init__(self, phan, angStart, angStop, angStep, sensStart, sensStop, pitch, sensDist):
        self.pha = phan # phantom
        self.nAngles = int((angStop - angStart)/float(angStep))
        self.a1 = angStart
        self.a2 = angStart + self.nAngles * angStep
        self.aStep = float(angStep)
        self.sen1 = sensStart
        self.nPix = int((sensStop - sensStart)/float(pitch))
        self.sen2 = sensStart + pitch * self.nPix
        self.senP = float(pitch)
        self.senD = sensDist # distance between sensor and phantom center

    def scan(self):
        detPos = []
        cont = []
        ang = []
        for iAng in range(self.nAngles):
            a = self.a1 + iAng * self.aStep
            print 'Scanning angle ', a
            self.pha.rotate(a)
            pos, con = self.calcProjection()
            for p, c in zip(pos, con):
                detPos.append(p)
                cont.append(c)
                ang.append(a)
        self.sinogram = np.histogram2d(ang, detPos, [self.nAngles, self.nPix], [[self.a1,  self.a2], [self.sen1, self.sen2]], False, cont)

    def calcProjection(self):
        linesPixel = 10 # number of line of response per pixel
        detPos = []
        cont = []
        for pix in range(self.nPix):
            # problem: the length of span is not constant!!!
            span = np.arange(self.sen1 + pix * self.senP, self.sen1 + (pix + 1) * self.senP, self.senP/linesPixel) # height of the lines pixel
            intensity = 0.
            for l in span:
                intensity += self.pha.calcSurvivors([-self.senD, l], [self.senD, l]) # line from -senD to senD with height l
            intensity /= len(span) # normalization, if no  material in phantom, this is 1
            detPos.append(self.sen1 + pix * self.senP )
            cont.append(intensity)
        return detPos, cont
