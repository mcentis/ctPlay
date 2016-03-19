import numpy as np

class phantom:
    def __init__(self, listOfShapes):
        self.los = listOfShapes

    def draw(self, p): # p is plt 
        for s in self.los:
            s.draw(p)

    def rotate(self, ang):
        for s in self.los:
            s.rotate(ang)

    def calcSurvivors(self, p1, p2):
        absorption = 0.
        for s in self.los:
            absorption += s.calcPath(p1, p2) * s.mu
        return np.exp(-1 * absorption)

class circle:
    phi = np.arange(0, 2*np.pi, 0.01)
    def __init__(self, center, radius, absCoeff):
        self.c0 = center # not to be rotated
        self.c = center
        self.r = radius
        self.lastAngle = 0.0 # last rotation angle
        self.mu = absCoeff

    def draw(self, p): # p is plt
        x = self.c[0] + self.r * np.cos(self.phi)
        y = self.c[1] + self.r * np.sin(self.phi)
        p.plot(x, y)

    def rotate(self, ang):
        self.lastAngle = ang
        m = np.matrix([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        self.c = m.dot(self.c0).getA1() # the getA1 gives the result as simple array

    def calcPath(self, p1, p2): # see http://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
        A = np.array(p1)
        B = np.array(p2)
        C = self.c
        R = self.r
        V = B - A # line parametrized as P = A + tV with 0 <= t <= 1
        U = A - C
        # solution of 2nd grade equation in t to find intersections (form at^2 + bt + c = 0)
        a = V.dot(V)
        b = 2 * V.dot(U)
        c = U.dot(U) - R * R
        discriminant = b * b - 4 * a * c
        if discriminant < 0:
            return 0
        t1 = (-b - np.sqrt(discriminant)) / (2 * a)
        t2 = (-b + np.sqrt(discriminant)) / (2 * a)
        # it follows that t2 >= t1

        # 3x HIT cases:
        #          -o->             --|-->  |            |  --|->
        # Impale(t1 hit,t2 hit), Poke(t1 hit,t2>1), ExitWound(t1<0, t2 hit), 

        # 3x MISS cases:
        #       ->  o                     o ->              | -> |
        # FallShort (t1>1,t2>1), Past (t1<0,t2<0), CompletelyInside(t1<0, t2>1)

        if 0 <= t1 and t1 <= 1 and 0 <= t2 and t2 <= 1: # the ray traverse the circle
            return np.sqrt((t2 * V - t1 * V).dot(t2 * V - t1 * V))

        if 0 <= t1 and t1 <= 1 and t2 > 1:
            return np.sqrt((V - t1 * V).dot(V - t1 * V))

        if t1 < 0 and 0 <= t2 and t2 <= 1:
            return np.sqrt((t2 * V).dot(t2 * V))

        # if not in the 3 above
        return 0
