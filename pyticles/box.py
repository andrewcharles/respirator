"""
    Extends the controller class to apply position determined constraints,
    usually in the form of box boundary collisions or periodic translation.

    This won't handle periodic boundary interactions, but will possibly
    contain the code for computing the nearest pair spacing.
    
    Copyright Andrew Charles 2008
    All rights reserved.
    This code is licensed under the terms of the new BSD license.

"""

import controller

class Box(controller.Controller):
    
    def __init__(self,p='none',xmax=64,ymax=48,zmax=100):
        self.p = p
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax

    def apply(self):
        print "Do nothing"
        return


class Box2d(controller.Controller):
    
    def __init__(self,p='none',xmax=64,ymax=48):
        self.p = p
        self.xmax = xmax
        self.ymax = ymax

    def apply(self):
        print "Do nothing"
        return


class PeriodicBox(Box):

    def min_image(r):
        print "does not exist yet"

    def apply(self,p):
         """ applies periodic boundaries """
         for i in range(p.n):
            if p.r[i,0] > self.xmax:
                p.r[i,0] = 0
            if p.r[i,0] < 0:
                p.r[i,0] = self.xmax 
            if p.r[i,1] > self.ymax:
                p.r[i,1] = 0
            if p.r[i,1] < 0:
                p.r[i,1] =self.ymax 
            if p.r[i,2] > self.zmax:
                p.r[i,2] = 0
            if p.r[i,2] < 0:
                p.r[i,2] =self.zmax


class PeriodicBox2d(Box2d):

    def min_image(r):
        print "does not exist yet"

    def apply(self,p):
         """ applies periodic boundaries """
         for i in range(p.n):
            if p.r[i,0] > self.xmax:
                p.r[i,0] = 0
            if p.r[i,0] < 0:
                p.r[i,0] = self.xmax 
            if p.r[i,1] > self.ymax:
                p.r[i,1] = 0
            if p.r[i,1] < 0:
                p.r[i,1] =self.ymax 


class MirrorBox(Box):

    def apply(self,p):
         """ applies mirror boundaries """
         for i in range(p.n):
            if p.r[i,0] > self.xmax:
                p.r[i,0] =self.xmax 
                p.v[i,0] = -p.v[i,0]
            if p.r[i,0] < 0:
                p.r[i,0] = 0
                p.v[i,0] = -p.v[i,0]

            if p.r[i,1] > self.ymax:
                p.r[i,1] = self.ymax 
                p.v[i,1] = -p.v[i,1]
            if p.r[i,1] < 0:
                p.r[i,1] = 0
                p.v[i,1] = -p.v[i,1]

            if p.r[i,2] > self.zmax:
                p.r[i,2] = self.zmax 
                p.v[i,2] = -p.v[i,2]
            if p.r[i,2] < 0:
                p.r[i,2] = 0
                p.v[i,2] = -p.v[i,2]


class MirrorBox2d(Box2d):

    def apply(self,p):
         """ applies mirror boundaries """
         for i in range(p.n):
            if p.r[i,0] > self.xmax:
                p.r[i,0] =self.xmax 
                p.v[i,0] = -p.v[i,0]
            if p.r[i,0] < 0:
                p.r[i,0] = 0
                p.v[i,0] = -p.v[i,0]

            if p.r[i,1] > self.ymax:
                p.r[i,1] = self.ymax 
                p.v[i,1] = -p.v[i,1]
            if p.r[i,1] < 0:
                p.r[i,1] = 0
                p.v[i,1] = -p.v[i,1]

