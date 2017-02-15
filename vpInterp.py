#!/usr/bin/env python

#This will be a general script that will have different interpolation and 
# extraplations schemes, these different schemes will be represented
# as different classes within this utilityModule


# Import modules
import math
import cmath
import sys
import string
import numpy

# scipy's inpterolation
from scipy import interpolate


class BilinearInterp2DCart:

    def __init__(self, x0, y0, fx ):
        
        # potential optimization
        # have the file read in here so for each
        # time vorpal calls the pyfunc this class
        # creates it only has to read in data once
        # so overall O(n) ~ complexity
        # but we will see if this is even needed
        # once i implement the full program flow
        # of micro ring resonator
        # ALSO: it would kill this interpolation 
        # class's generality
        
        # save input parameters
        self.X0 = x0
        self.Y0 = y0
        self.Fx = fx

        # derived parameters
        self.numNodesX = len(self.X0)
        self.numNodesY = len(self.Y0)
        self.numCellsX = self.numNodesX - 1
        self.numCellsY = self.numNodesY - 1
        self.lengthX   = self.X0[self.numNodesX - 1] - self.X0[0]
        self.lengthY   = self.Y0[self.numNodesY - 1] - self.Y0[0] 
        self.DX        = self.lengthX / self.numCellsX
        self.DY        = self.lengthY / self.numCellsY


    def getNearestX(self, x ):
        xx = ( x - self.X0[0] )
        index_x = int( xx / self.DX )
        return index_x

    def getNearestY(self, y ):
        yy = ( y - self.Y0[0] )
        index_y = int( yy / self.DY )
        return index_y

    def __call__(self, x, y):

        i = self.getNearestX(x)
        j = self.getNearestY(y)

        if( i == self.numCellsX ) or ( j == self.numCellsY ) :
            return self.Fx[i][j]

        if( i > self.numCellsX ) or ( j > self.numCellsY ) :
            return 0.00


        x0 = self.X0[i]
        x1 = self.X0[i+1]
        
        y0 = self.Y0[j]
        y1 = self.Y0[j+1]

        w11 = self.Fx[i][j]
        w12 = self.Fx[i][j+1]
        w21 = self.Fx[i+1][j]
        w22 = self.Fx[i+1][j+1]

        return ((w11* (x1 - x) * (y1 - y) +
                w12 * (x - x0) * (y1 - y) +
                w21 * (x1 - x) * (y - y0) +
                w22 * (x - x0) * (y - y0))/ ((x1 - x0)*(y1 - y0)))



class RectBivariateSpline:

        def __init__(self, x0, y0, fx ):
        
            # save input parameters
            self.X0 = x0
            self.Y0 = y0
            self.Fx = fx

            # derived parameters
            self.numNodesX = len(self.X0)
            self.numNodesY = len(self.Y0)
            self.numCellsX = self.numNodesX - 1
            self.numCellsY = self.numNodesY - 1
            self.lengthX   = self.X0[self.numNodesX - 1] - self.X0[0]
            self.lengthY   = self.Y0[self.numNodesY - 1] - self.Y0[0] 
            self.DX        = self.lengthX / self.numCellsX
            self.DY        = self.lengthY / self.numCellsY

            # scipy RectBivariateSpline
            self.interpMode = interpolate.RectBivariateSpline( self.X0, self.Y0, self.Fx )

        def __call__(self, x, y):

            return numpy.float(self.interpMode( x , y ))
        

        
