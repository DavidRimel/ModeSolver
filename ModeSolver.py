
import sys

#sys.path.append('../EMpy')# what is the best way to get EMpy?
import EMpy

import vpInterp #not working in some cases
# just going to use scipy's inpterp
# for now ......
#from scipy import interpolate


import numpy
import scipy


# ==================================================
# Utility Functions
# ==================================================

#-------------------------
# get edge values of nodal 1-D grid
def centered1d(x):
    return (x[1:] + x[:-1]) / 2.

#-------------------------
# get edge values of nodal 2-D grid
def centered2d(x):
    return (x[1:, 1:] + x[1:, :-1] + x[:-1, 1:] + x[:-1, :-1]) / 4.

#-------------------------
# This reads the modes in from a .txt file  
def readModeInFromTxt(modeNum = 0):
        
    # get grid info and define Nodal Grid 
    gridInfo = numpy.fromfile('Mode'+str(modeNum)+'_GridInfo.txt', sep = " ").reshape(3,2)
    numCells = gridInfo[0]
    lengths = gridInfo[1]
    startPositions = gridInfo[2]
    grid = ModeSolver(numCells, lengths, startPositions, None)
    x = grid.x
    y = grid.y
    
    # get Ex and Ey Mode Profiles from text files
    mode_Ex = numpy.fromfile('Mode'+str(modeNum)+'_Ex.txt', sep = " ").reshape(int(numCells[0]),int(numCells[1]-1))
    mode_Ey = numpy.fromfile('Mode'+str(modeNum)+'_Ey.txt', sep = " ").reshape(int(numCells[0]-1),int(numCells[1]))
    
    # create Ex yeeGrid from Nodal Grid definition 
    x_Ex_FDTD = EMpy.utils.centered1d(x)
    y_Ex_FDTD = y[1:-1]    

    # create Ey yeeGrid from Nodal Grid definition
    x_Ey_FDTD = x[1:-1]
    y_Ey_FDTD = EMpy.utils.centered1d(y)
    
    Ex = [x_Ex_FDTD, y_Ex_FDTD, mode_Ex]
    Ey = [x_Ey_FDTD, y_Ey_FDTD, mode_Ey]
    
    return  [ Ex, Ey ]





class ModeSolver:

# ========================================================
# Initialized the object based on input parameters
# defines the epsFunc that needs to be input into
# EMpy's mode solver, and defines grid in which
# mode is solved on.
# ========================================================
    def __init__(self, numCells, lengths, startPositions, epsFunc):

        # Save Input Parameters
        self.numCells = numCells.astype(int)
        self.lengths = lengths
        self.startPositions = startPositions

        # Derived Parameters
        self.numNodes = numCells + 1
        self.endPositions = self.startPositions + self.lengths
        self.DX = ( self.lengths[0] / (self.numCells[0]) )
        self.DY = ( self.lengths[1] / (self.numCells[1]) )
        self.x = numpy.linspace(self.startPositions[0], self.endPositions[0],
                           self.numNodes[0])
        self.y = numpy.linspace(self.startPositions[1], self.endPositions[1],
                           self.numNodes[1])
        
        # use these for a way to check 
        # if solve() has been performed
        # on this obj.
        self.modes = None
        self.neigs = None

        #EMpyEpsFunc : function
        #    This is a function that provides the relative permittivity
        #    matrix (square of the refractive index) as a function of its
        #    x and y numpy.arrays (the function's input parameters). 
        #    The function must be of the form:
        #    ``myRelativePermittivity(x,y)``
        #    The function returns a relative permittivity numpy.array 
        #    of shape( x.shape[0], y.shape[0] ) where each element of 
        #    the array can either be a single float, corresponding the
        #    an isotropic refractive index,or, a length-5 tuple. In the
        #    tuple case, the relative permittivity is given in the form
        #    (epsxx, epsxy, epsyx, epsyy, epszz).   
        def EMpyEpsFunc( x , y ):
            xLength = x.shape[0]
            yLength = y.shape[0]
            epsMatrix = numpy.zeros([xLength,yLength])
            index_x = 0
            for i in x:
                index_y = 0
                for j in y:
                    epsMatrix[index_x][index_y] = epsFunc( i, j )
                    index_y +=1
                index_x +=1
            return epsMatrix

        # stFunc that EMpy's mode solver needs
        self.epsFunc = EMpyEpsFunc



# ========================================================
# This function finds the eignemodes
#
#Parameters
#----------
#waveL : float
#    wave length of the mode you mode you want solved.
#neigs : int
#    number of eigenmodes to find
#tol : float
#    Relative accuracy for eigenvalues. value of 0 implies
#    machine precision.
#guess : float
#    a guess for the refractive index. Only finds eigenvectors with an 
#    effective refrative index
#    higher than this value.
#modeNum : int
#   the mode corresponding the the mode number you want returned
#
#Returns
#----------
#  self.modes = FDModes object defined in [EMpy => FD.py]
#      To get modes Example:
#             self.modes[i].Ex ~ i'th Ex mode
#             self.modes[i].Ey ~ i'th Ey mode
#             self.modes[i].Ez ~ i'th Ez mode
#   NOTE: not defined on Yee Grid!!!!
# ========================================================
    def solve( self, waveL, neigs, tol, guess=None ):
        self.neigs = neigs
        self.modes = EMpy.modesolvers.FD.VFDModeSolver(waveL,self.x,self.y,self.epsFunc,'0000').solve( neigs, tol ).modes
        return self.modes


# ========================================================
# EMpy interally interpolates the Ex,Ey,and Ez values of the
# modeNum'th mode to a standard yeeGrid
#
# Returns
# --------- )
# [ [ x_Ex, y_Ex, mode_Ex ], [x_Ey, y_Ey, mode_Ey ],...]
# ========================================================
    def interpolateModeToYeeGrid(self,modeNum):
        
        #need to test
        if( self.modes == None and self.neigs == None ):
            print "Need to solve for the modes before before inperpolating them to YeeGrid!\n "
            exit(0)

        #need to test
        if( self.neigs < modeNum ):
            print "Need to solve for more modes to get this mode number!\n "
            exit(0)
    
        mode = self.modes[modeNum]
        mode_yeeGrid = mode.get_fields_for_FDTD()
        #Ex
        x_Ex_yeeGrid = mode.x_Ex_FDTD
        y_Ex_yeeGrid = mode.y_Ex_FDTD
        Ex_yeeGrid = [x_Ex_yeeGrid, y_Ex_yeeGrid, mode_yeeGrid[0]]
        #Ey
        x_Ey_yeeGrid = mode.x_Ey_FDTD
        y_Ey_yeeGrid = mode.y_Ey_FDTD
        Ey_yeeGrid = [x_Ey_yeeGrid, y_Ey_yeeGrid, mode_yeeGrid[1]]
        #Ez
        x_Ez_yeeGrid = mode.x_Ez_FDTD
        y_Ez_yeeGrid = mode.y_Ez_FDTD
        Ez_yeeGrid = [x_Ez_yeeGrid, y_Ez_yeeGrid, mode_yeeGrid[2]]
        #Hx
        x_Hx_yeeGrid = mode.x_Hx_FDTD
        y_Hx_yeeGrid = mode.y_Hx_FDTD
        Hx_yeeGrid = [x_Hx_yeeGrid, y_Hx_yeeGrid, mode_yeeGrid[3]]
        #Hy
        x_Hy_yeeGrid = mode.x_Hy_FDTD
        y_Hy_yeeGrid = mode.y_Hy_FDTD
        Hy_yeeGrid = [x_Hy_yeeGrid, y_Hy_yeeGrid, mode_yeeGrid[4]]
        #Hz
        x_Hz_yeeGrid = mode.x_Hz_FDTD
        y_Hz_yeeGrid = mode.y_Hz_FDTD
        Hz_yeeGrid = [x_Hz_yeeGrid, y_Hz_yeeGrid, mode_yeeGrid[5]]

        
        return [ Ex_yeeGrid, Ey_yeeGrid, Ez_yeeGrid, Hx_yeeGrid, Hy_yeeGrid, Hz_yeeGrid ]


# ========================================================
# This Writes the modes out to a .txt file 
# ========================================================        
    def writeModeOutToTxt(self, modeNum = 0):
                    
        modesOnYeeGrid = self.interpolateModeToYeeGrid(modeNum)
        
        Ex = modesOnYeeGrid[0]
        Ey = modesOnYeeGrid[1]
        
        # Create Ex_Mode data file
        #x_Ex = Ex[0]
        #y_Ex = Ex[1]
        mode_Ex = Ex[2] 
        numpy.savetxt('Mode'+str(modeNum)+'_Ex.txt', numpy.real(mode_Ex), fmt='%.26e')

        # Create Ey_Mode data file
        #x_Ey = Ey[0]
        #y_Ey = Ey[1]
        mode_Ey = Ey[2]
        numpy.savetxt('Mode'+str(modeNum)+'_Ey.txt', numpy.real(mode_Ey), fmt='%.26e')

        # Create yeeGrid data file using: numCells, lengths, startPositions
        numCells = self.numCells
        lengths = self.lengths
        startPositions = self.startPositions
        gridInfo = numpy.array([ numCells, lengths, startPositions ])
        numpy.savetxt('Mode'+str(modeNum)+'_GridInfo.txt', gridInfo, fmt = '%.5e')

# ========================================================
# This turns a mode into an STPyFunc that vorpal can use
# ========================================================  
    def interpToStPyFunc(self,modeNum):

        mode = self.interpolateModeToYeeGrid(modeNum)

        # It seems that i cant keep the numpy.array
        # in memory when vorpal useses the stPyFunc
        # so i have to write it out and read it back 
        # in but it seems the writeout only happens
        # once so ..... i am confused
        #self.writeModeOutToTxt()
        #mode = readModeInFromTxt()

        x_Ex = mode[0][0]
        y_Ex = mode[0][1]
        mode_Ex = mode[0][2]

        x_Ey = mode[1][0]
        y_Ey = mode[1][1]
        mode_Ey = mode[1][2]
        
        # old interp script idk it is breaking
        #self.ExMode = vpInterp.BilinearInterp2DCart( x_Ex, y_Ex, mode_Ex )
        #self.EyMode = vpInterp.BilinearInterp2DCart( x_Ey, y_Ey, mode_Ey )

        self.ExMode = vpInterp.RectBivariateSpline( x_Ex, y_Ex, mode_Ex )
        self.EyMode = vpInterp.RectBivariateSpline( x_Ey, y_Ey, mode_Ey )


        return (self.ExMode, self.EyMode)








        
    

        

