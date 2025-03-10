# region imports
from scipy.integrate import solve_ivp
from math import sin
import math
import numpy as np
from matplotlib import pyplot as plt
# endregion
#region class definitions
class circuit():
    def __init__(self, R=20, L=20, C=0.05, A=20, w=20, p=0):
        '''
        #JES MISSING DOCSTRING
        :param R: #JES MISSING argument descriptions (give units)
        :param L:
        :param C:
        :param A:
        :param w:
        :param p:
        '''
        #region attributes
        #endregion

    #region methods
    def ode_system(self, t, X):
        """
        this is the odeSystem callback I'm using for solve_ivp().
        :param X: the current values of the state variables
        :param t: the current time
        :return: list of derivatives of state variables
        """
        pass

    def simulate(self, t=10, pts=500):
        """
        For simulating transient behavior of circuit.
        :param: time over which to carry out the simulation in seconds
        :param pts: number of points in the simulation
        :return: nothing, just store I
        """
        pass

    def doPlot(self, ax=None):
        """
        Re-written on 4/21/2022 to adapt to plotting on GUI if ax is not None
        :param args: contains ((R, list of time values, and results of solve_ivp))
        :param ax:
        :return:
        """
        if ax == None:
            ax = plt.subplot()
            QTPlotting = False  # actually, we are just using CLI and showing the plot
        else:
            QTPlotting = True

        #JES MISSING CODE for making the plot

        if not QTPlotting:
            plt.show()
    #endregion

#endregion
# region function definitions

def main():
    """
    For solving problem 2 on exam.
    :return:
    """
    goAgain = True
    Circuit = circuit(R=20, L=20, C=0.05, A=20,w=20,p=0)  # create a circuit object with default values
    while goAgain:
        #JES MISSING CODE for soliciting user input.
        Circuit.simulate(t=10, pts=500)
        Circuit.doPlot()
        #JES MISSING CODE for soliciting user input.
    pass
# endregion

# region function calls
if __name__ ==  "__main__":
    main()
# endregion