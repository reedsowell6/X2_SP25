# region imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad


# endregion

# region function definitions
def S(x):
    """
    Computes S(x) = ∫(from 0 to x) sin(t^2) dt using quad.
    """

    # integrand for S(x)
    def integrand(t):
        return np.sin(t ** 2)

    # quad returns (value, error_estimate); we just want the value
    s = quad(integrand, 0, x)
    return s[0]


def Exact(x):
    """
    Returns the exact solution y = 1/(2.5 - S(x)) + 0.01*x^2
    """
    return 1.0 / (2.5 - S(x)) + 0.01 * x ** 2


def ODE_System(x, y):
    """
    ODE system for y' = (y - 0.01x^2)^2 * sin(x^2) + 0.02x
    Here y is a list/array but we have only one state variable y[0].
    """
    Y = y[0]
    Ydot = (Y - 0.01 * x ** 2) ** 2 * np.sin(x ** 2) + 0.02 * x
    return [Ydot]


def PlotResults(*args):
    """
    Produce the plot according to the formatting criteria:
      - Exact solution: solid line
      - Numerical solution: upward facing triangles
      - x from 0.0 to 6.0, y from 0.0 to 1.0
      - tick marks inward, legend, one decimal format, etc.
    """
    xRange_Num, y_Num, xRange_Xct, y_Xct = args

    fig, ax = plt.subplots()

    # Plot exact (solid line)
    ax.plot(xRange_Xct, y_Xct, '-', label='Exact')

    # Plot numerical (upward triangles)
    ax.plot(xRange_Num, y_Num, '^', label='Numerical')

    # Set axis limits
    ax.set_xlim([0.0, 6.0])
    ax.set_ylim([0.0, 1.0])

    # Ticks pointed inward
    ax.tick_params(axis='both', which='both', direction='in')

    # Axis labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    # Format axis numbers with one decimal place
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))

    # Legend
    ax.legend()

    # Title
    ax.set_title("IVP: y'=(y-0.01x^2)^2 sin(x^2)+0.02x, y(0)=0.4")

    plt.show()


def main():
    """
    Solve the IVP:
      y' = (y - 0.01x^2)**2 * sin(x^2) + 0.02x
      y(0)=0.4
    over 0 <= x <= 5 with step size h=0.2,
    then plot the numerical and exact solutions.
    """
    # x values for the numerical solution (step size 0.2)
    xRange = np.arange(0, 5.2, 0.2)  # goes a bit past 5 so it includes 5.0
    # x values for the exact solution (fine grid up to x=6 to match plot range)
    xRange_xct = np.linspace(0, 6, 600)

    # Initial condition
    Y0 = [0.4]

    # Solve IVP numerically
    sln = solve_ivp(ODE_System, [0, 5], Y0, t_eval=xRange)

    # Compute exact solution
    xctSln = np.array([Exact(x) for x in xRange_xct])

    # Plot results
    PlotResults(xRange, sln.y[0], xRange_xct, xctSln)


# end region

# region function calls
if __name__ == "__main__":
    main()
# end region
