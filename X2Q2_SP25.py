# region imports
from scipy.integrate import solve_ivp # Import solve_ivp from SciPy for solving initial value problems for ODEs
from math import sin  # Import the sine function from the math module for computing trigonometric sine
import numpy as np  # Import NumPy for efficient numerical operations and array handling
import matplotlib.pyplot as plt # Import pyplot from Matplotlib for creating plots and visualizations

# endregion

class Circuit: # Define a Circuit class to model an electrical circuit
    def __init__(self, R=10, L=20, C=0.05, A=20, w=20, p=0): # Initialize a Circuit instance with default parameters: resistance (R), inductance (L), capacitance (C), amplitude (A), angular frequency (w), and phase (p)
        """
        Initialize the RLC circuit parameters.
        :param R: Resistance (Ohms)
        :param L: Inductance (Henry)
        :param C: Capacitance (Farads)
        :param A: Amplitude of input voltage (Volts)
        :param w: Frequency of input voltage (rad/s)
        :param p: Phase of input voltage (radians)
        """
        self.R = R # assign the passed resistance value to the instance variable R
        self.L = L # assign the passed inductance value to the instance variable L
        self.C = C # assign the passed capacitance value to the instance variable C
        self.A = A # assign the passed amplitude value to the instance variable A
        self.w = w  # assign the passed angular frequency value to the instance variable w
        self.p = p  # assign the passed phase value to the instance variable p
        self.solution = None  # initialize the solution attribute to None (to be computed later)
        self.t_vals = None # initialize the time values attribute to None (to be computed later)

    def ode_system(self, t, X):
        """
        Defines the system of differential equations.
        :param t: Time variable
        :param X: State variables [i1, i2]
        :return: Derivatives [di1/dt, di2/dt]
        """
        i1, i2 = X  # Current through L and R
        v_t = self.A * sin(self.w * t + self.p)  # Input voltage

        di1_dt = (v_t - self.R * i2 - i1 * self.L) / self.L # Compute the time derivative of i1 based on voltage, resistance, and inductance
        di2_dt = (i1 - i2) / self.C  # Compute the time derivative of i2 as the difference between i1 and i2 divided by the capacitance

        return [di1_dt, di2_dt]  # Return a list containing the derivatives [di1_dt, di2_dt]

    def simulate(self, t=10, pts=500):
        """
        Simulate the circuit response over a given time.
        :param t: Total simulation time (seconds)
        :param pts: Number of time points
        """
        t_span = (0, t)  # Define the time interval for integration from 0 to t
        t_eval = np.linspace(0, t, pts) # Generate an array of pts evenly spaced time points between 0 and t

        X0 = [0, 0]  # Initial conditions for i1 and i2

        self.solution = solve_ivp(self.ode_system, t_span, X0, t_eval=t_eval)  # Solve the ODE system using the specified time span, initial conditions, and evaluation points, then store the solution
        self.t_vals = t_eval # Save the array of time points to the object's attribute for later use

    def doPlot(self):
        """
        Generate the plot for i1, i2, and vc over time.
        """
        if self.solution is None: # Check if simulation data is available
            print("No simulation data. Run simulate() first.") # Inform the user that simulation data is missing
            return  # Exit the function early since there's no data to process

        t = self.t_vals  # Retrieve the simulation time values
        i1 = self.solution.y[0] # Extract the first current (i1) from the simulation results
        i2 = self.solution.y[1] # Extract the second current (i2) from the simulation results
        vc = (i1 - i2) * self.R   # Calculate the voltage across the capacitor based on the current difference and resistance

        fig, ax1 = plt.subplots() # Create a new figure and primary axis for plotting
        ax2 = ax1.twinx() # Create a secondary y-axis that shares the same x-axis as ax1

        ax1.plot(t, i1, label="$i_1(t)$ (A)", linestyle='solid')  # Plot i1 vs. time on the primary axis with a solid line
        ax1.plot(t, i2, label="$i_2(t)$ (A)", linestyle='dashed') # Plot i2 vs. time on the primary axis with a dashed line
        ax2.plot(t, vc, label="$v_c(t)$ (V)", linestyle='dotted', color='gray') # Plot capacitor voltage vs. time on the secondary axis with a dotted gray line

        ax1.set_xlabel("Time (s)")  # Set the label for the x-axis
        ax1.set_ylabel("Current (A)") # Set the label for the primary y-axis (current)
        ax2.set_ylabel("Voltage (V)")  # Set the label for the secondary y-axis (voltage)

        ax1.legend()  # Add a legend to the primary axis for current plots
        plt.title("RLC Circuit Response") # Set the overall title of the plot
        plt.show()  # Display the plot


def main():
    """
    Run the circuit simulation with user input.
    """
    while True:  # Start an infinite loop to allow repeated simulations until the user opts out
        try:  # Begin a try block to catch input conversion errors
            R = float(input("Enter resistance (Ohms, default 10): ") or 10)  # Prompt for resistance input; use 10 if input is empty
            L = float(input("Enter inductance (H, default 20): ") or 20) # Prompt for inductance input; use 20 if input is empty
            C = float(input("Enter capacitance (F, default 0.05): ") or 0.05) # Prompt for capacitance input; use 0.05 if input is empty
            A = float(input("Enter voltage amplitude (V, default 20): ") or 20) # Prompt for voltage amplitude; default to 20 if empty
            w = float(input("Enter voltage frequency (rad/s, default 20): ") or 20) # Prompt for voltage frequency; default to 20 if empty
            p = float(input("Enter voltage phase (radians, default 0): ") or 0) # Prompt for voltage phase; default to 0 if empty

            circuit = Circuit(R, L, C, A, w, p)  # Create a Circuit object with the user-provided parameters
            circuit.simulate() # Run the simulation for the created circuit
            circuit.doPlot()  # Generate and display a plot of the simulation results

            repeat = input("Run again? (y/n): ") # Ask the user if they want to run the simulation again
            if repeat.lower() != 'y':  # Check if the user's response is not 'y' (case-insensitive)
                break # Exit the loop if the user does not want to continue
        except ValueError: # Catch exceptions raised by invalid (non-numeric) input
            print("Invalid input. Please enter numeric values.")  # Inform the user of the invalid input error
            continue # Restart the loop to prompt the user for input again


if __name__ == "__main__":
    main()