# region imports
import numpy as np  # Import NumPy library and alias it as np for numerical operations and array manipulation
import math  # Import the math module for basic mathematical functions and constants
from scipy.optimize import fsolve  # Import fsolve from SciPy's optimize module to solve equations numerically
import random as rnd  # Import the random module and alias it as rnd for generating random numbers

# endregion

# region class definitions
class UC():  # a units conversion class
    def __init__(self):
        """
        This unit converter class is useful for the pipe network and perhaps other problems.
        The strategy is (number in current units)*(conversion factor)=(number desired units).
        """

    # region class constants
    ft_to_m = 1 / 3.28084 # Calculate the conversion factor from feet to meters (1 foot = 0.3048 m)
    ft2_to_m2 = ft_to_m ** 2 # Compute the conversion factor from square feet to square meters by squaring ft_to_m
    ft3_to_m3 = ft_to_m ** 3 # Compute the conversion factor from cubic feet to cubic meters by cubing ft_to_m
    ft3_to_L = ft3_to_m3 * 1000.0 # Convert cubic feet to liters (1 m³ = 1000 L, so 1 ft³ ≈ 28.3168 L)
    L_to_ft3 = 1.0 / ft3_to_L # Calculate the conversion factor from liters to cubic feet as the reciprocal of ft3_to_L
    in_to_m = ft_to_m / 12.0 # Calculate the conversion factor from inches to meters (1 inch = 1/12 foot)
    m_to_in = 1.0 / in_to_m # Calculate the conversion factor from meters to inches as the reciprocal of in_to_m
    in2_to_m2 = in_to_m ** 2 # Compute the conversion factor from square inches to square meters by squaring in_to_m
    g_SI = 9.80665 # Define gravitational acceleration in SI units (meters per second squared)
    g_EN = 32.174 # Define gravitational acceleration in English units (feet per second squared)
    gc_EN = 32.174 # Set the gravitational constant in English engineering units (lbm·ft/(lbf·s²))
    gc_SI = 1.0 # Set the gravitational constant in SI units (kg·m/(N·s²)), which is 1 by definition       # kg·m / (N·s^2)
    lbf_to_kg = 1 / 2.20462 # Compute the conversion factor from pounds (lb) to kilograms (kg)
    lbf_to_N = lbf_to_kg * g_SI # Convert pounds-force to newtons by multiplying the mass conversion factor by g_SI
    pa_to_psi = (1.0 / (lbf_to_N)) * in2_to_m2 # Compute the conversion factor from pascals to psi using lbf-to-newton and in²-to-m² conversions
    # endregion

    @classmethod
    def viscosityEnglishToSI(cls, mu, toSI=True):
        """
        Converts between lb*s/ft^2 and Pa*s
        :param mu: the viscosity in english units
        :param toSI:  True assumes english in, False assumes SI in
        :return: the viscosity in Pa*s if toSI=True, lb*s/ft^2 if toSI=False
        """
        # (lb*s)/ft^2 * (1/ft2_to_m2)*(lbf_to_kg)*g_SI => Pa*s
        cf = (1 / cls.ft2_to_m2) * (cls.lbf_to_kg) * cls.g_SI # Calculate the conversion factor using inverse ft²-to-m², lbf-to-kg, and SI gravity
        return mu * cf if toSI else mu / cf # Convert mu: multiply by cf to get SI units if toSI is True; otherwise, divide by cf for EN units

    @classmethod
    def densityEnglishToSI(cls, rho, toSI=True):
        """
        Converts between lb/ft^3 and kg/m^3
        :param rho: specific weight or density
        :param toSI:  True assumes english in, False assumes SI in
        :return: density in SI or EN
        """
        # (lb/ft^3)*(1/ft3_to_m3)*(lbf_to_kg) => kg/m^3
        cf = cls.lbf_to_kg / cls.ft3_to_m3  # Calculate density conversion factor from English to SI units (mass conversion divided by volume conversion)
        return rho * cf if toSI else rho / cf  # Convert density: multiply by cf for SI units if toSI is True; otherwise, divide by cf for English units

    @classmethod
    def psi_to_m(cls, p, rho):
        """
        For converting from psi to height of fluid in meters.
        p (psi) => Pa => h = p/(rho*g)
        """
        # first convert psi to Pa
        pa = p / cls.pa_to_psi # Convert pressure p from psi to pascals using the psi-to-pascal conversion factor
        h = pa / (rho * cls.g_SI) # Calculate pressure head h by dividing the pressure in pascals by (density * SI gravitational acceleration)
        return h # Return the computed pressure head h

    @classmethod
    def m_to_psi(cls, h, rho):
        """
        For converting from height of fluid in m to psi.
        """
        # p = rho*g*h (in Pa), then convert Pa => psi
        p_pa = rho * cls.g_SI * h # Compute pressure in pascals from density, gravitational acceleration, and head h
        return p_pa * cls.pa_to_psi # Convert the pressure from pascals to psi and return the result

class Fluid():
    def __init__(self, mu=0.00089, rho=1000, SI=True):
        """
        :param mu: dynamic viscosity (lb*s/ft^2 if SI=False, or Pa*s if SI=True)
        :param rho: density (lb/ft^3 if SI=False, or kg/m^3 if SI=True)
        :param SI: if False, convert from English to SI internally
        """
        # Convert to SI for internal usage
        if SI: # Check if the provided units are SI
            self.mu = mu # If SI, assign the dynamic viscosity directly
            self.rho = rho # If SI, assign the density directly
        else: # Otherwise, if units are not SI (assumed to be English)
            self.mu = UC.viscosityEnglishToSI(mu, toSI=True) # Convert viscosity from English to SI units using a helper function
            self.rho = UC.densityEnglishToSI(rho, toSI=True) # Convert density from English to SI units using a helper function
        self.nu = self.mu / self.rho # Calculate kinematic viscosity (nu) as dynamic viscosity divided by density (in m^2/s)

class Node(): # Define a Node class to represent a network node in a fluid system
    def __init__(self, Name='a', Pipes=None, ExtFlow=0): # Initialize a Node with a name, optional list of pipes, and external flow rate
        if Pipes is None: # Check if no pipes were provided (None)
            Pipes = [] # Set Pipes to an empty list if not provided
        self.name = Name # Assign the provided name to the node
        self.pipes = Pipes # Store the list of connected pipes for the node
        self.extFlow = ExtFlow # Set the external flow rate (in L/s; positive for inflow, negative for outflow)
        self.QNet = 0 # Initialize the net flow at the node to zero
        self.P = 0  # Initialize the pressure head at the node (in meters of fluid) to zero
        self.oCalculated = False # Flag indicating whether calculations for this node have been performed

    def getNetFlowRate(self):
        """
        Calculates the net flow rate into this node in L/s
        :return: the net flow rate into this node
        """
        Qtot = self.extFlow # Initialize the total flow (Qtot) with the node's external flow value
        for p in self.pipes: # Loop through each pipe connected to the node
            Qtot += p.getFlowIntoNode(self.name) # Add the flow from each pipe entering the node to Qtot
        self.QNet = Qtot # Set the node's net flow (QNet) to the computed total flow (Qtot)
        return self.QNet # Return the net flow value


    def setExtFlow(self, E, SI=True):
        """
        Sets the external flow rate for the node. If SI=False, E is in cfs => convert to L/s.
        """
        if SI: # Check if the system is using SI units
            self.extFlow = E # If SI, assign the external flow E directly to the node
        else: # If not SI (i.e., using English units)
            # cfs => L/s
            self.extFlow = E * UC.ft3_to_L # Convert external flow from cubic feet per second to liters per second and assign it

class Loop():
    def __init__(self, Name='A', Pipes=None): # Define the constructor for the node with a default name 'A' and optional pipes list
        if Pipes is None: # Check if no pipes list is provided
            Pipes = [] # Initialize pipes as an empty list if not provided
        self.name = Name # Set the node's name to the provided Name
        self.pipes = Pipes # Assign the provided pipes list to the node's pipes attribute

    def getLoopHeadLoss(self):
        """
        Calculates the net head loss as we traverse around the loop, in m of fluid.
        """
        deltaP = 0.0 # Initialize the total pressure drop (deltaP) to zero
        startNode = self.pipes[0].startNode # Set the starting node as the startNode of the first pipe in the list
        for p in self.pipes: # Loop over each pipe in the self.pipes list
            phl = p.getFlowHeadLoss(startNode) # Calculate the flow head loss (phl) for the current pipe starting from startNode
            deltaP += phl # Add the current pipe's head loss to the total pressure drop
            # move to the next node
            if startNode != p.endNode: # If the current startNode is not equal to the pipe's endNode
                startNode = p.endNode # then update startNode to be the pipe's endNode
            else: # Otherwise, if startNode is equal to p.endNode
                startNode = p.startNode # update startNode to be the pipe's startNode (reverse direction)
        return deltaP # Return the accumulated total pressure drop (deltaP)

class Pipe():
    def __init__(self, Start='A', End='B', L=100, D=200, r=0.00025, fluid=None, SI=True):
        """
        :param Start: the first endpoint name (string)
        :param End: the second endpoint name (string)
        :param L: pipe length (ft if SI=False, or m if SI=True)
        :param D: pipe diameter (in if SI=False, or mm if SI=True, etc. — but here we assume in if SI=False)
        :param r: pipe roughness (ft if SI=False, or m if SI=True)
        :param fluid: a Fluid() object
        :param SI: if False => interpret L in ft, D in in, r in ft => convert to internal SI
        """
        if fluid is None:
            fluid = Fluid()

        # Ensure alphabetical order for start/end node naming
        self.startNode = min(Start.lower(), End.lower())
        self.endNode = max(Start.lower(), End.lower())

        if SI:
            # already in meters
            self.length = L
            self.d = D / 1000.0
            self.rough = r
        else:
            # convert from ft to m
            self.length = L * UC.ft_to_m
            # convert from in to m
            self.d = D * UC.in_to_m
            # convert roughness from ft to m
            self.rough = r * UC.ft_to_m

        self.fluid = fluid
        # Relative roughness
        self.relrough = self.rough / self.d
        # Cross-sectional area (m^2)
        self.A = math.pi * (self.d ** 2) / 4.0

        # We store Q in L/s internally. Start with a guess of 10 L/s
        self.Q = 10.0
        # placeholders for velocity, Re, head loss
        self.vel = 0.0
        self.reynolds = 0.0
        self.hl = 0.0

    def V(self):
        """
        Velocity in m/s.  Q is in L/s => Q/1000 => m^3/s.  Divide by cross-sectional area in m^2 => m/s.
        """
        self.vel = (self.Q / 1000.0) / self.A # Compute velocity: convert flow rate Q from L/s to m³/s by dividing by 1000, then divide by cross-sectional area A
        return self.vel # Return the computed velocity

    def Re(self):
        """
        Reynolds number = (rho * V * d) / mu
        """
        v = self.V() # Retrieve velocity v by calling the method V()
        self.reynolds = (self.fluid.rho * v * self.d) / self.fluid.mu # Calculate Reynolds number: (density * velocity * diameter) divided by dynamic viscosity
        return self.reynolds # Return the computed Reynolds number

    def FrictionFactor(self):
        """
        Computes the Darcy-Weisbach friction factor using Colebrook for turbulent,
        64/Re for laminar, and a linear interpolation in the transition region.
        """
        Re = self.Re() # Get the Reynolds number by calling the method Re()
        rr = self.relrough # Retrieve the relative roughness of the pipe

        # Colebrook equation as a function
        def colebrook_eqn(f): # Define the Colebrook equation function for friction factor f (expected as an array)
            # f is an array (e.g. shape (1,))
            return 1.0 / np.sqrt(f) + 2.0 * np.log10(rr / 3.7 + 2.51 / (Re * np.sqrt(f))) # Return the residual of the Colebrook equation for f

        def colebrook(): # Define a function to compute the turbulent friction factor using the Colebrook equation
            sol = fsolve(colebrook_eqn, [0.02]) # Solve the Colebrook equation with an initial guess of 0.02 for f
            return sol[0] # Return the first element of the solution array as the friction factor

        def laminar(): # Define a function to compute the friction factor for laminar flow
            return 64.0 / Re # Return the friction factor using the laminar flow formula: 64 divided by Reynolds number

        if Re < 2000.0: # If the Reynolds number indicates laminar flow (Re < 2000)
            return laminar() # Return the laminar friction factor
        elif Re > 4000.0: # If the Reynolds number indicates turbulent flow (Re > 4000)
            return colebrook() # Return the friction factor computed via the Colebrook equation
        else: # For transitional flow (2000 <= Re <= 4000)
            # Transitional range => interpolate
            f_lam = laminar() # Calculate friction factor for laminar conditions
            f_turb = colebrook() # Calculate friction factor for turbulent conditions
            # Linear interpolation
            frac = (Re - 2000.0) / (4000.0 - 2000.0) # Determine the fractional position of Re within the transitional range
            f_mean = f_lam + frac * (f_turb - f_lam) # Interpolate linearly between laminar and turbulent friction factors
            # Return the smooth interpolation (no randomness)
            return f_mean # Return the interpolated friction factor for transitional flow

    def frictionHeadLoss(self):
        """
        Darcy-Weisbach head loss in m of fluid:
          hl = f * (L/d) * (V^2 / (2*g))
        """
        g = 9.81 # Set gravitational acceleration (m/s²)
        f = self.FrictionFactor() # Compute the friction factor using the object's method
        v = self.vel # Retrieve current velocity (ensure it is up-to-date)
        self.hl = f * (self.length / self.d) * (v ** 2 / (2.0 * g)) # Calculate head loss using friction factor, pipe geometry, and kinetic energy per unit weight
        return self.hl # Return the computed head loss

    def getFlowHeadLoss(self, s):
        """
        Signed head loss for loop traversal.  If we traverse the pipe in the same direction as
        the pipe's positive sense but the flow is negative, the sign flips, etc.
        :param s: the node we start from
        :return: signed head loss in m
        """
        nTraverse = 1 if s == self.startNode else -1 # Determine traversal direction: +1 if s is the start node, else -1
        nFlow = 1 if self.Q >= 0 else -1 # Determine flow direction: +1 for positive flow, -1 for negative flow
        # friction head loss is always positive magnitude, but we attach sign:
        return nTraverse * nFlow * self.frictionHeadLoss() # Return friction head loss with sign adjusted by traversal and flow directions

    def Name(self):
        """
        e.g. 'a-b'
        """
        return self.startNode + '-' + self.endNode # Return a string identifier combining the start and end node names with a hyphen

    def oContainsNode(self, node):
        return (self.startNode == node) or (self.endNode == node) # Return True if the given node matches either the start or end node; otherwise, False

    def printPipeFlowRate(self, SI=True):
        """
        Prints line like: The flow in segment a-b is 3.57 (cfs) and Re=286475.8
        """
        q_units = 'L/s' if SI else 'cfs' # Set q_units to 'L/s' if SI is True, otherwise 'cfs'
        if SI: # Check if using SI units
            q_val = self.Q # Use the flow value Q directly in SI units (L/s)
        else: # Otherwise, if using non-SI units
            # convert L/s => cfs
            q_val = self.Q * UC.L_to_ft3 # Convert flow from liters per second (L/s) to cubic feet per second (cfs) using the conversion factor

        print('The flow in segment {} is {:0.2f} ({}) and Re={:.1f}'.format(
            self.Name(), q_val, q_units, self.reynolds))

    def printPipeHeadLoss(self, SI=True):
        """
        Print line like: head loss in pipe a-b (L=1000.00 in, d=18.00 in) is 12.22 in of water
        matching the exact problem statement.
        """
        # For diameter and length printing:
        # If SI=False => length in inches, diameter in inches
        if SI: # Check if SI units are being used
            cfL = 1.0 # Set conversion factor for length (SI: no conversion needed)
            unitsL = 'm' # Set length unit label to meters
            cfd = 1000.0 # Set conversion factor for diameter from meters to millimeters (1 m = 1000 mm)
            unitsD = 'mm' # Set diameter unit label to millimeters
            # head loss in m of water
            valHL = self.hl # Use the head loss value directly in SI (meters of water)
            unitsHL = 'm of water' # Set head loss unit label to "m of water"
            L_print = self.length * cfL # Convert the pipe length for display (remains in meters)
            D_print = self.d * cfd # Convert the pipe diameter from meters to millimeters for display
        else: # If not using SI units (using English units)
            # we want to show length in inches, diameter in inches
            # 1 m = 39.3701 in
            cf_m_to_in = UC.m_to_in # Get conversion factor from meters to inches from the UC module
            L_print = self.length * cf_m_to_in # Convert the pipe length from meters to inches
            D_print = self.d * cf_m_to_in # Convert the pipe diameter from meters to inches
            unitsL = 'in' # Set length unit label to inches
            unitsD = 'in' # Set diameter unit label to inches
            # head loss in inches of water
            valHL = self.hl * UC.m_to_in # Convert head loss from meters to inches using the conversion factor
            unitsHL = 'in of water' # Set head loss unit label to "in of water"

        print("head loss in pipe {} (L={:.2f} {}, d={:.2f} {}) is {:.2f} {}".format(
            self.Name(), L_print, unitsL, D_print, unitsD, valHL, unitsHL)) # Print formatted output with pipe name, length, diameter, head loss, and their units

    def getFlowIntoNode(self, n):
        """
        returns +Q if flow enters node n, -Q if flow leaves node n
        """
        if n == self.startNode: # Check if the provided node n is the start node of the pipe
            return -self.Q # If it is, return the negative flow Q (indicating flow leaving the start node)
        else: # Otherwise, if n is not the start node
            return self.Q # Return the positive flow Q (indicating flow entering the node)

class PipeNetwork(): # Define a class to represent a network of pipes, loops, and nodes
    def __init__(self, Pipes=None, Loops=None, Nodes=None, fluid=None): # Constructor with optional parameters for pipes, loops, nodes, and fluid
        if Pipes is None: # Check if no pipes list is provided
            Pipes = [] # Initialize Pipes as an empty list
        if Loops is None: # Check if no loops list is provided
            Loops = [] # Initialize Loops as an empty list
        if Nodes is None: # Check if no nodes list is provided
            Nodes = [] # Initialize Nodes as an empty list
        if fluid is None: # Check if no fluid object is provided
            fluid = Fluid() # Create a default Fluid object (assuming Fluid() is defined elsewhere)
        self.pipes = Pipes # Assign the pipes list to the instance variable pipes
        self.loops = Loops # Assign the loops list to the instance variable loops
        self.nodes = Nodes # Assign the nodes list to the instance variable nodes
        self.Fluid = fluid # Assign the fluid object to the instance variable Fluid

    def getPipe(self, pipeName): # Define a method to retrieve a pipe by its name
        for p in self.pipes: # Iterate over each pipe in the network's pipes list
            if p.Name() == pipeName: # Check if the current pipe's name matches the provided pipeName
                return p # Return the matching pipe if found
        return None # If no matching pipe is found, return None

    def buildNodes(self):
        """
        Automatically create the node objects by looking at the pipe endpoints.
        """
        for p in self.pipes: # Iterate over each pipe in the network's list of pipes
            if not self.nodeBuilt(p.startNode): # If the start node of the pipe has not been created yet
                self.nodes.append(Node(p.startNode, self.getNodePipes(p.startNode))) # Create a new Node for the start node with its associated pipes and add it to the nodes list
            if not self.nodeBuilt(p.endNode): # If the end node of the pipe has not been created yet
                self.nodes.append(Node(p.endNode, self.getNodePipes(p.endNode))) # Create a new Node for the end node with its associated pipes and add it to the nodes list

    def nodeBuilt(self, nodeName): # Define a method to check if a node with the given name already exists
        for n in self.nodes: # Loop over each existing node in the network
            if n.name == nodeName: # Check if the current node's name matches the provided nodeName
                return True # Return True if a matching node is found
        return False # Return False if no node with the provided name exists

    def getNodePipes(self, node): # Define a method to get all pipes connected to a given node
        return [p for p in self.pipes if p.oContainsNode(node)] # Return a list of pipes that include the specified node using a list comprehension

    def getNode(self, name): # Define a method to retrieve a node by its name
        for n in self.nodes: # Loop over each node in the network
            if n.name == name: # Check if the current node's name matches the provided name
                return n # Return the node if a match is found
        return None # Return None if no node with the specified name exists

    def findFlowRates(self):
        """
        Use fsolve to find flows that make net node flow = 0 and net loop head loss = 0.
        """
        # number of equations = number of nodes + number of loops
        N = len(self.nodes) + len(self.loops)
        # initial guess for fsolve
        Q0 = np.full(N, 10.0)

        def fn(q):
            # put flows into each pipe
            for i in range(len(self.pipes)):
                self.pipes[i].Q = q[i]  # in L/s

            # Now build up the residuals: node mass balances + loop head losses
            # node flow rates
            res = self.getNodeFlowRates()
            # loop head losses
            res += self.getLoopHeadLosses()
            return res

        # solve
        sol = fsolve(fn, Q0)
        return sol

    def getNodeFlowRates(self):
        """
        For each node, compute net flow in.  We want that = 0 in the solution.
        """
        qNet = [] # Initialize an empty list to store net flow rates from each node
        for n in self.nodes: # Iterate over each node in the network's nodes list
            qNet.append(n.getNetFlowRate()) # Append the net flow rate (in L/s) of the current node to qNet
        return qNet # Return the list of net flow rates

    def getLoopHeadLosses(self):
        """
        For each loop, compute net head loss in m of fluid.  We want that = 0.
        """
        lhl = [] # Initialize an empty list to store head losses from each loop
        for l in self.loops: # Iterate over each loop in the network's loops list
            lhl.append(l.getLoopHeadLoss()) # Append the head loss of the current loop to lhl
        return lhl # Return the list of loop head losses

    def getNodePressures(self, knownNodeP, knownNode):
        """
        Once flows are known, we compute node 'head' (in m of fluid).
        Then we shift all so that knownNode has the specified pressure in psi.
        """
        # reset all
        for n in self.nodes: # Loop over each node in the network
            n.P = 0.0 # Reset the node's pressure head (P) to 0.0
            n.oCalculated = False # Mark the node as not yet calculated (oCalculated flag set to False)

        # Do a loop traversal to set node heads
        for l in self.loops:  # Loop over each loop in the network
            startNode = l.pipes[0].startNode # Set startNode as the start node of the first pipe in the loop
            nodeObj = self.getNode(startNode) # Retrieve the node object corresponding to startNode
            nodeObj.oCalculated = True # Mark this node as calculated (set oCalculated flag to True)
            CurrentP = nodeObj.P # Initialize the current pressure (CurrentP) with the pressure at startNode
            for p in l.pipes: # Loop over each pipe in the current loop
                phl = p.getFlowHeadLoss(startNode) # Compute the flow head loss (phl) for pipe p starting at startNode
                CurrentP -= phl # Decrease the current pressure by the flow head loss
                if startNode != p.endNode: # If startNode is not the end node of pipe p
                    startNode = p.endNode # Update startNode to be the end node of pipe p
                else: # Otherwise, if startNode is the end node
                    startNode = p.startNode # Update startNode to be the start node of pipe p (switch direction)
                self.getNode(startNode).P = CurrentP # Update the pressure head (P) of the new startNode with the updated CurrentP

        # shift so that knownNode is at knownNodeP (psi)
        kn = self.getNode(knownNode)
        # kn.P is in m of fluid; we want it to become knownNodeP (psi)
        # difference in psi
        desired_m = UC.psi_to_m(knownNodeP, self.Fluid.rho)  # Convert knownNodeP (in psi) to meters of water head using the fluid's density
        deltaP = desired_m - kn.P # Compute the pressure difference (deltaP) between the desired head and the known node's current pressure
        for n in self.nodes: # Iterate over each node in the network
            n.P += deltaP # Adjust each node's pressure by adding the computed pressure difference

    def printPipeFlowRates(self, SI=True): # Define a method to print flow rates for each pipe with an option for SI or English units
        for p in self.pipes: # Iterate over each pipe in the network
            p.printPipeFlowRate(SI=SI) # Call the pipe's method to print its flow rate using the specified unit system

    def printNetNodeFlows(self, SI=True): # Define a method to print net flow rates for each node with an option for SI or English units
        for n in self.nodes: # Iterate over each node in the network
            if SI: # If using SI units
                qval = n.QNet # Use the node's net flow value directly (in L/s)
                units = 'L/S' # Set the unit label to L/s
            else: # If not using SI units (i.e., using English units)
                qval = n.QNet * UC.L_to_ft3 # Convert the net flow from liters per second to cubic feet per second (cfs)
                units = 'cfs' # Set the unit label to cfs
            print('net flow into node {} is {:0.2f} ({})'.format(n.name, qval, units)) # Print the node's name and its net flow with units

    def printLoopHeadLoss(self, SI=True):
        """
        e.g.: head loss for loop A is 0.00 (psi)
        """
        for l in self.loops: # Iterate over each loop in the network
            val = l.getLoopHeadLoss() # Calculate the head loss for the current loop
            if SI: # Check if SI units should be used
                # in m => just print
                print('head loss for loop {} is {:0.2f} (m of water)'.format(l.name, val)) # Print the loop's head loss in meters of water
            else: # Otherwise, if non-SI units are used
                # convert m => psi
                psi_val = UC.m_to_psi(val, self.Fluid.rho) # Convert the head loss from meters to psi using the fluid's density
                print('head loss for loop {} is {:0.2f} (psi)'.format(l.name, psi_val)) # Print the loop's head loss in psi

    def printPipeHeadLoss(self, SI=True): # Define a method to print the head loss for each pipe, with an option for SI units
        for p in self.pipes: # Iterate over each pipe in the network
            p.printPipeHeadLoss(SI=SI) # Call each pipe's method to print its head loss in the specified unit system

    def printNodePressures(self, SI=True):
        """
        Print final pressures at each node in the requested units
        """
        for n in self.nodes: # Iterate over each node in the network
            if SI: # Check if SI units should be used
                # m of water
                print('Pressure at node {} = {:0.2f} m of water'.format(n.name, n.P)) # Print the node's pressure in meters of water
            else: # Otherwise, if non-SI units are used
                # convert m => psi
                psi_val = UC.m_to_psi(n.P, self.Fluid.rho) # Convert the node's pressure from meters to psi using the fluid's density
                print('Pressure at node {} = {:0.2f} psi'.format(n.name, psi_val)) # Print the node's pressure in psi

# endregion

# region function definitions
def main():
    """
    This program analyzes the given pipe network in English units (SI=False).
    We fill out all pipes, nodes, loops, external flows, then solve with fsolve.
    Finally, we print the results in exactly the format and numbers required.
    """
    SIUnits = False
    # Room temp water: mu=20.50e-6 lb*s/ft^2, gamma=62.3 lb/ft^3
    water = Fluid(mu=20.50e-6, rho=62.3, SI=False)

    # Roughness: 12" and 16" => cast iron => 0.00085 ft
    #            18" and 24" => concrete => 0.003 ft
    r_CI = 0.00085  # ft (cast iron)
    r_CN = 0.003    # ft (concrete)

    PN = PipeNetwork(fluid=water)

    # -----------------------------
    #  Add all the pipes to PN
    #  L in ft, D in inches, roughness in ft
    # -----------------------------
    # a-b: 1000 ft, 18 in, concrete
    PN.pipes.append(Pipe('a', 'b', L=1000, D=18, r=r_CN, fluid=water, SI=SIUnits))
    # a-h: 1600 ft, 24 in, concrete
    PN.pipes.append(Pipe('a', 'h', L=1600, D=24, r=r_CN, fluid=water, SI=SIUnits))
    # b-c: 500 ft, 18 in, concrete
    PN.pipes.append(Pipe('b', 'c', L=500, D=18, r=r_CN, fluid=water, SI=SIUnits))
    # b-e: 800 ft, 16 in, cast iron
    PN.pipes.append(Pipe('b', 'e', L=800, D=16, r=r_CI, fluid=water, SI=SIUnits))
    # c-d: 500 ft, 18 in, concrete
    PN.pipes.append(Pipe('c', 'd', L=500, D=18, r=r_CN, fluid=water, SI=SIUnits))
    # c-f: 800 ft, 16 in, cast iron
    PN.pipes.append(Pipe('c', 'f', L=800, D=16, r=r_CI, fluid=water, SI=SIUnits))
    # d-g: 800 ft, 16 in, cast iron
    PN.pipes.append(Pipe('d', 'g', L=800, D=16, r=r_CI, fluid=water, SI=SIUnits))
    # e-f: 500 ft, 12 in, cast iron
    PN.pipes.append(Pipe('e', 'f', L=500, D=12, r=r_CI, fluid=water, SI=SIUnits))
    # e-i: 800 ft, 18 in, concrete
    PN.pipes.append(Pipe('e', 'i', L=800, D=18, r=r_CN, fluid=water, SI=SIUnits))
    # f-g: 500 ft, 12 in, cast iron
    PN.pipes.append(Pipe('f', 'g', L=500, D=12, r=r_CI, fluid=water, SI=SIUnits))
    # g-j: 800 ft, 18 in, concrete
    PN.pipes.append(Pipe('g', 'j', L=800, D=18, r=r_CN, fluid=water, SI=SIUnits))
    # h-i: 1000 ft, 24 in, concrete
    PN.pipes.append(Pipe('h', 'i', L=1000, D=24, r=r_CN, fluid=water, SI=SIUnits))
    # i-j: 1000 ft, 24 in, concrete
    PN.pipes.append(Pipe('i', 'j', L=1000, D=24, r=r_CN, fluid=water, SI=SIUnits))

    # Build node objects automatically
    PN.buildNodes()

    # External flows:  +10 cfs at h, -3 cfs at e, -5 cfs at f, -2 cfs at d
    PN.getNode('h').setExtFlow(10, SI=SIUnits)
    PN.getNode('e').setExtFlow(-3, SI=SIUnits)
    PN.getNode('f').setExtFlow(-5, SI=SIUnits)
    PN.getNode('d').setExtFlow(-2, SI=SIUnits)

    # -----------------------------
    # Add loops
    # (We define 4 loops to match the final check of "loop A,B,C,D")
    # -----------------------------
    # Loop A: a->b->e->i->h->a
    PN.loops.append(
        Loop('A', [
            PN.getPipe('a-b'),
            PN.getPipe('b-e'),
            PN.getPipe('e-i'),
            PN.getPipe('h-i'),
            PN.getPipe('a-h')
        ])
    )
    # Loop B: b->c->f->e->b
    PN.loops.append(
        Loop('B', [
            PN.getPipe('b-c'),
            PN.getPipe('c-f'),
            PN.getPipe('e-f'),
            PN.getPipe('b-e')
        ])
    )
    # Loop C: c->d->g->f->c
    PN.loops.append(
        Loop('C', [
            PN.getPipe('c-d'),
            PN.getPipe('d-g'),
            PN.getPipe('f-g'),
            PN.getPipe('c-f')
        ])
    )
    # Loop D: i->j->g->f->e->i  (one possible traversal)
    PN.loops.append(
        Loop('D', [
            PN.getPipe('i-j'),
            PN.getPipe('g-j'),
            PN.getPipe('f-g'),
            PN.getPipe('e-f'),
            PN.getPipe('e-i')
        ])
    )

    # Solve for flows
    PN.findFlowRates()

    # Set node pressures so that node h = 80 psi
    PN.getNodePressures(knownNode='h', knownNodeP=80.0)

    # Print final results EXACTLY as in the problem statement
    PN.printPipeFlowRates(SI=SIUnits) # Call method to print each pipe's flow rate using SI units if SIUnits is True
    print() # Print a blank line for spacing
    print('Check node flows:') # Print a header message indicating that node flows will be checked
    PN.printNetNodeFlows(SI=SIUnits) # Call method to print net flow rates at each node using the SIUnits flag for unit selection
    print() # Print a blank line for spacing
    print('Check loop head loss:') # Print a header message indicating that loop head loss will be checked
    PN.printLoopHeadLoss(SI=SIUnits) # Call method to print the head loss for each loop using SI units if SIUnits is True
    print() # Print a blank line for spacing
    PN.printPipeHeadLoss(SI=SIUnits) # Call method to print the head loss for each pipe using the SIUnits flag
    print() # Print a blank line for spacing
    PN.printNodePressures(SI=SIUnits) # Call method to print the pressure at each node using SI units if SIUnits is True

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion