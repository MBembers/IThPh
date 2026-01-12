import sympy as sp
from sympy.physics.mechanics import dynamicsymbols as ds
import numpy as np
from animation import ConnectionType

# === SPRING PENDULUM CONFIGURATION ===

# === CONFIGURATION PARAMETERS ===
# Number of systems to run simultaneously
NUMBER_OF_OBJECTS = 2
# Number of points to plot per object (spring pendulum has 1 point - the mass)
POINTS_PER_OBJECT = 1
# Number of generalized coordinates per object
NUMBER_OF_COORDINATES = 2
# Timestep for the simulation
DT = 0.01  

# === PLOTTING CONFIG ===
CONNECTION_TYPE = ConnectionType.BETWEEN_POINTS_AND_ORIGIN

# === INITIAL CONDITIONS ===
# Define initial conditions in the same order as in lagrangian function
# r: radial distance (length of spring), theta: angle from vertical
INITIAL_q = [[1, np.pi+0.01], [1.2, np.pi/3]] 
INITIAL_dq = [[0.0, 0.0], [0.0, 0.0]]

# === LAGRANGIAN DEFINITION ===
# Define the Lagrangian of the system and return it along with the list of generalized coordinates
# q1, q2, ... are dynamicsymbols
# L = T - V
# return L, [q1, q2, ...]
def lagrangian():
    # --- Simple Pendulum on a Spring (r not fixed) ---
    # Coordinates
    r, theta = ds('r theta')
    r_d = r.diff()
    theta_d = theta.diff()

    # Constants
    m = 1.0        # mass
    k = 60.0       # spring constant
    r0 = 1.0       # equilibrium length of spring
    g = 9.81       # gravity

    # Kinetic Energy in polar coordinates
    # T = 0.5 * m * (r_dot^2 + r^2 * theta_dot^2)
    T = 0.5 * m * (r_d**2 + r**2 * theta_d**2)

    # Potential Energy
    # V = 0.5 * k * (r - r0)^2 - m*g*r*cos(theta)
    V = 0.5 * k * (r - r0)**2 - m * g * r * sp.cos(theta)

    L = T - V

    return L, [r, theta]

# === COORDINATE TRANSFORMATION FOR PLOTTING ===
# Function to convert generalized coordinates to Cartesian for plotting
# Polar coordinates (r, theta) -> (x, y) in Cartesian
def coordinates_transform(coordinates_list):
    coords = []
    for i in range(len(coordinates_list)):
        r = coordinates_list[i][0]
        theta = coordinates_list[i][1]
        # Convert to Cartesian coordinates
        # theta measured from vertical (downward)
        x = r * np.sin(theta)
        y = -r * np.cos(theta)
        coords.append([x, y])
    return coords
