import sympy as sp
from sympy.physics.mechanics import dynamicsymbols as ds
import numpy as np
from animation import ConnectionType

# === DOUBLE PENDULUM CONFIGURATION ===

# === CONFIGURATION PARAMETERS ===
# Number of systems to run simultaneously
NUMBER_OF_OBJECTS = 2
# Number of points to plot per object (e.g., double pendulum has 2 points) 
POINTS_PER_OBJECT = 2
# Number of generalized coordinates per object
NUMBER_OF_COORDINATES = 2
# Timestep for the simulation
DT = 0.01  

# === PLOTTING CONFIG ===
CONNECTION_TYPE = ConnectionType.BETWEEN_POINTS_AND_ORIGIN

# === INITIAL CONDITIONS ===
# Define initial conditions in the same order as in lagrangian function
# e.g., for double pendulum: [[theta1, theta2]] and [[dtheta1, dtheta2]]
INITIAL_q = [[0.0, np.pi/2], [np.pi/3, np.pi]] 
INITIAL_dq = [[0.0, 0.0], [0.0, 0.0]]

# === LAGRANGIAN DEFINITION ===
# Define the Lagrangian of the system and return it along with the list of generalized coordinates
# q1, q2, ... are dynamicsymbols
# L = T - V
# return L, [q1, q2, ...]
def lagrangian():
    # --- Double Pendulum ---
    # Coordinates
    q1, q2 = ds('q1 q2')
    q1d, q2d = q1.diff(), q2.diff()

    # Constants
    m1 = 1.0
    m2 = 1.0
    l1 = 1.0
    l2 = 1.0
    g = 9.81

    # Position of mass 1
    x1 = l1 * sp.sin(q1)
    y1 = -l1 * sp.cos(q1)

    # Position of mass 2
    x2 = x1 + l2 * sp.sin(q2)
    y2 = y1 - l2 * sp.cos(q2)

    # Kinetic Energy
    v1_sq = x1.diff(ds._t)**2 + y1.diff(ds._t)**2
    v2_sq = x2.diff(ds._t)**2 + y2.diff(ds._t)**2

    T = 0.5 * m1 * v1_sq + 0.5 * m2 * v2_sq

    # Potential Energy
    V = m1*g*y1 + m2*g*y2

    L = T - V

    return L, [q1, q2]

# === COORDINATE TRANSFORMATION FOR PLOTTING ===
# Function to convert generalized coordinates to Cartesian for plotting
# coordinates_list = [[q1, q2, ...], [q1, q2, ...], ...]
# returns list of [x, y] points 
# e.g. return coords = [[x1, y1], [x2, y2], ...]
def coordinates_transform(coordinates_list):
    r1 = 1.0  # Length of first rod
    r2 = 1.0  # Length of second rod
    coords = []
    for i in range(len(coordinates_list)):
        theta1 = coordinates_list[i][0]
        theta2 = coordinates_list[i][1]
        x1 = r1 * np.sin(theta1)
        y1 = -r1 * np.cos(theta1)
        x2 = x1 + r2 * np.sin(theta2)
        y2 = y1 - r2 * np.cos(theta2)
        coords.append([x1, y1])
        coords.append([x2, y2])
    return coords
