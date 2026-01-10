import sympy as sp
from sympy.physics.mechanics import dynamicsymbols as ds
import numpy as np
from animation import ConnectionType

# === CONFIGURATION PARAMETERS ===
# Number of systems to run simultaneously TODO: make this configurable
NUMBER_OF_OBJECTS = 1 
# Number of points to plot per object (e.g., double pendulum has 2 points) 
POINTS_PER_OBJECT = 2 
# Dimensionality of the system, 0 for non-cartesian coordinates, 
# 1 for 1D, 2 for 2D, 3 for 3D
# use 1D,2D,3D only when lagrangian is defined in cartesian coordinates
DIMENSIONS = 0
# Number of generalized coordinates per object (only for DIMENSIONS=0)
# Otherwise set to 0
NUMBER_OF_COORDINATES = 2
# Timestep for the simulation
DT = 0.01  

# === PLOTTING CONFIG ===
CONNECTION_TYPE = ConnectionType.BETWEEN_POINTS_AND_ORIGIN
# === INITIAL CONDITIONS ===
# WORKS ONLY IF DIMENSIONS = 0 OTHERWISE DEFINE MANUALLY IN PARTICLES.PY
# Define initial conditions in the same order as in lagrangian function
# e.g., for double pendulum: [theta1, theta2] and [dtheta1, dtheta2]
INITIAL_q = [np.pi/2, np.pi] 
INITIAL_dq = [0.0, 0.0]

# === LAGRANGIAN DEFINITION ===
# Define the Lagrangian of the system and return it along with the list of generalized coordinates
# q1, q1, ... are dynamicsymbols
# L = T - V
# return L, [q1, q2, ...]
# 
# I have no idea how to define a lagrangian using 2D or 3D vectors so the c_code generator
# correctly handles the coordinates and puts q[0].x in the code.  
def lagrangian():
    # --- Double Pendulum ---
    # Coordinates
    q1, q2 = ds('q1 q2')
    q1d, q2d = q1.diff(), q2.diff()

    # Constants
    m1, m2, l1, l2 = sp.symbols('m1 m2 l1 l2')
    g = sp.symbols('g')
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
# Needed only if DIMENSIONS = 0
# Function to convert generalized coordinates to Cartesian for plotting
# coordinates_list = [q1, q2, ...]
# returns list of [x, y] points 
# e.g. return coords = [[x1, y1], [x2, y2], ...]
def coordinates_transform(coordinates_list):
    
    r1 = 1.0  # Length of first rod
    r2 = 1.0  # Length of second rod
    coords = []
    for i in range(len(coordinates_list)//2):
        theta1 = coordinates_list[2*i]
        theta2 = coordinates_list[2*i + 1]
        x1 = r1 * np.sin(theta1)
        y1 = -r1 * np.cos(theta1)
        x2 = x1 + r2 * np.sin(theta2)
        y2 = y1 - r2 * np.cos(theta2)
        coords.append([x1, y1])
        coords.append([x2, y2])
    return coords
