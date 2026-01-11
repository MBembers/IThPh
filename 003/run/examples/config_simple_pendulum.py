import sympy as sp
from sympy.physics.mechanics import dynamicsymbols as ds
import numpy as np
from animation import ConnectionType

# === SIMPLE PENDULUM CONFIGURATION ===

# === CONFIGURATION PARAMETERS ===
# Number of systems to run simultaneously
NUMBER_OF_OBJECTS = 2
# Number of points to plot per object (simple pendulum has 1 point - the mass)
POINTS_PER_OBJECT = 1 
# Dimensionality of the system, 0 for non-cartesian coordinates, 
# 1 for 1D, 2 for 2D, 3 for 3D
# use 1D,2D,3D only when lagrangian is defined in cartesian coordinates
DIMENSIONS = 0
# Number of generalized coordinates per object (only for DIMENSIONS=0)
# Otherwise set to 0
NUMBER_OF_COORDINATES = 1
# Timestep for the simulation
DT = 0.01  

# === PLOTTING CONFIG ===
CONNECTION_TYPE = ConnectionType.BETWEEN_POINTS_AND_ORIGIN

# === INITIAL CONDITIONS ===
# WORKS ONLY IF DIMENSIONS = 0 OTHERWISE DEFINE MANUALLY IN PARTICLES.PY
# Define initial conditions in the same order as in lagrangian function
# e.g., for simple pendulum: [[theta]]
INITIAL_q = [[np.pi/3], [np.pi/2]] 
INITIAL_dq = [[0.0], [0.0]]

# === LAGRANGIAN DEFINITION ===
# Define the Lagrangian of the system and return it along with the list of generalized coordinates
# q1, q1, ... are dynamicsymbols
# L = T - V
# return L, [q1, q2, ...]
def lagrangian():
    # --- Simple Pendulum ---
    # Coordinates
    q = ds('q')
    qd = q.diff()

    # Constants
    m = 1.0
    l = 1.0
    g = 9.81

    # Position of mass
    x = l * sp.sin(q)
    y = -l * sp.cos(q)

    # Kinetic Energy
    v_sq = x.diff(ds._t)**2 + y.diff(ds._t)**2
    T = 0.5 * m * v_sq

    # Potential Energy
    V = m * g * y

    L = T - V

    return L, [q]

# === COORDINATE TRANSFORMATION FOR PLOTTING ===
# Needed only if DIMENSIONS = 0
# Function to convert generalized coordinates to Cartesian for plotting
# coordinates_list = [[q1, q2, ...], [q1, q2, ...], ...]
# returns list of [x, y] points 
def coordinates_transform(coordinates_list):
    r = 1.0  # Length of pendulum
    coords = []
    for i in range(len(coordinates_list)):
        theta = coordinates_list[i][0]
        x = r * np.sin(theta)
        y = -r * np.cos(theta)
        coords.append([x, y])
    return coords
