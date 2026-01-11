import sympy as sp
from sympy.physics.mechanics import dynamicsymbols as ds
import numpy as np
from animation import ConnectionType

# === SPHERICAL PENDULUM CONFIGURATION ===

# === CONFIGURATION PARAMETERS ===
# Number of systems to run simultaneously
NUMBER_OF_OBJECTS = 1
# Number of points to plot per object (spherical pendulum has 1 point - the mass)
POINTS_PER_OBJECT = 1 
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
# theta: polar angle from vertical, phi: azimuthal angle
INITIAL_q = [[np.pi/4, 0.0], [np.pi/6, np.pi/2]] 
INITIAL_dq = [[0.0, 0.0], [0.0, 0.0]]

# === LAGRANGIAN DEFINITION ===
# Define the Lagrangian of the system and return it along with the list of generalized coordinates
# q1, q1, ... are dynamicsymbols
# L = T - V
# return L, [q1, q2, ...]
def lagrangian():
    # --- Spherical Pendulum (viewed from above) ---
    # Coordinates
    theta, phi = ds('theta phi')
    theta_d = theta.diff()
    phi_d = phi.diff()

    # Constants
    m = 1.0
    l = 1.0
    g = 9.81

    # Kinetic Energy in spherical coordinates
    # T = 0.5 * m * l^2 * (theta_dot^2 + sin^2(theta) * phi_dot^2)
    T = 0.5 * m * l**2 * (theta_d**2 + sp.sin(theta)**2 * phi_d**2)

    # Potential Energy (height measured from pivot point)
    # V = -m*g*l*cos(theta)
    V = -m * g * l * sp.cos(theta)

    L = T - V

    return L, [theta, phi]

# === COORDINATE TRANSFORMATION FOR PLOTTING ===
# Needed only if DIMENSIONS = 0
# Function to convert generalized coordinates to Cartesian for plotting
# Spherical coordinates (theta, phi) -> (x, y) projection on horizontal plane
def coordinates_transform(coordinates_list):
    r = 1.0  # Length of pendulum
    coords = []
    for i in range(len(coordinates_list)):
        theta = coordinates_list[i][0]
        phi = coordinates_list[i][1]
        # Project onto horizontal plane (top-down view)
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        coords.append([x, y])
    return coords
