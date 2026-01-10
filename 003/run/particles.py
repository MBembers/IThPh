# === IMPORTS ===
# Standard library imports
from sys import argv

# Third party imports
import numpy as np

# Local imports
import cprototype as cp
import animation as anim
from ccompiler import CSharedLibraryCompiler
from lagrangian import LagrangianToC
import config

# === GENERATING C CODE FROM LAGRANGIAN ===
LagrangianToC.config_generate()

# === C LIBRARY LOADING ===
# Define the path to the compiled C library (.so file)
# This assumes 'libsolver.so' is in a 'solver' directory one level *up*
# from the directory containing this Python script.
ccompiler = CSharedLibraryCompiler(source_file="../solver/solver.c")
__solver_path = ccompiler.compile()

# === INITIAL CONDITIONS ===
# DOUBLE PENDULUM
_libsolver = cp.EOMSolver(__solver_path, NUM_OBJECTS=config.NUMBER_OF_OBJECTS, 
                          DIMENSIONS=config.DIMENSIONS, NUM_COORDS=config.NUMBER_OF_COORDINATES)
# AUTOMATIC GENERATION ONLY FOR GENERALIZED COORDINATES!!!
q = []
dq = []
for i in range(config.NUMBER_OF_COORDINATES * config.NUMBER_OF_OBJECTS):
    q.append(_libsolver.vector(config.INITIAL_q[i]))
    dq.append(_libsolver.vector(config.INITIAL_dq[i]))

# === PLOTTING SETUP ===
ani = anim.Animation2D(vector_factory=_libsolver.vector,
                        c_arr=_libsolver.c_arr,
                        next_step=_libsolver.next_step,
                        positions=q,
                        velocities=dq,
                        dt=config.DT,
                        DIMENSIONS=config.DIMENSIONS,
                        NUM_OBJECTS=config.NUMBER_OF_OBJECTS,
                        POINTS_PER_OBJECT=config.POINTS_PER_OBJECT,
                        change_coordinates=config.coordinates_transform)  # Each double pendulum has 2 points
ani.create_canvas()

# === RUN ANIMATION ===
print(f"Running simulation with {config.NUMBER_OF_OBJECTS} \
      objects and {ani.POINTS_PER_OBJECT * config.NUMBER_OF_OBJECTS} points.")
ani.run_animation()
