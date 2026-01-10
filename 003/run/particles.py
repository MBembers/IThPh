# === IMPORTS ===
# Standard library imports
from sys import argv

# Third party imports
import numpy as np

# Local imports
import cprototype as cp
import animation as anim
from ccompiler import CSharedLibraryCompiler

# === CONSTANTS ===
if len(argv) > 1:
    try:
        # Number of particles read from command line
        NUM_OBJECTS = int(argv[1])
    except ValueError as e:
        print(e)
        NUM_OBJECTS = np.random.randint(1, 23)
        print(
            f"Number of objects is now set to randomly chosen: {NUM_OBJECTS}")
else:
    NUM_OBJECTS = 1      # Number of particles in the simulation
RADIUS = 2.0    # Initial radius for particle placement
dt = 0.01   # Timestep for the simulation

print(f"Running simulation with {NUM_OBJECTS} objects.")

# === C LIBRARY LOADING ===
# Define the path to the compiled C library (.so file)
# This assumes 'libsolver.so' is in a 'solver' directory one level *up*
# from the directory containing this Python script.
ccompiler = CSharedLibraryCompiler(source_file="../solver/solver.c")
__solver_path = ccompiler.compile()

# === INITIAL CONDITIONS ===
# CIRCLE
# positions = [_libsolver.vector(x=RADIUS * np.cos(2 * np.pi * i / NUM_OBJECTS),
#                                y=RADIUS * np.sin(2 * np.pi * i / NUM_OBJECTS))
#              for i in range(NUM_OBJECTS)]
# velocities = [_libsolver.vector(x=0, y=0) for i in range(NUM_OBJECTS)]

# DOUBLE PENDULUM
_libsolver = cp.EOMSolver(__solver_path, NUM_OBJECTS, DIMENSIONS=0, NUM_COORDS=2)
q = [_libsolver.vector(np.pi/2), _libsolver.vector(np.pi)]
dq = [_libsolver.vector(0), _libsolver.vector(0)]

# === PLOTTING SETUP ===
ani = anim.Animation2D(vector_factory=_libsolver.vector,
                        c_arr=_libsolver.c_arr,
                        next_step=_libsolver.next_step,
                        positions=q,
                        velocities=dq,
                        dt=dt,
                        DIMENSIONS=0,
                        NUM_OBJECTS=NUM_OBJECTS,
                        POINTS_PER_OBJECT=2,
                        change_coordinates=anim.Animation2D.double_pendulum_coordinates)  # Each double pendulum has 2 points
ani.create_canvas()

# === RUN ANIMATION ===
ani.run_animation()
