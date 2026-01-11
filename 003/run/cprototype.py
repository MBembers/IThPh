"""
A prototype Python script that demonstrates how to interface with a C library
"""

# === IMPORTS ===
from ctypes import c_float, c_int, POINTER, cdll


class EOMSolver:
    def __init__(self, path, NUM_OBJECTS=1, NUM_COORDS=0):
        """
        Load a C shared library from the specified path.
        """
        self.lib = cdll.LoadLibrary(path)
        self.NUM_OBJECTS = NUM_OBJECTS
        self.NUM_COORDS = NUM_COORDS
        self.c_vec_ptr = POINTER(c_float)  # Alias for pointer to float
        self._prototype_N()

        self.next_step.restype = None

    def _prototype_N(self):
        """
        Prototype the general N-coordinate next step function from the C library.
        Assuming function
        `void next_N(float* q, float* dq, float* new_q, float* new_dq, float dt, size_t N);`
        exists in the C library.
        (IN q, IN dq, OUT new_(q|dq), IN dt, IN N)
        """
        self.next_step = self.lib.next_N
        self.c_arr = c_float * self.NUM_COORDS # Alias for pointer to an array of floats
        self.next_step.argtypes = [self.c_vec_ptr, self.c_vec_ptr,
                                   self.c_vec_ptr, self.c_vec_ptr, c_float, c_int]

    def c_type(self, value):
        """
        Return c_float instance for the given value.
        """
        return c_float(value)