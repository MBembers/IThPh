# === IMPORTS ===
# Standard library imports
import itertools as it

# Numpy (https://numpy.org/)
# and ctypes (https://docs.python.org/3/library/ctypes.html)
import numpy as np

# Matplotlib (https://matplotlib.org/) 
# imports for plotting and animation
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.animation as animation


class Animation2D:
    def __init__(self,vector_factory=None, c_arr=None,
                 next_step=None,positions=None, velocities=None,
                 dt=0.01, DIMENSIONS=2, NUM_OBJECTS=1, POINTS_PER_OBJECT=1, change_coordinates=None):
        self.data = np.zeros((NUM_OBJECTS * POINTS_PER_OBJECT, 2))

        self.vector = vector_factory
        self.c_arr = c_arr
        self.next_step = next_step
        self.change_coordinates = change_coordinates
        self.set_positions(positions)
        self.set_velocities(velocities)
        self.colours = it.cycle(mcolors.TABLEAU_COLORS)
        self.dt = dt
        self.DIMENSIONS = DIMENSIONS
        self.NUM_OBJECTS = NUM_OBJECTS
        self.POINTS_PER_OBJECT = POINTS_PER_OBJECT
        self.flag = True

    def set_positions(self, positions):
        self.positions = positions
                # 3. Update the data plotting array
        if self.change_coordinates is None:
            self.update_data(self.positions)
        else:
            self.update_data(self.change_coordinates(self.positions))
        

    def set_velocities(self, velocities):
        self.velocities = velocities

    def create_canvas(self,**kwargs):
        """
        This function sets up the Matplotlib figure and axes for the animation.
        It initializes the plot elements that will be updated in each frame.
        """
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal', adjustable='box')  # Ensure x and y axes have the same scale
        # Set plot limits and labels
        if 'xlim' not in kwargs:
            kwargs['xlim']=[-2.1, 2.1]
        if 'ylim' not in kwargs:
            kwargs['ylim']=[-2.1, 2.1]
        if 'xlabel' not in kwargs:
            kwargs['xlabel']='x'
        if 'ylabel' not in kwargs:
            kwargs['ylabel']='y'
        self.ax.set(**kwargs)
        
        if self.POINTS_PER_OBJECT > 1:
            # Append first point to close the loop
            self.lines = self.ax.plot(np.append(self.data[:, 0], self.data[0, 0]),
                                      np.append(self.data[:, 1], self.data[0, 1]), lw=1)[0]
        # 'points' is a scatter plot of the particles themselves
        self.points = self.ax.scatter(self.data[:, 0], self.data[:, 1],
                      c=[clr for clr, _ in zip(self.colours, range(self.NUM_OBJECTS))], s=57)
    
    def update_lines(self):
        """
        Update the lines connecting points if there are multiple points per object.
        """
        # TODO: Connect first to origin (optional) and generalize for multiple objects
            
        if self.POINTS_PER_OBJECT > 1:
            self.lines.set_data(np.append(self.data[:, 0], self.data[0, 0]),
                                np.append(self.data[:, 1], self.data[0, 1]))

    def update_frame(self, frame):
        """
        This function is called for each frame of the animation.
        It calculates the new state of the simulation and updates the plot.
        """
        global positions,velocities

        # Create empty Vector2D objects to hold the C function results
        new_positions, new_velocities = [], []
        if self.DIMENSIONS == 0 or self.DIMENSIONS == 1:
            new_positions  = [self.vector(0) for _ in range(self.NUM_OBJECTS * self.POINTS_PER_OBJECT)]
            new_velocities = [self.vector(0) for _ in range(self.NUM_OBJECTS * self.POINTS_PER_OBJECT)]
        elif self.DIMENSIONS == 2:
            new_positions  = [self.vector(x=0, y=0) for _ in range(self.NUM_OBJECTS * self.POINTS_PER_OBJECT)]
            new_velocities = [self.vector(x=0, y=0) for _ in range(self.NUM_OBJECTS * self.POINTS_PER_OBJECT)]
        elif self.DIMENSIONS == 3:
            new_positions  = [self.vector(x=0, y=0, z=0) for _ in range(self.NUM_OBJECTS * self.POINTS_PER_OBJECT)]
            new_velocities = [self.vector(x=0, y=0, z=0) for _ in range(self.NUM_OBJECTS * self.POINTS_PER_OBJECT)]

        # Convert Python lists to C arrays
        if self.flag:
            print("Python arrays:")
            print("Positions:", self.positions)
        c_positions       = self.c_arr(*self.positions)
        c_velocities      = self.c_arr(*self.velocities)
        c_new_positions   = self.c_arr(*new_positions)
        c_new_velocities  = self.c_arr(*new_velocities)

        # 1. Calculate the new positions and velocities
        self.next_step(c_positions, c_velocities,
                       c_new_positions, c_new_velocities,
                       self.dt, self.POINTS_PER_OBJECT)
        
        self.positions  = c_positions[:]
        self.velocities = c_velocities[:]
        new_positions   = c_new_positions[:]
        new_velocities  = c_new_velocities[:]

        if self.flag:
            print("Python arrays after step:")
            print("Positions:", self.positions)
            print("New positions:", new_positions)
            self.flag = False

        # 2. Update the master Python lists with the new state
        for i,new_position in enumerate(new_positions):
            self.positions[i]  = new_position
            self.velocities[i] = new_velocities[i]

        # 3. Update the data plotting array
        if self.change_coordinates is None:
            self.update_data(self.positions)
        else:
            self.update_data(self.change_coordinates(self.positions))
        

        # --- Update Matplotlib elements ---
        # Update the positions of the scattered points
        self.points.set_offsets(self.data)

        self.update_lines()

    def run_animation(self,frames=60,interval=30):
        self.ani = animation.FuncAnimation(fig=self.fig, func=self.update_frame,
                                           frames=frames , interval=interval)
        plt.show()

    def update_data(self, new_data):
        for i, coord in enumerate(new_data):
            self.data[i, :] = [coord[0], coord[1]]

    # Coordinate conversion functions:
    # TODO: add more conversion functions as needed
    # TODO: add config options to select conversion functions
    @staticmethod
    def double_pendulum_coordinates(positions):
        """
        Convert the list of angles (theta1, theta2) into (x, y) coordinates
        for plotting a double pendulum.
        Returns a list of [x, y] for each mass.
        """
        r1 = 1.0  # Length of first rod
        r2 = 1.0  # Length of second rod
        coords = []
        for i in range(len(positions)//2):
            theta1 = positions[2*i]
            theta2 = positions[2*i + 1]
            x1 = r1 * np.sin(theta1)
            y1 = -r1 * np.cos(theta1)
            x2 = x1 + r2 * np.sin(theta2)
            y2 = y1 - r2 * np.cos(theta2)
            coords.append([x1, y1])
            coords.append([x2, y2])
        return coords
