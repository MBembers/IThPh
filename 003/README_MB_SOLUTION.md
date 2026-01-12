# Finished project from Introduction to Theoretical Physics

As a sidenote:
I probably changed a little too much, but the project was very enjoyable to work on. In its current state, the program is far from being optimized in terms of performance because I focused on improving usability and user-friendliness.

## Setup

### 1. Create a Python Virtual Environment

Navigate to the project directory and create a virtual environment:

```bash
cd 003
python3 -m venv venv
```

### 2. Activate the Virtual Environment

On Linux/Mac:
```bash
source venv/bin/activate
```

On Windows:
```bash
venv\Scripts\activate
```

### 3. Install Dependencies

Install all required packages from the requirements file:

```bash
pip install -r run/requirements.txt
```

## Configuration

Edit the `run/config.py` file to configure your simulation with your own desired configuration:

- **System parameters**: Modify `NUMBER_OF_OBJECTS`, `POINTS_PER_OBJECT`, `NUMBER_OF_COORDINATES`
- **Simulation timestep**: Adjust `DT` for accuracy vs. performance
- **Initial conditions**: Set `INITIAL_q` and `INITIAL_dq` arrays
- **Lagrangian function**: Define your system's Lagrangian in the `lagrangian()` function
- **Coordinate transformation**: Implement `coordinates_transform()` to convert generalized coordinates to Cartesian for visualization

## Running the Simulation

From the `run/` directory, execute:

```bash
cd run
python3 -m particles
```

## Examples

The `run/examples/` directory contains several pre-configured example systems:

- **Spring Pendulum**: A mass on a spring that can swing and oscillate
- **Double Pendulum**: Classic chaotic double pendulum system
- **Other examples**: Check the `run/examples/` directory for more configurations

To use an example configuration, copy the desired example file to `run/config.py`:

```bash
cp examples/config_spring_pendulum.py config.py
```

## Configuration Reference

### Constants

- **`NUMBER_OF_OBJECTS`** (int): Number of independent systems to simulate simultaneously
- **`POINTS_PER_OBJECT`** (int): Number of points to plot per object (e.g., 1 for single mass, 2 for double pendulum)
- **`NUMBER_OF_COORDINATES`** (int): Number of generalized coordinates per object (degrees of freedom)
- **`DT`** (float): Timestep for numerical integration (smaller = more accurate but slower)
- **`CONNECTION_TYPE`** (ConnectionType): How to draw connections between points (e.g., `ConnectionType.BETWEEN_POINTS_AND_ORIGIN`)

### Initial Conditions

**`INITIAL_q`**: Nested list of initial generalized coordinates
- Type: `list[list[float]]`
- Shape: `[NUMBER_OF_OBJECTS][NUMBER_OF_COORDINATES]`
- Example: `[[1.0, np.pi], [1.2, np.pi/3]]` for 2 objects with 2 coordinates each

**`INITIAL_dq`**: Nested list of initial generalized velocities
- Type: `list[list[float]]`
- Shape: `[NUMBER_OF_OBJECTS][NUMBER_OF_COORDINATES]`
- Example: `[[0.0, 0.0], [0.0, 0.0]]` for 2 objects starting at rest

### Function Specifications

#### `lagrangian()`

**Returns**: `tuple[sympy.Expr, list[sympy.Function]]`
- `L`: The Lagrangian (T - V) as a SymPy expression
- `coordinates`: List of generalized coordinates as dynamicsymbols (e.g., `[r, theta]`)

**Example**:
```python
def lagrangian():
    r, theta = ds('r theta')
    # ... define T and V ...
    L = T - V
    return L, [r, theta]
```

#### `coordinates_transform(coordinates_list)`

**Arguments**: 
- `coordinates_list`: List of coordinate arrays, one per object
  - Type: `list[list[float]]`
  - Shape: `[NUMBER_OF_OBJECTS][NUMBER_OF_COORDINATES]`

**Returns**: `list[list[float]]`
- Cartesian coordinates for visualization
- Shape: `[NUMBER_OF_OBJECTS * POINTS_PER_OBJECT][2]`
- Each point is `[x, y]` 

**Example**:
```python
def coordinates_transform(coordinates_list):
    coords = []
    for i in range(len(coordinates_list)):
        r = coordinates_list[i][0]
        theta = coordinates_list[i][1]
        x = r * np.sin(theta)
        y = -r * np.cos(theta)
        coords.append([x, y])
    return coords
```

