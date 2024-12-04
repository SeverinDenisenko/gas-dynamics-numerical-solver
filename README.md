# HLL solver for gas dynamics

Uses HLL method for solving 1D gas dynamics tasks

## Build

```
./build.sh
./setup_venv.sh
```

## Usage as Python library

```Python
from pygasdynamics import solver

gas_state_left = solver.gas_state(1, 0, 3) # density, velocity, pressure
gas_state_right = solver.gas_state(1, 0, 1)
gas_state = solver.gas_discontinuity_state(
    gas_state_left, gas_state_right, 5 / 3 # gamma
)

solver.GRID_DIVISIONS = 320

# Numeric solution
x, density, velocity, pressure = solver.solve_discontinuity_numeric(
    test_1_gas_state, 0.0, -0.5, 0.5, 0.1, 0.9 # x_0, x_min, x_max, time, CFL
)

# Exact solution
x, density, velocity, pressure = solver.solve_discontinuity_exact(
    test_3_gas_state, 0.0, -0.5, 0.5, 0.1
)
```

## Usage as C++ library

```C++
span<real> density(/* your data */);
span<real> velocity(/* your data */);
span<real> pressure(/* your data */);
span<real> argument(/* your data */);

hll_solver solver(density, velocity, pressure, argument, gamma, c);

while (solver.time() < time) {
    solver.step();
}
```
