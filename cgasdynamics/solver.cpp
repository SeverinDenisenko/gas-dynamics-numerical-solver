#include "solver.hpp"
#include "riemann.hpp"

extern "C" {

#include "solver.h"

void gs_solve_discontinuity_exact(
    float density_left,
    float velocity_left,
    float pressure_left,
    float density_right,
    float velocity_right,
    float pressure_right,
    float gamma,
    float discontinuity_position,
    float left_border,
    float right_border,
    float time,
    int grid_divisions,
    float* density_out,
    float* velocity_out,
    float* pressure_out,
    float* argument_out)
{
    using namespace riemann;

    gas_discontinuity discontinuity {
        .left  = gas_state { .density  = density_left,
                             .velocity = velocity_left,
                             .pressure = pressure_left },
        .right = gas_state { .density  = density_right,
                             .velocity = velocity_right,
                             .pressure = pressure_right },
        .gamma = gamma
    };

    solver solver(discontinuity);
    auto solution = solver.solve();

    real x1 = left_border + discontinuity_position;
    real x2 = right_border + discontinuity_position;
    real dx = (x2 - x1) / grid_divisions;

    for (integer i = 0; i < (integer)grid_divisions; ++i) {
        real argument   = x1 + dx * i;
        gas_state state = solution(time, argument);
        density_out[i]  = state.density;
        velocity_out[i] = state.velocity;
        pressure_out[i] = state.pressure;
        argument_out[i] = argument;
    }
}

void gs_solve_discontinuity_numeric_hll(
    float density_left,
    float velocity_left,
    float pressure_left,
    float density_right,
    float velocity_right,
    float pressure_right,
    float gamma,
    float discontinuity_position,
    float left_border,
    float right_border,
    float time,
    float c,
    int grid_divisions,
    float* density_out,
    float* velocity_out,
    float* pressure_out,
    float* argument_out)
{
    using namespace gs;
    using riemann::integer;

    span<real> density(density_out, grid_divisions);
    span<real> velocity(velocity_out, grid_divisions);
    span<real> pressure(pressure_out, grid_divisions);
    span<real> argument(argument_out, grid_divisions);

    real x1 = left_border + discontinuity_position;
    real x2 = right_border + discontinuity_position;
    real dx = (x2 - x1) / grid_divisions;

    for (integer i = 0; i < (integer)grid_divisions; ++i) {
        real arg    = x1 + dx * i;
        argument[i] = arg;
        if (arg <= discontinuity_position) {
            density[i]  = density_left;
            velocity[i] = velocity_left;
            pressure[i] = pressure_left;
        } else {
            density[i]  = density_right;
            velocity[i] = velocity_right;
            pressure[i] = pressure_right;
        }
    }

    hll_solver solver(density, velocity, pressure, argument, gamma, c);

    while (solver.time() < time) {
        solver.step();
    }
}
}
