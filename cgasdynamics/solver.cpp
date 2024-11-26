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

    for (integer i = 0; i < grid_divisions; ++i) {
        real argument   = x1 + dx * i;
        gas_state state = solution(time, argument);
        density_out[i]  = state.density;
        velocity_out[i] = state.velocity;
        pressure_out[i] = state.pressure;
        argument_out[i] = argument;
    }
}
}
