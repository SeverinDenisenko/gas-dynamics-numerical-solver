#pragma once

/*
 * Solves gas discontinuity using my other project:
 *   github.com/SeverinDenisenko/riemann-problem-solver
 *
 * `*_out` must contain `grid_divisions` pre-allocated floats.
 */

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
    float* argument_out);

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
    float* argument_out);
