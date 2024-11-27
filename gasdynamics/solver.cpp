#include "solver.hpp"

namespace gs {

hll_solver::hll_solver(
    span<real> density,
    span<real> velocity,
    span<real> pressure,
    span<real> argument,
    real gamma,
    real c)
    : density_(density)
    , velocity_(velocity)
    , pressure_(pressure)
    , argument_(argument)
    , gamma_(gamma)
    , c_(c)
    , t_(0)
{
}

real hll_solver::estimate_dt()
{
    real dt   = 1e+10;
    real dx   = argument_[1] - argument_[0];
    integer n = argument_.size();

    for (integer i = 0; i < n; i++) {
        real dt_i = dx
            / std::abs(velocity_[i]
                       + sqrt(gamma_ * pressure_[i] / density_[i]));
        dt = std::min(dt_i, dt);
    }

    return dt * c_;
}

void hll_solver::step()
{
    real dt = estimate_dt();

    t_ += dt;
}

real hll_solver::time()
{
    return t_;
}

}
