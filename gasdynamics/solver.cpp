#include "solver.hpp"

#include <algorithm>
#include <array>

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
    , dx_(argument_[1] - argument_[0])
    , n_(argument_.size())
    , momentum_(n_)
    , energy_(n_)
{
    update_conservative();
}

real hll_solver::estimate_dt()
{
    real vmax = 0.0;
    for (integer i = 0; i < n_; i++) {
        real v = std::abs(velocity_[i])
            + std::sqrt(gamma_ * pressure_[i] / density_[i]);
        vmax = std::max(v, vmax);
    }

    real dt = c_ * dx_ / vmax;

    return dt;
}

void hll_solver::step()
{
    real dt = estimate_dt();

    std::vector<real> flux_dencity(n_);
    std::vector<real> flux_momentum(n_);
    std::vector<real> flux_energy(n_);

    // indecies here are half-numbers
    for (integer i = 1; i < n_ - 2; ++i) {
        double c_left  = std::sqrt(gamma_ * pressure_[i] / density_[i]);
        double c_rigth = std::sqrt(gamma_ * pressure_[i + 1] / density_[i + 1]);

        real dl = std::min(velocity_[i], velocity_[i + 1])
            - std::max(c_left, c_rigth);
        real dr = std::max(velocity_[i], velocity_[i + 1])
            + std::max(c_left, c_rigth);

        real flux_left_dencity  = momentum_[i];
        real flux_left_momentum = momentum_[i] * velocity_[i] + pressure_[i];
        real flux_left_energy   = velocity_[i] * (energy_[i] + pressure_[i]);

        real flux_right_dencity = momentum_[i + 1];
        real flux_right_momentum
            = momentum_[i + 1] * velocity_[i + 1] + pressure_[i + 1];
        real flux_right_energy
            = velocity_[i + 1] * (energy_[i + 1] + pressure_[i + 1]);

        real flux_center_dencity
            = (dr * flux_left_dencity - dl * flux_right_dencity
               + dl * dr * (density_[i + 1] - density_[i]))
            / (dr - dl);
        real flux_center_momentum
            = (dr * flux_left_momentum - dl * flux_right_momentum
               + dl * dr * (momentum_[i + 1] - momentum_[i]))
            / (dr - dl);
        real flux_center_energy
            = (dr * flux_left_energy - dl * flux_right_energy
               + dl * dr * (energy_[i + 1] - energy_[i]))
            / (dr - dl);

        if (dl > 0) {
            flux_dencity[i]  = flux_left_dencity;
            flux_momentum[i] = flux_left_momentum;
            flux_energy[i]   = flux_left_energy;
        } else if (dl <= 0 && dr >= 0) {
            flux_dencity[i]  = flux_center_dencity;
            flux_momentum[i] = flux_center_momentum;
            flux_energy[i]   = flux_center_energy;
        } else {
            flux_dencity[i]  = flux_right_dencity;
            flux_momentum[i] = flux_right_momentum;
            flux_energy[i]   = flux_right_energy;
        }
    }

    for (integer i = 2; i < n_ - 2; ++i) {
        density_[i] -= dt / dx_ * (flux_dencity[i] - flux_dencity[i - 1]);
        momentum_[i] -= dt / dx_ * (flux_momentum[i] - flux_momentum[i - 1]);
        energy_[i] -= dt / dx_ * (flux_energy[i] - flux_energy[i - 1]);
    }

    density_[1]  = density_[0];
    momentum_[1] = momentum_[0];
    energy_[1]   = energy_[0];

    density_[n_ - 2]  = density_[n_ - 1];
    momentum_[n_ - 2] = momentum_[n_ - 1];
    energy_[n_ - 2]   = energy_[n_ - 1];

    t_ += dt;

    update_primitive();
}

void hll_solver::update_conservative()
{
    for (integer i = 0; i < n_; ++i) {
        momentum_[i] = density_[i] * velocity_[i];
        energy_[i]   = pressure_[i] / (gamma_ - 1)
            + 0.5 * density_[i] * std::pow(velocity_[i], 2);
    }
}

void hll_solver::update_primitive()
{
    for (integer i = 0; i < n_; ++i) {
        velocity_[i] = momentum_[i] / density_[i];
        pressure_[i]
            = (energy_[i] - 0.5 * density_[i] * std::pow(velocity_[i], 2))
            * (gamma_ - 1);
    }
}

real hll_solver::time()
{
    return t_;
}

}
