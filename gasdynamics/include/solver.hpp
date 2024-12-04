#pragma once

#include <span>
#include <vector>

namespace gs {

template <typename T>
using span = std::span<T>;

using real    = float;
using integer = float;

class hll_solver {
public:
    hll_solver(
        span<real> density,
        span<real> velocity,
        span<real> pressure,
        span<real> argument,
        real gamma,
        real c);

    void step();
    real time();

private:
    real estimate_dt();
    void update_conservative();
    void update_primitive();

    span<real> density_;
    span<real> velocity_;
    span<real> pressure_;
    span<real> argument_;
    real gamma_;
    real c_;
    real t_;
    real dx_;
    size_t n_;
    std::vector<real> momentum_;
    std::vector<real> energy_;
};

}
