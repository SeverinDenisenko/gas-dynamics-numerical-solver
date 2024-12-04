from . import library
import ctypes

import numpy as np
from dataclasses import dataclass


GRID_DIVISIONS = 1000


@dataclass
class gas_state:
    density: float
    velocity: float
    pressure: float


@dataclass
class gas_discontinuity_state:
    left: gas_state
    right: gas_state
    gamma: float


def solve_discontinuity_exact(
    state: gas_discontinuity_state,
    discontinuity_position: float,
    left_border: float,
    right_border: float,
    time: float,
):
    density = np.zeros(GRID_DIVISIONS, dtype=np.float32)
    velocity = np.zeros(GRID_DIVISIONS, dtype=np.float32)
    pressure = np.zeros(GRID_DIVISIONS, dtype=np.float32)
    argument = np.zeros(GRID_DIVISIONS, dtype=np.float32)

    library.gs_solve_discontinuity_exact(
        state.left.density,
        state.left.velocity,
        state.left.pressure,
        state.right.density,
        state.right.velocity,
        state.right.pressure,
        state.gamma,
        discontinuity_position,
        left_border,
        right_border,
        time,
        ctypes.c_int(GRID_DIVISIONS),
        density.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        velocity.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        argument.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    )

    return argument, density, velocity, pressure


def solve_discontinuity_numeric(
    state: gas_discontinuity_state,
    discontinuity_position: float,
    left_border: float,
    right_border: float,
    time: float,
    c: float,
):
    density = np.zeros(GRID_DIVISIONS, dtype=np.float32)
    velocity = np.zeros(GRID_DIVISIONS, dtype=np.float32)
    pressure = np.zeros(GRID_DIVISIONS, dtype=np.float32)
    argument = np.zeros(GRID_DIVISIONS, dtype=np.float32)

    library.gs_solve_discontinuity_numeric_hll(
        state.left.density,
        state.left.velocity,
        state.left.pressure,
        state.right.density,
        state.right.velocity,
        state.right.pressure,
        state.gamma,
        discontinuity_position,
        left_border,
        right_border,
        time,
        c,
        ctypes.c_int(GRID_DIVISIONS),
        density.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        velocity.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        pressure.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        argument.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
    )

    return argument, density, velocity, pressure
