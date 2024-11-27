from test.typinganndata.ann_module3 import C_OK
from . import library
import ctypes

import numpy as np


GRID_DIVISIONS = 1000


def solve_discontinuity_exact(
    density_left: float,
    velocity_left: float,
    pressure_left: float,
    density_right: float,
    velocity_right: float,
    pressure_right: float,
    gamma: float,
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
        density_left,
        velocity_left,
        pressure_left,
        density_right,
        velocity_right,
        pressure_right,
        gamma,
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
    density_left: float,
    velocity_left: float,
    pressure_left: float,
    density_right: float,
    velocity_right: float,
    pressure_right: float,
    gamma: float,
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
        density_left,
        velocity_left,
        pressure_left,
        density_right,
        velocity_right,
        pressure_right,
        gamma,
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
