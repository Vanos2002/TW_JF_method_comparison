# This code is called SATURDAY.py in VS Code
# THIS CODE WAS CREATED ON THURSDAY JUNE 4, 2026
# THIS FILE COMPARES THE FEREISL-TUCKERWILL 4.5PN EXPRESSIONS TO THE RK4 ADAPTIVE CODE
#
# From the article "Residual eccentricity of inspiralling orbits at the
# gravitational-wave detection threshold: Accurate estimates using
# post-Newtonian theory" by Alexandria Tucker and Clifford M. Will
# (arXiv:2108.12210v2 [gr-qc] 15 Nov 2021)
#
# We compare the transformed 4.5PN contributions of dp/dtheta and de/dtheta
# from the Fereisl-Tucker-Will (FTW) paper against the numerical
# orbit-averaged QLT results.

from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field
from typing import Callable, List, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants (geometrized units)
# ---------------------------------------------------------------------------
G = 1.0           # Gravitational constant in geometrized units
c = 1.0           # Speed of light in geometrized units

M_sun       = 1.98847e30    # Solar mass in kg
G_SI        = 6.6743e-11    # Gravitational constant  m^3/kg/s^2
c_SI        = 299792458.0   # Speed of light  m/s
PI          = math.pi
PI2         = PI * PI

# Masses
m1_solar      = 1.0
m2_solar      = 1.0
M_total_solar = m1_solar + m2_solar

m1  = m1_solar / M_total_solar
m2  = m2_solar / M_total_solar
M   = m1 + m2
mu  = m1 * m2 / M
eta = mu / M
Mc  = (m1 * m2) ** 0.6 * M ** (-0.2)
M_kg = M_total_solar * M_sun

time_unit_seconds = G_SI * M_kg / (c_SI ** 3)
sep_unit_meters   = G_SI * M_kg / (c_SI ** 2)

# ---------------------------------------------------------------------------
# PN coefficient tables
# ---------------------------------------------------------------------------
@dataclass
class PNCoeffs:
    """Holds 3-D coefficient tables a, b, c, d of shape (4, 4, 4)."""
    a: np.ndarray = field(default_factory=lambda: np.zeros((4, 4, 4)))
    b: np.ndarray = field(default_factory=lambda: np.zeros((4, 4, 4)))
    c: np.ndarray = field(default_factory=lambda: np.zeros((4, 4, 4)))
    d: np.ndarray = field(default_factory=lambda: np.zeros((4, 4, 4)))


def build_coefficients(include_4p5PN: bool = True) -> PNCoeffs:
    """Populate all PN coefficient tables and return them."""
    e2 = eta * eta
    e3 = e2  * eta

    K = PNCoeffs()
    a, b, cc, d = K.a, K.b, K.c, K.d

    # ── 1 PN ──────────────────────────────────────────────────────────────
    a[1,0,0] =  2.0 * (2.0 + eta)
    a[0,1,0] =  1.5 * eta
    a[0,0,1] = -(1.0 + 3.0 * eta)

    b[1,0,0] =  0.0
    b[0,1,0] =  2.0 * (2.0 - eta)
    b[0,0,1] =  0.0

    # ── 2 PN ──────────────────────────────────────────────────────────────
    a[2,0,0] = -0.75 * (12.0 + 29.0 * eta)
    a[0,2,0] = -1.875 * eta * (1.0 - 3.0 * eta)
    a[0,0,2] = -eta * (3.0 - 4.0 * eta)
    a[1,1,0] =  2.0 + 25.0 * eta + 2.0 * e2
    a[1,0,1] =  0.5 * eta * (13.0 - 4.0 * eta)
    a[0,1,1] =  1.5 * eta * (3.0 - 4.0 * eta)

    b[2,0,0] =  0.0
    b[0,2,0] = -1.5 * eta * (3.0 + 2.0 * eta)
    b[0,0,2] =  0.0
    b[1,1,0] = -0.5 * (4.0 + 41.0 * eta + 8.0 * e2)
    b[1,0,1] =  0.0
    b[0,1,1] =  0.5 * eta * (15.0 + 4.0 * eta)

    # ── 3 PN ──────────────────────────────────────────────────────────────
    a[3,0,0] =  16.0 + eta * (5596.0 - 123.0 * PI2 + 1704.0 * eta) / 48.0
    a[0,3,0] =  (35.0 / 16.0) * eta * (1.0 - 5.0 * eta + 5.0 * e2)
    a[0,0,3] = -(1.0 / 4.0) * eta * (11.0 - 49.0 * eta + 52.0 * e2)
    a[2,1,0] = (-1.0 - (22717.0 / 168.0 + 615.0 / 64.0 * PI2) * eta
                - (11.0 / 8.0) * e2 + 7.0 * e3)
    a[2,0,1] =  eta * (20827.0 / 840.0 + 123.0 / 64.0 * PI2 - e2)
    a[1,2,0] = -0.5 * eta * (158.0 - 69.0 * eta - 60.0 * e2)
    a[1,1,1] =  eta * (121.0 - 16.0 * eta - 20.0 * e2)
    a[1,0,2] = -(1.0 / 4.0) * eta * (75.0 + 32.0 * eta - 40.0 * e2)
    a[0,2,1] = -(15.0 / 8.0) * eta * (4.0 - 18.0 * eta + 17.0 * e2)
    a[0,1,2] =  (3.0 / 8.0) * eta * (20.0 - 79.0 * eta + 60.0 * e2)

    b[3,0,0] =  0.0
    b[0,3,0] =  (15.0 / 8.0) * eta * (3.0 - 8.0 * eta - 2.0 * e2)
    b[0,0,3] =  0.0
    b[2,1,0] =  (4.0 + (5849.0 / 840.0 + 123.0 / 32.0 * PI2) * eta
                 - 25.0 * e2 - 8.0 * e3)
    b[2,0,1] =  0.0
    b[1,2,0] = -(1.0 / 6.0) * eta * (329.0 + 177.0 * eta + 108.0 * e2)
    b[1,1,1] =  eta * (15.0 + 27.0 * eta + 10.0 * e2)
    b[1,0,2] =  0.0
    b[0,2,1] = -(3.0 / 4.0) * eta * (16.0 - 37.0 * eta - 16.0 * e2)
    b[0,1,2] =  (1.0 / 8.0) * eta * (65.0 - 152.0 * eta - 48.0 * e2)

    # ── 2.5 PN  (radiation reaction) ──────────────────────────────────────
    cc[1,0,0] =  17.0 / 3.0
    cc[0,1,0] =  0.0
    cc[0,0,1] =  3.0

    d[1,0,0] = -3.0
    d[0,1,0] =  0.0
    d[0,0,1] = -1.0

    # ── 3.5 PN  (radiation reaction) ──────────────────────────────────────
    cc[2,0,0] = -(23.0 / 14.0) * (43.0 + 14.0 * eta)
    cc[0,2,0] = -70.0
    cc[0,0,2] = -(3.0 / 28.0) * (61.0 + 70.0 * eta)
    cc[1,1,0] = -(1.0 / 4.0) * (147.0 + 188.0 * eta)
    cc[1,0,1] = -(1.0 / 42.0) * (519.0 - 1267.0 * eta)
    cc[0,1,1] =  (15.0 / 4.0) * (19.0 + 2.0 * eta)

    d[2,0,0] =  (1.0 / 42.0) * (1325.0 + 546.0 * eta)
    d[0,2,0] =  75.0
    d[0,0,2] =  (1.0 / 28.0) * (313.0 + 42.0 * eta)
    d[1,1,0] =  (1.0 / 12.0) * (205.0 + 424.0 * eta)
    d[1,0,1] = -(1.0 / 42.0) * (205.0 + 777.0 * eta)
    d[0,1,1] = -(3.0 / 4.0) * (113.0 + 2.0 * eta)

    # ── 4.5 PN  (radiation reaction — Fereisl term) ────────────────────────
    if include_4p5PN:
        cc[3,0,0] =  (1.0 / 756.0) * (289079.0 + 284127.0 * eta + 22632.0 * e2)
        cc[0,3,0] =  0.0
        cc[0,0,3] =  (1.0 / 168.0) * (779.0 + 604.0 * eta - 7090.0 * e2)
        cc[2,1,0] =  (1.0 / 756.0) * (250221.0 - 6032.0 * eta + 74134.0 * e2)
        cc[2,0,1] = -(1.0 / 252.0) * (20916.0 - 24324.0 * eta + 23483.0 * e2)
        cc[1,2,0] =  (1.0 / 252.0) * (108322.0 - 43996.0 * eta + 12839.0 * e2)
        cc[1,1,1] = -(1.0 / 504.0) * (218401.0 - 160227.0 * eta + 95987.0 * e2)
        cc[1,0,2] =  (1.0 / 504.0) * (40758.0 - 88311.0 * eta + 43474.0 * e2)
        cc[0,2,1] =  (5.0 / 18.0) * (87.0 - 215.0 * eta - 97.0 * e2)
        cc[0,1,2] = -(1.0 / 84.0) * (1205.0 - 260.0 * eta - 8785.0 * e2)

        d[3,0,0] = -(1.0 / 2268.0) * (395929.0 + 398700.0 * eta + 87048.0 * e2)
        d[0,3,0] =  (5.0 / 18.0) * (291.0 - 919.0 * eta + 97.0 * e2)
        d[0,0,3] = -(1.0 / 56.0) * (834.0 - 1956.0 * eta - 1743.0 * e2)
        d[2,1,0] = -(1.0 / 252.0) * (37992.0 + 62832.0 * eta + 9649.0 * e2)
        d[2,0,1] =  (1.0 / 252.0) * (26703.0 + 21304.0 * eta + 28486.0 * e2)
        d[1,2,0] = -(1.0 / 252.0) * (99499.0 + 24002.0 * eta + 33443.0 * e2)
        d[1,1,1] =  (1.0 / 504.0) * (200244.0 + 65460.0 * eta + 83501.0 * e2)
        d[1,0,2] = -(1.0 / 504.0) * (16731.0 + 24785.0 * eta + 41471.0 * e2)
        d[0,2,1] = -(5.0 / 168.0) * (6889.0 - 21631.0 * eta + 2380.0 * e2)
        d[0,1,2] =  (1.0 / 168.0) * (21280.0 - 60733.0 * eta - 11999.0 * e2)

    return K


# ---------------------------------------------------------------------------
# Orbital kinematics
# ---------------------------------------------------------------------------

def norm_r(p, e, alpha, beta, phi):
    return p / (1.0 + e * (alpha * math.cos(phi) + beta * math.sin(phi)))


def r_dot(p, e, alpha, beta, phi):
    return math.sqrt(G * M * p) * e / p * (alpha * math.sin(phi) - beta * math.cos(phi))


def norm_v2(p, e, alpha, beta, phi):
    r  = norm_r(p, e, alpha, beta, phi)
    rd = r_dot(p, e, alpha, beta, phi)
    return rd * rd + G * M * p / (r * r)


# ---------------------------------------------------------------------------
# Monomial sum over one coefficient table
# ---------------------------------------------------------------------------

def sum_table(tab: np.ndarray, N: int,
              rd2: float, v2: float, gmr: float, c_val: float) -> float:
    s = 0.0
    c2N = c_val ** (2 * N)
    for l in range(N + 1):
        for m in range(N - l + 1):
            n = N - l - m
            if n < 0:
                continue
            s += tab[l, m, n] * (rd2 ** m) * (v2 ** n) * (gmr ** l) / c2N
    return s


# ---------------------------------------------------------------------------
# QLT equations of motion
# ---------------------------------------------------------------------------

@dataclass
class QLTrhs:
    dp: float = 0.0
    dalpha: float = 0.0
    dbeta: float = 0.0


def compute_qlt(K: PNCoeffs, p: float, alpha: float, beta: float,
                phi: float, c_val: float, pn_order: int) -> QLTrhs:
    e   = math.sqrt(alpha * alpha + beta * beta)
    r   = norm_r(p, e, alpha, beta, phi)
    rd  = r_dot(p, e, alpha, beta, phi)
    v2  = norm_v2(p, e, alpha, beta, phi)
    rd2 = rd * rd
    gmr = G * M / r
    eps = 1.0 / c_val

    Atot = 0.0
    Btot = 0.0
    for N in range(1, pn_order + 1):
        Atot += sum_table(K.a, N, rd2, v2, gmr, c_val)
        Btot += sum_table(K.b, N, rd2, v2, gmr, c_val)

    Crr = 0.0
    Drr = 0.0
    for N in range(1, pn_order + 1):
        Crr += sum_table(K.c, N, rd2, v2, gmr, c_val)
        Drr += sum_table(K.d, N, rd2, v2, gmr, c_val)

    eps3 = eps ** 3
    rr_R_pre = (8.0 / 5.0) * eta * eps3 * (G * M) ** 2 / r ** 3 * rd
    rr_S_pre = (8.0 / 5.0) * eta * eps3 * (G * M / r ** 2) ** 2 * math.sqrt(G * M * p)

    ScR = G * M / r ** 2 * (Atot + Btot) + rr_R_pre * (Crr + Drr)
    ScS = G * M / r ** 3 * math.sqrt(G * M * p) * rd * Btot + rr_S_pre * Drr

    dp_dphi = 2.0 * r ** 3 / (G * M) * ScS
    dalpha  = (r ** 2 / (G * M) * (ScR * math.sin(phi)
               + ScS * (alpha + math.cos(phi)) * (1.0 + r / p)
               - ScS * alpha))
    dbeta   = (r ** 2 / (G * M) * (-ScR * math.cos(phi)
               + ScS * (beta + math.sin(phi)) * (1.0 + r / p)
               - ScS * beta))

    return QLTrhs(dp_dphi, dalpha, dbeta)


# ---------------------------------------------------------------------------
# Physical structs
# ---------------------------------------------------------------------------

@dataclass
class PhysicalParams:
    G:   float = 1.0
    M:   float = 1.0
    eta: float = 0.25
    phi: float = 0.0
    eps: float = 1.0


@dataclass
class BinaryState:
    p:     float = 0.0
    alpha: float = 0.0
    beta:  float = 0.0


SecularRHS = List[float]   # [dp, dalpha, dbeta]

PHI_AVERAGE_SAMPLES = 64


# ---------------------------------------------------------------------------
# Dissipative-only QLT for a single order
# ---------------------------------------------------------------------------

def compute_qlt_dissipative_order(K: PNCoeffs, state: BinaryState,
                                  params: PhysicalParams,
                                  dissipative_order: int) -> QLTrhs:
    p, alpha, beta = state.p, state.alpha, state.beta
    phi   = params.phi
    c_val = 1.0 / params.eps
    e     = math.sqrt(alpha * alpha + beta * beta)
    r     = norm_r(p, e, alpha, beta, phi)
    rd    = r_dot(p, e, alpha, beta, phi)
    v2    = norm_v2(p, e, alpha, beta, phi)
    rd2   = rd * rd
    GM    = params.G * params.M
    gmr   = GM / r

    Crr = sum_table(K.c, dissipative_order, rd2, v2, gmr, c_val)
    Drr = sum_table(K.d, dissipative_order, rd2, v2, gmr, c_val)

    eps  = params.eps
    eps3 = eps ** 3
    rr_R_pre = (8.0 / 5.0) * params.eta * eps3 * GM * GM / r ** 3 * rd
    rr_S_pre = (8.0 / 5.0) * params.eta * eps3 * (GM / r ** 2) ** 2 * math.sqrt(GM * p)

    ScR = rr_R_pre * (Crr + Drr)
    ScS = rr_S_pre * Drr

    dp_dphi = 2.0 * r ** 3 / GM * ScS
    dalpha  = (r ** 2 / GM * (ScR * math.sin(phi)
               + ScS * (alpha + math.cos(phi)) * (1.0 + r / p)
               - ScS * alpha))
    dbeta   = (r ** 2 / GM * (-ScR * math.cos(phi)
               + ScS * (beta + math.sin(phi)) * (1.0 + r / p)
               - ScS * beta))

    return QLTrhs(dp_dphi, dalpha, dbeta)


def average_qlt_rhs(f: Callable[[float], QLTrhs],
                    samples: int = PHI_AVERAGE_SAMPLES) -> QLTrhs:
    h = 2.0 * PI / samples
    dp_sum = dalpha_sum = dbeta_sum = 0.0
    for i in range(samples):
        phi = (i + 0.5) * h
        val = f(phi)
        dp_sum    += val.dp
        dalpha_sum += val.dalpha
        dbeta_sum  += val.dbeta
    scale = h / (2.0 * PI)
    return QLTrhs(dp_sum * scale, dalpha_sum * scale, dbeta_sum * scale)


# ---------------------------------------------------------------------------
# Secular PN terms (Fereisl)
# ---------------------------------------------------------------------------

def secular_2_5PN_explicit(K: PNCoeffs, state: BinaryState,
                            params: PhysicalParams) -> SecularRHS:
    def f(phi):
        pp = PhysicalParams(params.G, params.M, params.eta, phi, params.eps)
        return compute_qlt_dissipative_order(K, state, pp, 1)
    avg = average_qlt_rhs(f)
    return [avg.dp, avg.dalpha, avg.dbeta]


def secular_3_5PN_explicit(K: PNCoeffs, state: BinaryState,
                            params: PhysicalParams) -> SecularRHS:
    h = 2.0 * PI / PHI_AVERAGE_SAMPLES
    dp = dalpha = dbeta = 0.0
    for i in range(PHI_AVERAGE_SAMPLES):
        phi = (i + 0.5) * h
        pp  = PhysicalParams(params.G, params.M, params.eta, phi, params.eps)
        F35 = compute_qlt_dissipative_order(K, state, pp, 2)
        dp    += F35.dp
        dalpha += F35.dalpha
        dbeta  += F35.dbeta
    scale = h / (2.0 * PI)
    return [dp * scale, dalpha * scale, dbeta * scale]


def secular_4_5PN_explicit(K: PNCoeffs, state: BinaryState,
                            params: PhysicalParams) -> SecularRHS:
    h = 2.0 * PI / PHI_AVERAGE_SAMPLES
    dp = dalpha = dbeta = 0.0
    for i in range(PHI_AVERAGE_SAMPLES):
        phi = (i + 0.5) * h
        pp  = PhysicalParams(params.G, params.M, params.eta, phi, params.eps)
        F45 = compute_qlt_dissipative_order(K, state, pp, 3)
        dp    += F45.dp
        dalpha += F45.dalpha
        dbeta  += F45.dbeta
    scale = h / (2.0 * PI)
    return [dp * scale, dalpha * scale, dbeta * scale]


def secular_1PN(state: BinaryState, params: PhysicalParams) -> SecularRHS:
    p     = state.p
    GM    = params.G * params.M
    eps2  = params.eps ** 2
    return [
        0.0,
        eps2 * (-(3.0 * GM * state.beta) / p),
        eps2 * ( (3.0 * GM * state.alpha) / p),
    ]


def oscillatory_1PN(state: BinaryState, params: PhysicalParams) -> List[float]:
    p, alpha, beta = state.p, state.alpha, state.beta
    phi = params.phi
    GM  = params.G * params.M

    cos_phi  = math.cos(phi)
    cos_2phi = math.cos(2.0 * phi)
    cos_3phi = math.cos(3.0 * phi)
    sin_phi  = math.sin(phi)
    sin_2phi = math.sin(2.0 * phi)
    sin_3phi = math.sin(3.0 * phi)

    Y0 = 4.0 * GM * (-2.0 + eta) * (alpha * cos_phi + beta * sin_phi)

    Y1 = (1.0 / (8.0 * p)) * GM * (
        (8.0 * (-3.0 + eta) + beta * beta * (8.0 + 21.0 * eta)
         + alpha * alpha * (-56.0 + 47.0 * eta)) * cos_phi
        + 4.0 * alpha * (-5.0 + 4.0 * eta) * cos_2phi
        + (alpha * alpha - beta * beta) * eta * cos_3phi
        + 2.0 * alpha * beta * (-32.0 + 13.0 * eta) * sin_phi
        + 4.0 * beta * (-5.0 + 4.0 * eta) * sin_2phi
        + 2.0 * alpha * beta * eta * sin_3phi
    )

    Y2 = (1.0 / (8.0 * p)) * GM * (
        2.0 * alpha * beta * (-32.0 + 13.0 * eta) * cos_phi
        + 4.0 * beta * (5.0 - 4.0 * eta) * cos_2phi
        - 2.0 * alpha * beta * eta * cos_3phi
        + (8.0 * (-3.0 + eta) + alpha * alpha * (8.0 + 21.0 * eta)
           + beta * beta * (-56.0 + 47.0 * eta)) * sin_phi
        + 4.0 * alpha * (-5.0 + 4.0 * eta) * sin_2phi
        + (alpha * alpha - beta * beta) * eta * sin_3phi
    )

    return [Y0, Y1, Y2]


def secular_2PN(state: BinaryState, params: PhysicalParams) -> SecularRHS:
    p, alpha, beta = state.p, state.alpha, state.beta
    GM   = params.G * params.M
    GM2  = GM * GM
    eps4 = params.eps ** 4

    bracket = (-10.0 + beta * beta - 4.0 * eta + 10.0 * beta * beta * eta
               + alpha * alpha * (1.0 + 10.0 * eta))

    return [
        0.0,
        eps4 * (-(3.0 * GM2 * beta * bracket) / (4.0 * p * p)),
        eps4 * ( (3.0 * GM2 * alpha * bracket) / (4.0 * p * p)),
    ]


def oscillatory_2PN(state: BinaryState, params: PhysicalParams) -> List[float]:
    p, alpha, beta = state.p, state.alpha, state.beta
    phi = params.phi
    GM  = params.G * params.M
    GM2 = GM * GM

    cos_phi  = math.cos(phi)
    cos_2phi = math.cos(2.0 * phi)
    cos_3phi = math.cos(3.0 * phi)
    cos_4phi = math.cos(4.0 * phi)
    cos_5phi = math.cos(5.0 * phi)
    sin_phi  = math.sin(phi)
    sin_2phi = math.sin(2.0 * phi)
    sin_3phi = math.sin(3.0 * phi)
    sin_4phi = math.sin(4.0 * phi)
    sin_5phi = math.sin(5.0 * phi)

    a2 = alpha * alpha
    b2 = beta  * beta
    a4 = a2 * a2
    b4 = b2 * b2

    Y0 = (1.0 / (4.0 * p)) * GM2 * (
        alpha * (-160.0 + eta * (256.0 - 33.0 * a2 - 33.0 * b2
                 + 2.0 * (-8.0 + a2 + b2) * eta)) * cos_phi
        + (a2 - b2) * (68.0 + 3.0 * eta * (-15.0 + 4.0 * eta)) * cos_2phi
        - alpha * (a2 - 3.0 * b2) * eta * (3.0 + 2.0 * eta) * cos_3phi
        + beta * (-160.0 + eta * (256.0 - 33.0 * a2 - 33.0 * b2
                  + 2.0 * (-8.0 + a2 + b2) * eta)) * sin_phi
        + 2.0 * alpha * beta * (68.0 + 3.0 * eta * (-15.0 + 4.0 * eta)) * sin_2phi
        + beta * (-3.0 * a2 + b2) * eta * (3.0 + 2.0 * eta) * sin_3phi
    )

    Y1 = (1.0 / (128.0 * p * p)) * GM2 * (
        -2.0 * (5.0 * b4 * eta * (-27.0 + 41.0 * eta)
                + a4 * eta * (477.0 + 161.0 * eta)
                + b2 * (-944.0 + 92.0 * eta + 64.0 * eta * eta)
                - 16.0 * (60.0 + eta * (17.0 + 8.0 * eta))
                + a2 * (3056.0 + 2.0 * eta * (-1958.0 + 224.0 * eta
                        + 3.0 * b2 * (57.0 + 61.0 * eta)))) * cos_phi
        + 16.0 * alpha * (28.0 * (-1.0 + 5.0 * eta)
                          + b2 * (-162.0 + eta * (57.0 + 28.0 * eta))
                          + a2 * (106.0 + eta * (-163.0 + 60.0 * eta))) * cos_2phi
        + (-a4 * eta * (73.0 + 53.0 * eta)
           + b2 * (-480.0 - 8.0 * eta * (9.0 + 16.0 * eta)
                   + b2 * eta * (-33.0 + 19.0 * eta))
           + 2.0 * a2 * (240.0 + eta * (36.0 + 64.0 * eta
                         + 3.0 * b2 * (53.0 + 17.0 * eta)))) * cos_3phi
        - 8.0 * alpha * (a2 - 3.0 * b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * cos_4phi
        - 3.0 * (a4 - 6.0 * a2 * b2 + b4) * eta * (-1.0 + 3.0 * eta) * cos_5phi
        + 8.0 * alpha * beta * (-1000.0 + eta * (1002.0 - 96.0 * eta
                                + a2 * (-153.0 + 11.0 * eta)
                                + b2 * (-153.0 + 11.0 * eta))) * sin_phi
        + 16.0 * beta * (-28.0 * (1.0 + b2 - 5.0 * eta)
                         + b2 * eta * (-53.0 + 44.0 * eta)
                         + a2 * (240.0 + eta * (-273.0 + 76.0 * eta))) * sin_2phi
        - 4.0 * alpha * beta * (-240.0 + eta * (b2 * (-43.0 + eta)
                                + 7.0 * a2 * (9.0 + 5.0 * eta)
                                - 4.0 * (9.0 + 16.0 * eta))) * sin_3phi
        + 8.0 * beta * (-3.0 * a2 + b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * sin_4phi
        - 12.0 * alpha * (alpha - beta) * beta * (alpha + beta) * eta * (-1.0 + 3.0 * eta) * sin_5phi
    )

    Y2 = (1.0 / (128.0 * p * p)) * GM2 * (
        8.0 * alpha * beta * (-1000.0 + eta * (1002.0 - 96.0 * eta
                              + a2 * (-153.0 + 11.0 * eta)
                              + b2 * (-153.0 + 11.0 * eta))) * cos_phi
        - 16.0 * beta * (28.0 * (-1.0 + 5.0 * eta)
                         + a2 * (-162.0 + eta * (57.0 + 28.0 * eta))
                         + b2 * (106.0 + eta * (-163.0 + 60.0 * eta))) * cos_2phi
        + 4.0 * alpha * beta * (-240.0 + eta * (a2 * (-43.0 + eta)
                                + 7.0 * b2 * (9.0 + 5.0 * eta)
                                - 4.0 * (9.0 + 16.0 * eta))) * cos_3phi
        - 8.0 * beta * (-3.0 * a2 + b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * cos_4phi
        + 12.0 * alpha * (alpha - beta) * beta * (alpha + beta) * eta * (-1.0 + 3.0 * eta) * cos_5phi
        - 2.0 * (5.0 * a4 * eta * (-27.0 + 41.0 * eta)
                 + b4 * eta * (477.0 + 161.0 * eta)
                 - 16.0 * (60.0 + eta * (17.0 + 8.0 * eta))
                 + 4.0 * b2 * (764.0 + eta * (-979.0 + 112.0 * eta))
                 + a2 * (-944.0 + 2.0 * eta * (46.0 + 32.0 * eta
                          + 3.0 * b2 * (57.0 + 61.0 * eta)))) * sin_phi
        + 16.0 * alpha * (28.0 * (-1.0 + 5.0 * eta)
                          + a2 * (-28.0 + eta * (-53.0 + 44.0 * eta))
                          + b2 * (240.0 + eta * (-273.0 + 76.0 * eta))) * sin_2phi
        + (a4 * (33.0 - 19.0 * eta) * eta
           + b2 * (-480.0 - 8.0 * eta * (9.0 + 16.0 * eta)
                   + b2 * eta * (73.0 + 53.0 * eta))
           + 2.0 * a2 * (240.0 + eta * (36.0 + 64.0 * eta
                         - 3.0 * b2 * (53.0 + 17.0 * eta)))) * sin_3phi
        - 8.0 * alpha * (a2 - 3.0 * b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * sin_4phi
        - 3.0 * (a4 - 6.0 * a2 * b2 + b4) * eta * (-1.0 + 3.0 * eta) * sin_5phi
    )

    return [Y0, Y1, Y2]


def secular_2_5PN(state: BinaryState, params: PhysicalParams) -> SecularRHS:
    p, alpha, beta = state.p, state.alpha, state.beta
    GM = params.G * params.M

    sqrt_GMp  = math.sqrt(GM * p)
    GMp_5_2   = sqrt_GMp ** 5
    a2 = alpha * alpha
    b2 = beta  * beta

    dp    = -(8.0 * GMp_5_2 * (8.0 + 7.0 * a2 + 7.0 * b2) * params.eta) / (5.0 * p ** 4)
    dalpha = -(GMp_5_2 * alpha * (304.0 + 121.0 * a2 + 121.0 * b2) * params.eta) / (15.0 * p ** 5)
    dbeta  = -(GMp_5_2 * beta  * (304.0 + 121.0 * a2 + 121.0 * b2) * params.eta) / (15.0 * p ** 5)
    return [dp, dalpha, dbeta]


def secular_3_5PN(state: BinaryState, params: PhysicalParams) -> SecularRHS:
    p, alpha, beta = state.p, state.alpha, state.beta
    GM  = params.G * params.M
    GM3 = GM ** 3

    sqrt_GMp = math.sqrt(GM * p)
    GMp_3_2  = sqrt_GMp ** 3

    a2 = alpha * alpha
    b2 = beta  * beta
    a4 = a2 * a2
    b4 = b2 * b2

    term_p = (-8.0 * (2759.0 + 252.0 * params.eta)
              + 8.0 * b2 * (758.0 + 889.0 * params.eta)
              + a4 * (1483.0 + 4424.0 * params.eta)
              + b4 * (1483.0 + 4424.0 * params.eta)
              + a2 * (8.0 * (758.0 + 889.0 * params.eta)
                      + b2 * (2966.0 + 8848.0 * params.eta)))

    term_alpha = (-8.0 * (18049.0 + 4452.0 * params.eta)
                  + 4.0 * b2 * (8692.0 + 12803.0 * params.eta)
                  + a4 * (2251.0 + 15064.0 * params.eta)
                  + b4 * (2251.0 + 15064.0 * params.eta)
                  + a2 * (34768.0 + 51212.0 * params.eta
                          + b2 * (4502.0 + 30128.0 * params.eta)))

    term_beta = term_alpha  # same polynomial
    eps2 = params.eps ** 2

    dp     = -eps2 * (GM3 * GMp_3_2 * params.eta * term_p) / (210.0 * p ** 3)
    dalpha = -eps2 * (GM3 * GMp_3_2 * alpha * params.eta * term_alpha) / (840.0 * p ** 4)
    dbeta  = -eps2 * (GM3 * GMp_3_2 * beta  * params.eta * term_beta) / (840.0 * p ** 4)
    return [dp, dalpha, dbeta]


def secular_4_5PN(state: BinaryState, params: PhysicalParams) -> SecularRHS:
    """Fereisl 4.5PN secular terms."""
    p, alpha, beta = state.p, state.alpha, state.beta
    GM   = params.G * params.M
    GM4  = GM ** 4

    sqrt_GMp = math.sqrt(GM * p)
    GMp_5_2  = sqrt_GMp ** 5

    a2 = alpha * alpha
    b2 = beta  * beta
    a4 = a2 * a2
    b4 = b2 * b2
    a6 = a4 * a2
    b6 = b4 * b2
    eta2 = params.eta ** 2

    term_p = (8.0 * (-1034075.0 - 261369.0 * params.eta + 36288.0 * eta2)
              + 9.0 * a6 * (527.0 - 6300.0 * params.eta + 53088.0 * eta2)
              + 9.0 * b6 * (527.0 - 6300.0 * params.eta + 53088.0 * eta2)
              + 12.0 * b2 * (-64831.0 + 1186209.0 * params.eta + 64575.0 * eta2)
              + b4 * (947991.0 - 5386365.0 * params.eta + 3988656.0 * eta2)
              + 3.0 * a4 * (315997.0 - 1795455.0 * params.eta + 1329552.0 * eta2
                            + 9.0 * b2 * (527.0 - 6300.0 * params.eta + 53088.0 * eta2))
              + 3.0 * a2 * (9.0 * b4 * (527.0 - 6300.0 * params.eta + 53088.0 * eta2)
                            + 4.0 * (-64831.0 + 1186209.0 * params.eta + 64575.0 * eta2)
                            + b2 * (631994.0 - 3590910.0 * params.eta + 2659104.0 * eta2)))

    term_alpha = (16.0 * (-2739835.0 - 1394559.0 * params.eta + 145152.0 * eta2)
                  + 3.0 * a6 * (-25845.0 - 78380.0 * params.eta + 361536.0 * eta2)
                  + 3.0 * b6 * (-25845.0 - 78380.0 * params.eta + 361536.0 * eta2)
                  + 12.0 * b2 * (-354911.0 + 4848903.0 * params.eta + 511413.0 * eta2)
                  + 2.0 * b4 * (605645.0 - 8079297.0 * params.eta + 5758704.0 * eta2)
                  + a4 * (9.0 * b2 * (-25845.0 - 78380.0 * params.eta + 361536.0 * eta2)
                          + 2.0 * (605645.0 - 8079297.0 * params.eta + 5758704.0 * eta2))
                  + a2 * (9.0 * b4 * (-25845.0 - 78380.0 * params.eta + 361536.0 * eta2)
                          + 12.0 * (-354911.0 + 4848903.0 * params.eta + 511413.0 * eta2)
                          + 4.0 * b2 * (605645.0 - 8079297.0 * params.eta + 5758704.0 * eta2)))
    eps4 = params.eps ** 4

    dp     = eps4 * (GM4 * GMp_5_2 * params.eta * term_p) / (11340.0 * p ** 4)
    dalpha = eps4 * (GM4 * GMp_5_2 * alpha * params.eta * term_alpha) / (30240.0 * p ** 5)
    dbeta  = eps4 * (GM4 * GMp_5_2 * beta  * params.eta * term_alpha) / (30240.0 * p ** 5)
    return [dp, dalpha, dbeta]


# ---------------------------------------------------------------------------
# Composite secular RHS  (Fereisl)
# ---------------------------------------------------------------------------

_K_global: PNCoeffs | None = None

def _get_K() -> PNCoeffs:
    global _K_global
    if _K_global is None:
        _K_global = build_coefficients(True)
    return _K_global


def compute_secular_rhs(state: BinaryState, params: PhysicalParams,
                        max_pn_order: int) -> SecularRHS:
    K = _get_K()
    total = [0.0, 0.0, 0.0]

    def _add(rhs):
        total[0] += rhs[0]; total[1] += rhs[1]; total[2] += rhs[2]

    if max_pn_order >= 1:
        _add(secular_1PN(state, params))
    if max_pn_order >= 2:
        _add(secular_2PN(state, params))
    if max_pn_order >= 3:   # 2.5 PN
        _add(secular_2_5PN_explicit(K, state, params))
    if max_pn_order >= 4:   # 3.5 PN
        _add(secular_3_5PN_explicit(K, state, params))
    if max_pn_order >= 5:   # 4.5 PN
        _add(secular_4_5PN_explicit(K, state, params))

    return total


# ---------------------------------------------------------------------------
# Tucker-Will 4.5PN secular terms (differs only at 4.5PN)
# ---------------------------------------------------------------------------

def secular_4_5PN_TW(state: BinaryState, params: PhysicalParams) -> SecularRHS:
    p, alpha, beta = state.p, state.alpha, state.beta
    GM   = params.G * params.M
    GM4  = GM ** 4

    sqrt_GMp = math.sqrt(GM * p)
    GMp_5_2  = sqrt_GMp ** 5

    a2 = alpha * alpha
    b2 = beta  * beta
    a4 = a2 * a2
    b4 = b2 * b2
    a6 = a4 * a2
    b6 = b4 * b2
    eta2 = params.eta ** 2

    term_p = (-8272600.0
              + 72.0 * params.eta * (-29041.0 + 4032.0 * params.eta)
              + 9.0 * a6 * (527.0 + 84.0 * params.eta * (-75.0 + 632.0 * params.eta))
              + 9.0 * b6 * (527.0 + 84.0 * params.eta * (-75.0 + 632.0 * params.eta))
              + 12.0 * b2 * (-64831.0 + 9.0 * params.eta * (131801.0 + 7175.0 * params.eta))
              + b4 * (947991.0 + 27.0 * params.eta * (-199495.0 + 147728.0 * params.eta))
              + a4 * (9.0 * b2 * (527.0 + 84.0 * params.eta * (-75.0 + 632.0 * params.eta))
                      + (947991.0 + 27.0 * params.eta * (-199495.0 + 147728.0 * params.eta)))
              + a2 * (12.0 * (-64831.0 + 9.0 * params.eta * (131801.0 + 7175.0 * params.eta))
                      + 9.0 * b4 * (527.0 + 84.0 * params.eta * (-75.0 + 632.0 * params.eta))
                      + b2 * (2.0 * (947991.0 + 27.0 * params.eta * (-199495.0 + 147728.0 * params.eta)))))

    term_alpha = (43837360.0
                  + 144.0 * params.eta * (154951.0 - 16128.0 * params.eta)
                  + 9.0 * a6 * (8615.0 + 4.0 * params.eta * (6565.0 - 30128.0 * params.eta))
                  + 9.0 * b6 * (8615.0 + 4.0 * params.eta * (6565.0 - 30128.0 * params.eta))
                  - 12.0 * b2 * (-354911.0 + 4848303.0 * params.eta + 511413.0 * eta2)
                  - 2.0 * b4 * (605645.0 + 9.0 * params.eta * (-898433.0 + 639856.0 * params.eta))
                  + a4 * (9.0 * b2 * (8615.0 + 4.0 * params.eta * (6565.0 - 30128.0 * params.eta))
                          + (-2.0) * (605645.0 + 9.0 * params.eta * (-898433.0 + 639856.0 * params.eta)))
                  + a2 * (9.0 * b4 * (8615.0 + 4.0 * params.eta * (6565.0 - 30128.0 * params.eta))
                          + (-12.0) * (-354911.0 + 4848303.0 * params.eta + 511413.0 * eta2)
                          + (-4.0) * b2 * (605645.0 + 9.0 * params.eta * (-898433.0 + 639856.0 * params.eta))))
    eps4 = params.eps ** 4

    dp     = eps4 * (GM4 * GMp_5_2 * params.eta * term_p) / (11340.0 * p ** 4)
    dalpha = -eps4 * (alpha * params.eta * term_alpha * GM4 * GMp_5_2) / (30240.0 * p ** 5)
    dbeta  = -eps4 * (beta  * params.eta * term_alpha * GM4 * GMp_5_2) / (30240.0 * p ** 5)
    return [dp, dalpha, dbeta]


def compute_secular_rhs_TW(state: BinaryState, params: PhysicalParams,
                            max_pn_order: int) -> SecularRHS:
    total = compute_secular_rhs(state, params, min(max_pn_order, 4))
    if max_pn_order >= 5:
        rhs_tw = secular_4_5PN_TW(state, params)
        total[0] += rhs_tw[0]
        total[1] += rhs_tw[1]
        total[2] += rhs_tw[2]
    return total


def _scale_for_phi(rhs: SecularRHS, eps: float) -> SecularRHS:
    factor = 1.0 / eps if eps != 0.0 else 0.0
    return [rhs[0] * factor, rhs[1] * factor, rhs[2] * factor]


def compute_secular_rhs_phi(state: BinaryState, params: PhysicalParams,
                             max_pn_order: int) -> SecularRHS:
    return _scale_for_phi(compute_secular_rhs(state, params, max_pn_order), params.eps)


def compute_secular_rhs_TW_phi(state: BinaryState, params: PhysicalParams,
                                max_pn_order: int) -> SecularRHS:
    return _scale_for_phi(compute_secular_rhs_TW(state, params, max_pn_order), params.eps)


def transform_tilde_state_to_actual(tilde: BinaryState,
                                    params: PhysicalParams,
                                    use_TW: bool = False) -> BinaryState:
    eps2 = params.eps ** 2
    eps4 = eps2 * eps2

    Y2 = oscillatory_1PN(tilde, params)
    Y4 = oscillatory_2PN(tilde, params)

    return BinaryState(
        p     = tilde.p     + eps2 * Y2[0] + eps4 * Y4[0],
        alpha = tilde.alpha + eps2 * Y2[1] + eps4 * Y4[1],
        beta  = tilde.beta  + eps2 * Y2[2] + eps4 * Y4[2],
    )


def _approximate_tilde_from_actual(actual: BinaryState,
                                   params: PhysicalParams,
                                   use_TW: bool = False,
                                   iterations: int = 8) -> BinaryState:
    """Fixed-point inversion of the near-identity PN map actual = T(tilde)."""
    guess = BinaryState(actual.p, actual.alpha, actual.beta)
    for _ in range(iterations):
        mapped = transform_tilde_state_to_actual(guess, params, use_TW)
        guess = BinaryState(
            p=guess.p + (actual.p - mapped.p),
            alpha=guess.alpha + (actual.alpha - mapped.alpha),
            beta=guess.beta + (actual.beta - mapped.beta),
        )
    return guess


# ---------------------------------------------------------------------------
# Adaptive RK4 Integrator
# ---------------------------------------------------------------------------

@dataclass
class IntegrationResult:
    theta:           List[float]       = field(default_factory=list)
    states:          List[BinaryState] = field(default_factory=list)
    error_estimates: List[float]       = field(default_factory=list)
    num_steps:       int               = 0


def _interpolate_state_at_x(result: IntegrationResult, x: float) -> BinaryState:
    if not result.theta:
        return BinaryState()
    if x <= result.theta[0]:
        return result.states[0]
    if x >= result.theta[-1]:
        return result.states[-1]

    # binary search
    lo, hi = 0, len(result.theta) - 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if result.theta[mid] <= x:
            lo = mid
        else:
            hi = mid

    x0, x1 = result.theta[lo], result.theta[hi]
    t = (x - x0) / (x1 - x0) if x1 != x0 else 0.0
    s0, s1 = result.states[lo], result.states[hi]
    return BinaryState(
        p     = s0.p     + t * (s1.p     - s0.p),
        alpha = s0.alpha + t * (s1.alpha - s0.alpha),
        beta  = s0.beta  + t * (s1.beta  - s0.beta),
    )


RHSFunc = Callable[[BinaryState, PhysicalParams, float], SecularRHS]


class AdaptiveRK4Integrator:
    def __init__(self, tolerance: float = 1e-4,
                 min_step: float = 1e-5,
                 max_step: float = 0.1,
                 max_steps: int = 100):
        self.tolerance = tolerance
        self.min_step  = min_step
        self.max_step  = max_step
        self.max_steps = max_steps

    # -- low-level RK4 step --------------------------------------------------

    def _rk4_step(self, state: BinaryState, params: PhysicalParams,
                  rhs_func: RHSFunc, x: float, h: float) -> BinaryState:
        k1 = rhs_func(state, params, x)
        s2 = BinaryState(state.p     + 0.5*h*k1[0],
                         state.alpha + 0.5*h*k1[1],
                         state.beta  + 0.5*h*k1[2])
        k2 = rhs_func(s2, params, x + 0.5*h)
        s3 = BinaryState(state.p     + 0.5*h*k2[0],
                         state.alpha + 0.5*h*k2[1],
                         state.beta  + 0.5*h*k2[2])
        k3 = rhs_func(s3, params, x + 0.5*h)
        s4 = BinaryState(state.p     + h*k3[0],
                         state.alpha + h*k3[1],
                         state.beta  + h*k3[2])
        k4 = rhs_func(s4, params, x + h)

        return BinaryState(
            p     = state.p     + (h/6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]),
            alpha = state.alpha + (h/6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]),
            beta  = state.beta  + (h/6.0) * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]),
        )

    # -- adaptive step -------------------------------------------------------

    def adaptive_step(self, state: BinaryState, params: PhysicalParams,
                      rhs_func: RHSFunc, x: float,
                      h: float) -> Tuple[BinaryState, float]:
        y1     = self._rk4_step(state, params, rhs_func, x, h)
        y_half = self._rk4_step(state, params, rhs_func, x, h / 2.0)
        y2     = self._rk4_step(y_half, params, rhs_func, x + 0.5*h, h / 2.0)

        max_rel_error = 0.0
        for v1, v2 in [(y1.p, y2.p), (y1.alpha, y2.alpha), (y1.beta, y2.beta)]:
            if abs(v2) > 1e-15:
                max_rel_error = max(max_rel_error, abs(v1 - v2) / abs(v2))

        return y2, max_rel_error

    # -- full integration ----------------------------------------------------

    def integrate(self, initial_state: BinaryState, params: PhysicalParams,
                  rhs_func: RHSFunc,
                  theta_start: float, theta_end: float) -> IntegrationResult:
        result = IntegrationResult()
        current = initial_state
        x = theta_start
        h = self.max_step if theta_end > theta_start else -self.max_step
        step_count = 0

        result.theta.append(x)
        result.states.append(current)
        result.error_estimates.append(0.0)

        while (h > 0 and x < theta_end) or (h < 0 and x > theta_end):
            if step_count >= self.max_steps:
                break

            h_attempt = h
            if h > 0 and x + h > theta_end:
                h_attempt = theta_end - x
            elif h < 0 and x + h < theta_end:
                h_attempt = theta_end - x

            try:
                next_state, error = self.adaptive_step(current, params, rhs_func, x, h_attempt)
            except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
                # Non-physical trial states can appear during adaptive RK substeps.
                # Reject this step and force a smaller retry.
                next_state, error = current, float("inf")

            if error < self.tolerance:
                current = next_state
                x += h_attempt
                result.theta.append(x)
                result.states.append(current)
                result.error_estimates.append(error)
                step_count += 1

                scale = 0.9 * (self.tolerance / (error + 1e-16)) ** 0.2
                h *= max(0.1, min(2.0, scale))
                h = math.copysign(min(abs(h), self.max_step), h)
            else:
                h *= 0.5
                if abs(h) < self.min_step:
                    # Cannot reduce further; return the partial trajectory.
                    break

        result.num_steps = step_count
        return result


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _resolve_output_path(output_file: str) -> str:
    if not output_file:
        return output_file
    if os.path.isabs(output_file):
        return output_file
    base = os.path.dirname(os.path.abspath(sys.argv[0])) if sys.argv[0] else os.getcwd()
    return os.path.join(base, output_file)


# ---------------------------------------------------------------------------
# Fast RHS-only convergence analysis
# ---------------------------------------------------------------------------

def compute_rhs_fractional_differences(
        initial_state: BinaryState, params: PhysicalParams,
        max_pn_order: int,
        epsilon_values: List[float]
    ) -> List[Tuple[float, float, float]]:
    """
    For each epsilon, evaluate a one-step instantaneous rate probe and compute
    the fractional difference in de/dtheta in a common (physical) coordinate system:
    
        Δ(de/dtheta) = |de/dtheta_QLT − de/dtheta_Fereisl| / |de/dtheta_QLT|
        
    Returns list of tuples: (eps, frac_diff_fereisl, frac_diff_tw)
    
    INTERPRETATION (log-log plot):
    - For Fereisl vs QLT (4.5PN comparison): Expected slope +2 since the next uncaptured 
      order is 5.5PN, which is eps² above 4.5PN. The error should scale as eps².
    - For Fereisl vs TW (both 4.5PN, differ only in polynomial constants): Expected 
      slope ≈ 0 since they differ by a constant fraction independent of eps; the two 
      term_alpha polynomials differ by a tiny numerical constant, not a PN-order difference.
    """
    print("\n=== Instantaneous RHS Convergence Analysis (de/dtheta) ===")
    print(f"Initial state: p={initial_state.p}, alpha={initial_state.alpha}, beta={initial_state.beta}")
    print(f"Evaluating RHS at fixed initial state for each epsilon...")
    print(f"Observable: Δ(de/dtheta) = |de/dtheta_QLT − de/dtheta_method| / |de/dtheta_QLT|")
    print(f"On log-log plot: slope +2 indicates next missing PN order is eps² away\n")
    
    rows: List[Tuple[float, float, float]] = []
    probe_h = 1e-8
    
    for eps in epsilon_values:
        sp = PhysicalParams(params.G, params.M, params.eta, params.phi, eps)
        sp0 = PhysicalParams(params.G, params.M, params.eta, 0.0, eps)
        sp1 = PhysicalParams(params.G, params.M, params.eta, probe_h, eps)

        # Recover approximate tilde states that map to the same physical initial state.
        tilde_fereisl = _approximate_tilde_from_actual(initial_state, sp0, use_TW=False)
        tilde_tw = _approximate_tilde_from_actual(initial_state, sp0, use_TW=True)
        
        # Evaluate RHS at initial state for all three methods
        rhs_qlt = compute_qlt_rhs_phi(initial_state, sp, max_pn_order, 0.0)
        rhs_fereisl = compute_secular_rhs_phi(tilde_fereisl, sp, max_pn_order)
        rhs_tw = compute_secular_rhs_TW_phi(tilde_tw, sp, max_pn_order)

        # QLT is already physical.
        qlt_next = BinaryState(
            p=initial_state.p + probe_h * rhs_qlt[0],
            alpha=initial_state.alpha + probe_h * rhs_qlt[1],
            beta=initial_state.beta + probe_h * rhs_qlt[2],
        )
        e_dot_qlt = (_eccentricity(qlt_next) - _eccentricity(initial_state)) / probe_h

        # Fereisl/TW are in tilde coordinates; transform probe endpoints to physical first.
        fereisl_next_tilde = BinaryState(
            p=tilde_fereisl.p + probe_h * rhs_fereisl[0],
            alpha=tilde_fereisl.alpha + probe_h * rhs_fereisl[1],
            beta=tilde_fereisl.beta + probe_h * rhs_fereisl[2],
        )
        fereisl_actual_0 = transform_tilde_state_to_actual(tilde_fereisl, sp0, use_TW=False)
        fereisl_actual_1 = transform_tilde_state_to_actual(fereisl_next_tilde, sp1, use_TW=False)
        e_dot_fereisl = (_eccentricity(fereisl_actual_1) - _eccentricity(fereisl_actual_0)) / probe_h

        tw_next_tilde = BinaryState(
            p=tilde_tw.p + probe_h * rhs_tw[0],
            alpha=tilde_tw.alpha + probe_h * rhs_tw[1],
            beta=tilde_tw.beta + probe_h * rhs_tw[2],
        )
        tw_actual_0 = transform_tilde_state_to_actual(tilde_tw, sp0, use_TW=True)
        tw_actual_1 = transform_tilde_state_to_actual(tw_next_tilde, sp1, use_TW=True)
        e_dot_tw = (_eccentricity(tw_actual_1) - _eccentricity(tw_actual_0)) / probe_h
        
        # Compute fractional difference: |QLT - method| / |QLT|
        denom = abs(e_dot_qlt) + 1e-30  # avoid division by zero
        frac_diff_fereisl = abs(e_dot_qlt - e_dot_fereisl) / denom
        frac_diff_tw = abs(e_dot_qlt - e_dot_tw) / denom
        
        rows.append((eps, frac_diff_fereisl, frac_diff_tw))
        
        print(
            f"  eps={eps:.4g}:  Fereisl Δ(de/dtheta)={frac_diff_fereisl:.6e}  |  "
            f"TW Δ(de/dtheta)={frac_diff_tw:.6e}"
        )
    
    return rows


def plot_rhs_fractional_differences(rows: List[Tuple[float, float, float]], 
                                     output_image: str) -> None:
    """
    Plot fractional RHS differences vs epsilon on a log-log scale.
    
    The slope on this plot reveals the PN convergence order:
    - Slope +2: next missing order is eps² away (e.g., Fereisl 4.5PN vs QLT with 5.5PN terms)
    - Slope 0: constant fractional difference independent of eps (e.g., two implementations 
      of same PN order differing only in polynomial coefficients)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping fractional difference plot generation.")
        return
    
    out_path = _resolve_output_path(output_image)
    
    eps_vals = []
    frac_fereisl_vals = []
    frac_tw_vals = []
    
    for eps, frac_fereisl, frac_tw in rows:
        if eps > 0.0 and math.isfinite(eps) and math.isfinite(frac_fereisl) and math.isfinite(frac_tw):
            eps_vals.append(eps)
            frac_fereisl_vals.append(frac_fereisl)
            frac_tw_vals.append(frac_tw)
    
    def _log_series(y_vals: List[float]) -> Tuple[np.ndarray, np.ndarray]:
        x = []
        y = []
        for eps, y_val in zip(eps_vals, y_vals):
            if y_val > 0.0:
                x.append(math.log10(eps))
                y.append(math.log10(y_val))
        return np.array(x), np.array(y)
    
    x_fereisl, y_fereisl = _log_series(frac_fereisl_vals)
    x_tw, y_tw = _log_series(frac_tw_vals)
    
    plt.figure(figsize=(10, 6))
    
    if x_fereisl.size > 0:
        plt.plot(x_fereisl, y_fereisl, "o-", linewidth=2.5, markersize=8, 
                 label="Fereisl vs QLT (expect slope +2)")
    if x_tw.size > 0:
        plt.plot(x_tw, y_tw, "s-", linewidth=2.5, markersize=8, 
                 label="TW vs Fereisl (expect slope ≈0)")
    
    plt.xlabel("log₁₀(ε)", fontsize=12, fontweight='bold')
    plt.ylabel("log₁₀(|Δ(de/dθ)| / |de/dθ|_QLT)", fontsize=12, fontweight='bold')
    plt.title("Instantaneous RHS Convergence: Fractional Difference vs PN Expansion Parameter", 
              fontsize=13, fontweight='bold')
    plt.grid(True, alpha=0.4, linestyle='--')
    plt.legend(fontsize=11, loc='best')
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    
    print(f"\nFractional RHS difference plot written to: {out_path}")


# ---------------------------------------------------------------------------
# Comparison functions
# ---------------------------------------------------------------------------

def _map_secular_order_to_qlt_order(max_pn_order: int) -> int:
    if max_pn_order <= 2:
        return max_pn_order
    if max_pn_order <= 4:
        return 2
    return 3


def compute_qlt_rhs_phi(state: BinaryState, params: PhysicalParams,
                         max_pn_order: int, phi: float) -> SecularRHS:
    K = _get_K()
    pp      = PhysicalParams(params.G, params.M, params.eta, phi, params.eps)
    c_val   = (1.0 / pp.eps) if pp.eps != 0.0 else float("inf")
    qlt_ord = _map_secular_order_to_qlt_order(max_pn_order)
    rhs     = compute_qlt(K, state.p, state.alpha, state.beta, phi, c_val, qlt_ord)
    return [rhs.dp, rhs.dalpha, rhs.dbeta]


def _eccentricity(state: BinaryState) -> float:
    return math.sqrt(state.alpha * state.alpha + state.beta * state.beta)


def _first_phi_below_e_threshold(result: IntegrationResult, e_threshold: float) -> float:
    if not result.theta or not result.states:
        return float("nan")

    prev_phi = result.theta[0]
    prev_e = _eccentricity(result.states[0])
    if prev_e <= e_threshold:
        return prev_phi

    for i in range(1, len(result.theta)):
        curr_phi = result.theta[i]
        curr_e = _eccentricity(result.states[i])
        if curr_e <= e_threshold:
            de = curr_e - prev_e
            if de == 0.0:
                return curr_phi
            t = (e_threshold - prev_e) / de
            return prev_phi + t * (curr_phi - prev_phi)
        prev_phi, prev_e = curr_phi, curr_e

    return float("nan")


def compute_log_phi_vs_log_epsilon_data(
        initial_state: BinaryState, params: PhysicalParams,
        max_pn_order: int,
        epsilon_values: List[float],
        e_target_fraction: float,
        phi_end: float,
    tolerance: float) -> List[Tuple[float, float, float, float]]:
    """
    For each epsilon, integrate QLT/Fereisl/TW in phi and record the first phi
    where eccentricity e = sqrt(alpha^2 + beta^2) drops below a target fraction
    of the initial eccentricity.
    """
    scan_min_step = 1e-10
    scan_max_step = 1e-2
    scan_max_steps = 100_000
    integrator = AdaptiveRK4Integrator(tolerance, scan_min_step, scan_max_step, scan_max_steps)
    e0 = _eccentricity(initial_state)
    e_target = e0 * e_target_fraction
    print(
        "[eps scan config] "
        f"tol={tolerance:.1e}, min_step={scan_min_step:.1e}, "
        f"max_step={scan_max_step:.1e}, max_steps={scan_max_steps}, "
        f"e_target={e_target:.12e}"
    )

    rows: List[Tuple[float, float, float, float]] = []
    for eps in epsilon_values:
        sp = PhysicalParams(params.G, params.M, params.eta, params.phi, eps)
        sp0 = PhysicalParams(params.G, params.M, params.eta, 0.0, eps)

        tilde_fereisl_init = _approximate_tilde_from_actual(
            initial_state,
            sp0,
            use_TW=False,
        )
        tilde_tw_init = _approximate_tilde_from_actual(
            initial_state,
            sp0,
            use_TW=True,
        )

        def rhs_qlt(s, p, phi):
            return compute_qlt_rhs_phi(s, p, max_pn_order, phi)

        def rhs_fereisl(s, p, _phi):
            return compute_secular_rhs_phi(s, p, max_pn_order)

        def rhs_tw(s, p, _phi):
            return compute_secular_rhs_TW_phi(s, p, max_pn_order)

        result_qlt = integrator.integrate(initial_state, sp, rhs_qlt, 0.0, phi_end)
        result_fereisl = integrator.integrate(tilde_fereisl_init, sp, rhs_fereisl, 0.0, phi_end)
        result_tw = integrator.integrate(tilde_tw_init, sp, rhs_tw, 0.0, phi_end)

        # Extract end-state orbital elements in a common (physical) coordinate system.
        s_qlt = result_qlt.states[-1] if result_qlt.states else BinaryState()
        s_fereisl_tilde = result_fereisl.states[-1] if result_fereisl.states else BinaryState()
        s_tw_tilde = result_tw.states[-1] if result_tw.states else BinaryState()

        # Get final phi values for reference
        phi_qlt = result_qlt.theta[-1] if result_qlt.theta else phi_end
        phi_fereisl = result_fereisl.theta[-1] if result_fereisl.theta else phi_end
        phi_tw = result_tw.theta[-1] if result_tw.theta else phi_end

        # Compare all methods at a fixed orbital phase in physical coordinates.
        compare_params = PhysicalParams(params.G, params.M, params.eta, 0.0, eps)

        s_fereisl = transform_tilde_state_to_actual(
            s_fereisl_tilde,
            compare_params,
            use_TW=False,
        )
        s_tw = transform_tilde_state_to_actual(
            s_tw_tilde,
            compare_params,
            use_TW=True,
        )

        # Compute end-state differences
        dp_qf = abs(s_qlt.p - s_fereisl.p)
        da_qf = abs(s_qlt.alpha - s_fereisl.alpha)
        db_qf = abs(s_qlt.beta - s_fereisl.beta)

        dp_qt = abs(s_qlt.p - s_tw.p)
        da_qt = abs(s_qlt.alpha - s_tw.alpha)
        db_qt = abs(s_qlt.beta - s_tw.beta)

        dp_ft = abs(s_fereisl.p - s_tw.p)
        da_ft = abs(s_fereisl.alpha - s_tw.alpha)
        db_ft = abs(s_fereisl.beta - s_tw.beta)

        rows.append((eps, phi_qlt, phi_fereisl, phi_tw))

        print(
            f"[eps scan] eps={eps:.4g} | phi_end: qlt={phi_qlt:.6e}, "
            f"fereisl={phi_fereisl:.6e}, tw={phi_tw:.6e}"
        )
        print(
            f"           end-state: QLT vs Fereisl: "
            f"Δp={dp_qf:.6e}, Δα={da_qf:.6e}, Δβ={db_qf:.6e}"
        )
        print(
            f"           end-state: QLT vs TW: "
            f"Δp={dp_qt:.6e}, Δα={da_qt:.6e}, Δβ={db_qt:.6e}"
        )
        print(
            f"           end-state: Fereisl vs TW: "
            f"Δp={dp_ft:.6e}, Δα={da_ft:.6e}, Δβ={db_ft:.6e}"
        )

    return rows


def plot_log_phi_vs_log_epsilon(rows: List[Tuple[float, float, float, float]], output_image: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed; skipping log-log plot generation.")
        return

    out_path = _resolve_output_path(output_image)

    eps_vals = []
    phi_qlt_vals = []
    phi_fereisl_vals = []
    phi_tw_vals = []

    for eps, p_qlt, p_fereisl, p_tw in rows:
        eps_vals.append(eps)
        phi_qlt_vals.append(p_qlt)
        phi_fereisl_vals.append(p_fereisl)
        phi_tw_vals.append(p_tw)

    def _log_series(y_vals: List[float]) -> Tuple[np.ndarray, np.ndarray]:
        x = []
        y = []
        for eps, phi in zip(eps_vals, y_vals):
            if eps > 0.0 and phi > 0.0 and math.isfinite(phi):
                x.append(math.log10(eps))
                y.append(math.log10(phi))
        return np.array(x), np.array(y)

    x_q, y_q = _log_series(phi_qlt_vals)
    x_f, y_f = _log_series(phi_fereisl_vals)
    x_t, y_t = _log_series(phi_tw_vals)

    plt.figure(figsize=(8, 5))
    if x_q.size:
        plt.plot(x_q, y_q, "o-", label="QLT")
    if x_f.size:
        plt.plot(x_f, y_f, "s-", label="Fereisl")
    if x_t.size:
        plt.plot(x_t, y_t, "^-", label="TW")

    plt.xlabel("log10(epsilon)")
    plt.ylabel("log10(phi)")
    plt.title("log(phi) vs log(epsilon)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    print(f"log(phi) vs log(epsilon) plot written to: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    initial_state = BinaryState(p=10.0, alpha=0.1, beta=0.1)
    params        = PhysicalParams(G=1.0, M=1.0, eta=0.25, phi=0.0, eps=1.0)

    theta_start  = 0.0
    theta_end    = 10.0
    tolerance    = 1e-9
    max_pn_order = 5

    print("=== Adaptive RK4 Integration Comparison ===")
    print(f"PN order: {max_pn_order} (4.5PN)")
    print(f"Tolerance: {tolerance}")
    print(f"Integration: theta = [{theta_start}, {theta_end}]")

    phi_start = 0.0
    phi_end   = params.eps * theta_end
    print(f"Phi-transformed integration: phi = [{phi_start}, {phi_end}]")

    # PNG-only mode: skip CSV-producing routines.

    eps_scan_values = [0.032, 0.016, 0.008, 0.004, 0.002]
    
    # Fast RHS-only convergence analysis (milliseconds)
    rhs_frac_rows = compute_rhs_fractional_differences(
        initial_state=initial_state,
        params=params,
        max_pn_order=max_pn_order,
        epsilon_values=eps_scan_values,
    )
    plot_rhs_fractional_differences(rhs_frac_rows, "rhs_fractional_differences.png")
    
    # Full integration epsilon scan (slow, ~20 minutes)
    eps_scan_phi_end = 100.0
    print(f"\nEpsilon scan phi_end: {eps_scan_phi_end}")
    eps_phi_rows = compute_log_phi_vs_log_epsilon_data(
        initial_state=initial_state,
        params=params,
        max_pn_order=max_pn_order,
        epsilon_values=eps_scan_values,
        e_target_fraction=0.95,
        phi_end=eps_scan_phi_end,
        tolerance=tolerance,
    )
    plot_log_phi_vs_log_epsilon(eps_phi_rows, "log_phi_vs_log_epsilon_frompy_rk4.png")


if __name__ == "__main__":
    main()
