# 4.5PN Secular Evolution: Feireisl vs Tucker-Will

Numerical verification of the 4.5 post-Newtonian secular radiation-reaction equations for inspiralling eccentric binaries, comparing two independent derivations against a reference quasi-Lie transform (QLT) integration.

---

## Background

At 4.5PN order, two independent derivations exist for the orbit-averaged (secular) evolution of the orbital elements under gravitational-wave radiation reaction:

- **Feireisl** — the expressions being validated in this work
- **Tucker-Will (TW)** — from Tucker & Will, [arXiv:2108.12210](https://arxiv.org/abs/2108.12210)

Both give secular equations for $d\tilde{p}/d\theta$, $d\tilde{\alpha}/d\theta$, $d\tilde{\beta}/d\theta$, where $(\tilde{p}, \tilde{\alpha}, \tilde{\beta})$ are the **tilde (secular) orbital elements** — the quasi-Keplerian orbital parameters with fast orbital oscillations removed. The physical orbital elements are recovered via the near-identity transformation

$$p(\varphi) = \tilde{p}(\varphi) + \varepsilon^2\, Y_2^p(\tilde{p},\tilde{\alpha},\tilde{\beta},\varphi) + \varepsilon^4\, Y_4^p(\tilde{p},\tilde{\alpha},\tilde{\beta},\varphi)$$

and analogously for $\alpha$ and $\beta$, where $\varepsilon = 1/c$ is the PN expansion parameter and $Y_2$, $Y_4$ are explicit Fourier series in the true anomaly $\varphi$ encoding the 1PN and 2PN oscillatory corrections.

---

## Orbital Parametrization

The orbit is described by the triplet $(p, \alpha, \beta)$, where

$$e = \sqrt{\alpha^2 + \beta^2}$$

is the eccentricity, related to the standard quasi-Keplerian elements by $\alpha = e\cos\omega$, $\beta = e\sin\omega$. The instantaneous radius is

$$r(\varphi) = \frac{p}{1 + e\cos\varphi}$$

---

## Reference Solution: QLT

The reference against which both secular methods are tested is the **quasi-Lie transform (QLT)** integration, which evolves the full instantaneous PN equations of motion in $\varphi$ without any orbit-averaging:

- Conservative: 1PN + 2PN + 3PN
- Dissipative: 2.5PN + 3.5PN + 4.5PN

This gives the physical $(p, \alpha, \beta)$ directly at every $\varphi$, with no near-identity transformation needed.

---

## What is Tested

Two questions are addressed:

**1. Do Feireisl and TW both correctly implement 4.5PN?**

All three systems (QLT, Feireisl, TW) are integrated from identical physical initial conditions $(p_0, \alpha_0, \beta_0)$ to a fixed $\varphi_\mathrm{end}$. The Feireisl and TW tilde end-states are mapped back to physical coordinates via the near-identity reconstruction before differencing against QLT. If 4.5PN is correctly captured, the error should scale as

$$|\Delta p(\varphi_\mathrm{end})| \propto \varepsilon^2$$

since the next missing order is 5.5PN.

**2. Are Feireisl and TW mutually consistent at 4.5PN?**

The difference between the two implementations is measured directly. If both correctly implement 4.5PN and differ only in their treatment of higher-order remainders, the Feireisl–TW difference should scale as

$$|\Delta p_\mathrm{Feireisl-TW}| \propto \varepsilon^5$$

---

## Numerical Methods

To ensure all results are independent of the integration scheme, two entirely independent integrators are used and their outputs compared:

| Method | Type | Order | Step control | Step size limits | Max steps |
|---|---|---|---|---|---|
| Adaptive RK4 | Explicit | 4 | Step-doubling, scale $(0.9 \cdot (\mathrm{tol}/\mathrm{err})^{0.2})$ | $[10^{-13},\ 10^{-3}]$ | $10^7$ |
| Implicit Gauss collocation | Implicit | High | Fixed tolerance | $[10^{-13},\ 10^{-3}]$ | $10^7$ |


Both integrators use a tolerance $10^{-14}$. Their end-states agree to machine precision across all $\varepsilon$ values tested, confirming that all observed differences between methods are physical, not numerical.

---
## Results
TW is marginally closer to QLT by 4.4% at ε = 0.032, but this advantage is ε³-suppressed and becomes completely negligible at physically relevant ε values — at ε = 0.004 (already a modest PN regime) the difference is below 0.009%. 
The figure below shows $|\Delta p(\varphi_\mathrm{end})|$ vs $\varepsilon$ for the three pairwise comparisons.
<img width="1009" height="771" alt="Screenshot 2026-06-30 at 13 20 57" src="https://github.com/user-attachments/assets/c8e42fc3-7641-402b-942b-b6059ce68f93" />

## Numerical Comparison: Feireisl vs Tucker-Will at 4.5PN

Adaptive RK4 (tol = 1e-14) integration of the orbit-averaged secular equations over φ ∈ [0, 100], initial state p = 10, α = 0.1, β = 0.1, five ε values ∈ {0.032, 0.016, 0.008, 0.004, 0.002}.

### Convergence to QLT

| Comparison | Variable | Log-log slope | Expected |
|---|---|---|---|
| QLT vs Feireisl | p | 2.001 | ε² |
| QLT vs TW | p | 1.987 | ε² |
| QLT vs Feireisl | α, β | 2.25 – 2.39 | ε² (approached from above) |
| Feireisl vs TW | p | 4.82 | ε⁵ |
| Feireisl vs TW | α, β | 5.11 | ε⁵ |

Both Feireisl and TW converge to QLT with the same asymptotic ε² rate and nearly identical prefactors, confirming correct independent implementations of the 4.5PN secular equations.

The Feireisl–TW mutual difference scales as ε⁵ (doubling ratio ≈ 32 = 2⁵ across three clean steps), indicating they first disagree at 5PN order — consistent with identical 4.5PN content but different resummation choices generating distinct 5PN remainders.

### Which method is closer to QLT

The differences are tiny and inconsistent in direction:

| Variable | Closer method | Margin at ε = 2×10⁻³ |
|---|---|---|
| p | TW | ~1.000019x (0.0019%) |
| α | Feireisl | ~1.0000017x (0.00017%) |
| β | TW | ~1.0000009x (0.000086%) |

### Practical magnitude

| ε | TW closer to QLT (in p) by |
|---|---|
| 0.032 | 4.40% |
| 0.016 | 0.55% |
| 0.008 | 0.069% |
| 0.004 | 0.009% |
| 0.002 | 0.002% |

TW has a marginal advantage in p at large ε, but the gap is ε³-suppressed (drops by ≈ 8³ = 512× over the tested range) and is physically negligible for ε ≲ 0.01. For α and β the two methods are effectively tied across all tested values.

**Conclusion:** Feireisl and TW are numerically equivalent implementations of the 4.5PN secular equations. Any preference between them is irrelevant at physically realistic ε values.

---

## Initial Conditions

| Parameter | Value |
|---|---|
| $p_0$ | $10.0$ |
| $\alpha_0$ | $0.1$ |
| $\beta_0$ | $0.1$ |
| $e_0$ | $0.1\sqrt{2}$ |
| $\eta$ | $0.25$ (equal mass) |
| $\varphi_\mathrm{end}$ | $100$ |
| $\varepsilon$ scan | $\{0.032,\ 0.016,\ 0.008,\ 0.004,\ 0.002\}$ |

Feireisl and TW are initialized from the tilde variables obtained by inverting the near-identity map at each $\varepsilon$, so that all three methods share the same physical initial conditions.

---

## Reference

Tucker, A. & Will, C. M., *Residual eccentricity of inspiralling orbits at the gravitational-wave detection threshold*, [arXiv:2108.12210](https://arxiv.org/abs/2108.12210) (2021).
