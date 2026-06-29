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

| Method | Type | Order | Step control |
|---|---|---|---|
| Adaptive RK4 | Explicit | 4 | Step-doubling, scale $(0.9 \cdot (\mathrm{tol}/\mathrm{err})^{0.2})$ |
| Implicit Gauss collocation | Implicit | High | Fixed tolerance |

Both integrators use tolerance $10^{-14}$. Their end-states agree to machine precision across all $\varepsilon$ values tested, confirming that all observed differences between methods are physical, not numerical.
| Method | Type | Order | Step control | Step size limits | Max steps |
|---|---|---|---|---|---|
| Adaptive RK4 | Explicit | 4 | Step-doubling, scale $(0.9 \cdot (\mathrm{tol}/\mathrm{err})^{0.2})$ | $[10^{-10},\ 10^{-2}]$ | $10^5$ |
| Implicit Gauss collocation | Implicit | High | Fixed tolerance | $[10^{-10},\ 10^{-2}]$ | $10^5$ |

---

## Results

The figure below shows $|\Delta p(\varphi_\mathrm{end})|$ vs $\varepsilon$ for the three pairwise comparisons.
<img width="1088" height="840" alt="Screenshot 2026-06-29 at 11 54 50" src="https://github.com/user-attachments/assets/d097e7fa-c866-4a78-91d2-616512a6129b" />



| Comparison | Observed scaling | Interpretation |
|---|---|---|
| QLT vs Feireisl | $\varepsilon^2$ | Correct 4.5PN implementation ✓ |
| QLT vs TW | $\varepsilon^2$ | Correct 4.5PN implementation ✓ |
| Feireisl vs TW | $\varepsilon^5$ | Agree at 4.5PN, differ at 5.5PN |

**Key finding:** the QLT vs Feireisl and QLT vs TW curves are indistinguishable at all tested $\varepsilon$ values. Neither implementation converges to QLT more accurately than the other within the tested parameter range. The $\varepsilon^5$ scaling of the Feireisl–TW difference confirms that the two expressions encode different 5.5PN remainders — a consequence of different resummation choices in the two derivations — but are equivalent at 4.5PN order.

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
