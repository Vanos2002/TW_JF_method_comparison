4.5PN Secular Evolution: Feireisl vs Tucker-Will Comparison
Overview
This code compares two independent derivations of the 4.5 post-Newtonian (PN) secular orbital evolution equations for inspiralling binary systems: the Feireisl expressions and the Tucker-Will (TW) expressions from arXiv:2108.12210. Both derive secular equations for the orbit-averaged evolution of the quasi-Keplerian orbital elements — semi-latus rectum pp
p and eccentricity vector (α,β)(\alpha, \beta)
(α,β) — under gravitational-wave radiation reaction up to and including 4.5PN order.
Physical Setup
The orbit is parametrized in quasi-Keplerian form by the triplet (p,α,β)(p, \alpha, \beta)
(p,α,β), where e=α2+β2e = \sqrt{\alpha^2 + \beta^2}
e=α2+β2​ is the eccentricity. The secular equations govern the slow evolution of the tilde variables (p~,α~,β~)(\tilde{p}, \tilde{\alpha}, \tilde{\beta})
(p~​,α~,β~​) — the orbital elements with fast orbital oscillations removed via a near-identity transformation:
p(φ)=p~(φ)+ε2Y2p+ε4Y4pp(\varphi) = \tilde{p}(\varphi) + \varepsilon^2 Y_2^p + \varepsilon^4 Y_4^pp(φ)=p~​(φ)+ε2Y2p​+ε4Y4p​
where ε=1/c\varepsilon = 1/c
ε=1/c is the PN expansion parameter and Y2Y_2
Y2​, Y4Y_4
Y4​ are explicit Fourier series in the true anomaly φ\varphi
φ encoding the 1PN and 2PN oscillatory corrections.
The reference solution is provided by the quasi-Lie transform (QLT) equations of motion, which integrate the full instantaneous PN acceleration (conservative 1PN–3PN plus dissipative 2.5PN–4.5PN) orbit by orbit in φ\varphi
φ, without any averaging.
What is Being Tested
The central question is whether the Feireisl 4.5PN secular expressions converge to the QLT reference more accurately than the TW expressions, and whether the two implementations are mutually consistent.
The comparison is structured in two layers:
Convergence to QLT. For each value of ε\varepsilon
ε, all three systems (QLT, Feireisl, TW) are integrated from the same physical initial conditions (p0,α0,β0)(p_0, \alpha_0, \beta_0)
(p0​,α0​,β0​) to a fixed φend\varphi_\mathrm{end}
φend​. The Feireisl and TW tilde end-states are mapped back to physical coordinates via the near-identity reconstruction before differencing against QLT. The error ∣Δp(φend)∣|\Delta p(\varphi_\mathrm{end})|
∣Δp(φend​)∣ is expected to scale as ε2\varepsilon^2
ε2 if the 4.5PN secular equations are correctly capturing the leading dissipative evolution, since the next missing order is 5.5PN.
Feireisl vs TW consistency. The difference between the two 4.5PN implementations is measured directly. If both correctly implement 4.5PN, their difference should appear at 5.5PN order, i.e. scale as ε5\varepsilon^5
ε5 relative to the leading dissipative term.
Numerical Methods
To ensure results are independent of the integration scheme, two entirely independent integrators are used:

Adaptive RK4 — classical fourth-order Runge-Kutta with step-doubling error control and step-size rescaling via (0.9⋅(tol/err)0.2)(0.9 \cdot (\mathrm{tol}/\mathrm{err})^{0.2})
(0.9⋅(tol/err)0.2)
Implicit Gauss collocation — a high-order implicit method with superior stability properties for stiff or oscillatory systems

Both integrators use tolerance 10−1410^{-14}
10−14 and produce bit-for-bit identical end-states, confirming that all observed differences are physical rather than numerical artifacts.
Results
The end-state convergence plot shows three curves of ∣Δp(φend)∣|\Delta p(\varphi_\mathrm{end})|
∣Δp(φend​)∣ vs ε\varepsilon
ε:

QLT vs Feireisl and QLT vs TW both follow ε2\varepsilon^2
ε2 scaling, confirming correct 4.5PN implementation in both cases
Feireisl vs TW follows ε5\varepsilon^5
ε5 scaling, confirming the two expressions agree at 4.5PN and differ only in their 5.5PN remainder — a consequence of different resummation choices in the two derivations
The two upper curves are indistinguishable from each other at all tested ε\varepsilon
ε values, meaning neither implementation can be shown to converge to QLT more accurately than the other within the tested parameter range
