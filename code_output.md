=== Adaptive RK4 Integration Comparison ===
PN order: 5 (4.5PN)
Tolerance: 1e-14
Integration: theta = [0, 10]
Phi-transformed integration: phi = [0, 10]

=== Instantaneous RHS Convergence Analysis (de/dphi) ===
Initial state: p=10, alpha=0.1, beta=0.1
Evaluating RHS at fixed initial state for each epsilon...
Observable: |de/dphi_QLT - de/dphi_method| / |de/dphi_QLT|
  eps=3.2000e-02: Feireisl frac=3.240191e-01 | TW frac=3.240468e-01
  eps=1.6000e-02: Feireisl frac=3.237546e-01 | TW frac=3.237546e-01
  eps=8.0000e-03: Feireisl frac=3.237187e-01 | TW frac=3.237187e-01
  eps=4.0000e-03: Feireisl frac=3.238687e-01 | TW frac=3.238687e-01
  eps=2.0000e-03: Feireisl frac=3.226950e-01 | TW frac=3.226950e-01
RHS fractional-difference plot written to: /Users/matejvana/PycharmProjects/PythonProject/.venv/bin/bbh_animations/4_5pn_comparison/rhs_fractional_differences_fromcpp_rk4.png
Epsilon scan phi_end: 100
[eps scan config] tol=1.0e-14, min_step=1.0e-13, max_step=1.0e-03, max_steps=10000000
[eps scan] eps=3.2000e-02 | phi_end: qlt=1.0000e+02, feireisl=1.0000e+02, tw=1.0000e+02
           end-state: QLT vs Feireisl: dp=4.633369e-04, da=1.384253e-04, db=1.814154e-04
           end-state: QLT vs TW: dp=4.429408e-04, da=1.388571e-04, db=1.809852e-04
           end-state: Feireisl vs TW: dp=2.039609e-05, da=4.317704e-07, db=4.302742e-07
[eps scan] eps=1.6000e-02 | phi_end: qlt=1.0000e+02, feireisl=1.0000e+02, tw=1.0000e+02
           end-state: QLT vs Feireisl: dp=1.154780e-04, da=2.235341e-05, db=3.298447e-05
           end-state: QLT vs TW: dp=1.148403e-04, da=2.236687e-05, db=3.297102e-05
           end-state: Feireisl vs TW: dp=6.377035e-07, da=1.346069e-08, db=1.344647e-08
[eps scan] eps=8.0000e-03 | phi_end: qlt=1.0000e+02, feireisl=1.0000e+02, tw=1.0000e+02
           end-state: QLT vs Feireisl: dp=2.885598e-05, da=4.054469e-06, db=6.706328e-06
           end-state: QLT vs TW: dp=2.883609e-05, da=4.054889e-06, db=6.705907e-06
           end-state: Feireisl vs TW: dp=1.989513e-08, da=4.203204e-10, db=4.203042e-10
[eps scan] eps=4.0000e-03 | phi_end: qlt=1.0000e+02, feireisl=1.0000e+02, tw=1.0000e+02
           end-state: QLT vs Feireisl: dp=7.213427e-06, da=8.217392e-07, db=1.484394e-06
           end-state: QLT vs TW: dp=7.212716e-06, da=8.217531e-07, db=1.484380e-06
           end-state: Feireisl vs TW: dp=7.105410e-10, da=1.387807e-11, db=1.387765e-11
[eps scan] eps=2.0000e-03 | phi_end: qlt=1.0000e+02, feireisl=1.0000e+02, tw=1.0000e+02
           end-state: QLT vs Feireisl: dp=1.803330e-06, da=1.814421e-07, db=3.470884e-07
           end-state: QLT vs TW: dp=1.803330e-06, da=1.814421e-07, db=3.470884e-07
           end-state: Feireisl vs TW: dp=0.000000e+00, da=0.000000e+00, db=0.000000e+00
