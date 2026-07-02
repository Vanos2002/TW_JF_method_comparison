// This file is called "what.cpp" in VS Code
// THIS CODE WAS CREATED ON THURSDAY JUNE 4, 2026
// THIS FILE COMPARES THE Feireisl-TUCKERWILL 4.5PN EXPRESSIONS TO THE RK4 ADAPTIVE CODE

// From the article "Residual eccentricity of inspiralling orbits at the gravitational-wave detection threshold: Accurate estimates using post-Newtonian theory"
// by Alexandria Tucker and Clifford M. Will (arXiv:2108.12210v2 [gr-qc] 15 Nov 2021)
// We compare the transformed 4.5PN contributions of dp/dtheta and de/dtheta from the Feireisl-Tucker-Will (FTW) paper against the numerical orbit-averaged QLT results

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <functional>
#include <iomanip>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdio>
#include <unistd.h>
#include <limits.h>
#include <chrono>
#include <sys/wait.h>
using namespace std;

// Constants (geometrized units)
constexpr double G = 1.0;               // Gravitational constant in geometrized units
constexpr double c = 1.0;               // Speed of light in geometrized units    

constexpr double M_sun = 1.98847e30;                // 1 Solar mass in kg
constexpr double G_SI = 6.6743e-11;                 // Gravitational constant in m^3/kg/s^2
constexpr double c_SI = 299792458;                  // Speed of light in m/s
constexpr double PI = 3.14159265358979323846;       // M_PI is not constexpr in <cmath>, so define our own
constexpr double PI2 = PI * PI;                     // π²

// Masses
constexpr double m1_solar = 1.0;                        // Mass of black hole 1 in Solar masses
constexpr double m2_solar = 1.0;                        // Mass of black hole 2 in Solar masses
constexpr double M_total_solar = m1_solar + m2_solar;   // Total mass of the system in Solar masses

constexpr double m1 = m1_solar / M_total_solar;                 // Mass of black hole 1
constexpr double m2 = m2_solar / M_total_solar;                 // Mass of black hole 2
constexpr double M = m1 + m2;                                   // Total mass of the system, defined as M=1
constexpr double mu = m1 * m2 / M;                              // Reduced mass in code units   
constexpr double eta = mu / M;                                   // Symmetric mass ratio in code units
const double Mc = std::pow(m1 * m2, 0.6) * std::pow(M, -0.2);         // Chirp mass in code units
constexpr double M_kg = M_total_solar * M_sun;                  // Total mass in kg

constexpr double time_unit_seconds = G_SI * M_kg / (c_SI * c_SI * c_SI);  // 1 code unit in seconds
constexpr double sep_unit_meters = G_SI * M_kg / (c_SI * c_SI);           // 1 code unit in meters

// PN coefficient tables
struct PNCoeffs {
    using Table = std::array<std::array<std::array<double,4>,4>,4>;
    Table a{}, b{}, c{}, d{};
    double A(int l,int m,int n) const { return a[l][m][n]; }
    double B(int l,int m,int n) const { return b[l][m][n]; }
    double C(int l,int m,int n) const { return c[l][m][n]; }
    double D(int l,int m,int n) const { return d[l][m][n]; }
};

// Defining the PN coefficients (both conservative and dissipative)
PNCoeffs buildCoefficients(const bool include_4p5PN = true)
{
    const double e2 = eta * eta;
    const double e3 = e2  * eta;

    PNCoeffs K{};

    // ── 1 PN ──────────────────────────────────────────────────────────────
    K.a[1][0][0] =  2.0*(2.0 + eta);
    K.a[0][1][0] =  1.5*eta;
    K.a[0][0][1] = -(1.0 + 3.0*eta);

    K.b[1][0][0] = 0.0;
    K.b[0][1][0] =  2.0*(2.0 - eta);
    K.b[0][0][1] =  0.0;

    // ── 2 PN ──────────────────────────────────────────────────────────────
    K.a[2][0][0] = -0.75*(12.0 + 29.0*eta);
    K.a[0][2][0] = -1.875*eta*(1.0 - 3.0*eta);
    K.a[0][0][2] = -eta*(3.0 - 4.0*eta);
    K.a[1][1][0] =  2.0 + 25.0*eta + 2.0*e2;
    K.a[1][0][1] =  0.5*eta*(13.0 - 4.0*eta);
    K.a[0][1][1] =  1.5*eta*(3.0 - 4.0*eta);

    K.b[2][0][0] = 0.0;
    K.b[0][2][0] = -1.5*eta*(3.0 + 2.0*eta);
    K.b[0][0][2] = 0.0;
    K.b[1][1][0] = -0.5*(4.0 + 41.0*eta + 8.0*e2);
    K.b[1][0][1] = 0.0;
    K.b[0][1][1] =  0.5*eta*(15.0 + 4.0*eta);

    // ── 3 PN ──────────────────────────────────────────────────────────────
    K.a[3][0][0] =  16.0 + eta*(5596.0 - 123.0*PI2 + 1704.0*eta)/48.0;
    K.a[0][3][0] =  (35.0/16.0)*eta*(1.0 - 5.0*eta + 5.0*e2);
    K.a[0][0][3] = -(1.0/4.0)*eta*(11.0 - 49.0*eta + 52.0*e2);
    K.a[2][1][0] = -1.0 - (22717.0/168.0 + 615.0/64.0*PI2)*eta
                        - (11.0/8.0)*e2 + 7.0*e3;
    K.a[2][0][1] =  eta*(20827.0/840.0 + 123.0/64.0*PI2 - e2);
    K.a[1][2][0] = -0.5*eta*(158.0 - 69.0*eta - 60.0*e2);
    K.a[1][1][1] =  eta*(121.0 - 16.0*eta - 20.0*e2);
    K.a[1][0][2] = -(1.0/4.0)*eta*(75.0 + 32.0*eta - 40.0*e2);
    K.a[0][2][1] = -(15.0/8.0)*eta*(4.0 - 18.0*eta + 17.0*e2);
    K.a[0][1][2] =  (3.0/8.0)*eta*(20.0 - 79.0*eta + 60.0*e2);

    K.b[3][0][0] = 0.0;
    K.b[0][3][0] =  (15.0/8.0)*eta*(3.0 - 8.0*eta - 2.0*e2);
    K.b[0][0][3] = 0.0;
    K.b[2][1][0] =  4.0 + (5849.0/840.0 + 123.0/32.0*PI2)*eta
                        - 25.0*e2 - 8.0*e3;
    K.b[2][0][1] = 0.0;
    K.b[1][2][0] = -(1.0/6.0)*eta*(329.0 + 177.0*eta + 108.0*e2);
    K.b[1][1][1] =  eta*(15.0 + 27.0*eta + 10.0*e2);
    K.b[1][0][2] = 0.0;
    K.b[0][2][1] = -(3.0/4.0)*eta*(16.0 - 37.0*eta - 16.0*e2);
    K.b[0][1][2] =  (1.0/8.0)*eta*(65.0 - 152.0*eta - 48.0*e2);

    // ── 2.5 PN  (radiation reaction) ──────────────────────────────────────
    K.c[1][0][0] =  17.0/3.0;
    K.c[0][1][0] =  0.0;
    K.c[0][0][1] =  3.0;

    K.d[1][0][0] = -3.0;
    K.d[0][1][0] =  0.0;
    K.d[0][0][1] = -1.0;

    // ── 3.5 PN  (radiation reaction) ──────────────────────────────────────
    K.c[2][0][0] = -(23.0/14.0)*(43.0 + 14.0*eta);
    K.c[0][2][0] = -70.0;
    K.c[0][0][2] = -(3.0/28.0)*(61.0 + 70.0*eta);
    K.c[1][1][0] = -(1.0/4.0)*(147.0 + 188.0*eta);
    K.c[1][0][1] = -(1.0/42.0)*(519.0 - 1267.0*eta);
    K.c[0][1][1] =  (15.0/4.0)*(19.0 + 2.0*eta);

    K.d[2][0][0] =  (1.0/42.0)*(1325.0 + 546.0*eta);
    K.d[0][2][0] =  75.0;
    K.d[0][0][2] =  (1.0/28.0)*(313.0 + 42.0*eta);
    K.d[1][1][0] =  (1.0/12.0)*(205.0 + 424.0*eta);
    K.d[1][0][1] = -(1.0/42.0)*(205.0 + 777.0*eta);
    K.d[0][1][1] = -(3.0/4.0)*(113.0 + 2.0*eta);

    // ── 4.5 PN  (radiation reaction — disputed term) ───────────────────────
    if (include_4p5PN) {
        K.c[3][0][0] =  (1.0/756.0)*(289079.0 + 284127.0*eta + 22632.0*e2);
        K.c[0][3][0] =  0.0;
        K.c[0][0][3] =  (1.0/168.0)*(779.0 + 604.0*eta - 7090.0*e2);
        K.c[2][1][0] =  (1.0/756.0)*(250221.0 - 6032.0*eta + 74134.0*e2);
        K.c[2][0][1] = -(1.0/252.0)*(20916.0 - 24324.0*eta + 23483.0*e2);
        K.c[1][2][0] =  (1.0/252.0)*(108322.0 - 43996.0*eta + 12839.0*e2);
        K.c[1][1][1] = -(1.0/504.0)*(218401.0 - 160227.0*eta + 95987.0*e2);
        K.c[1][0][2] =  (1.0/504.0)*(40758.0 - 88311.0*eta + 43474.0*e2);
        K.c[0][2][1] =  (5.0/18.0)*(87.0 - 215.0*eta - 97.0*e2);
        K.c[0][1][2] = -(1.0/84.0)*(1205.0 - 260.0*eta - 8785.0*e2);

        K.d[3][0][0] = -(1.0/2268.0)*(395929.0 + 398700.0*eta + 87048.0*e2);
        K.d[0][3][0] =  (5.0/18.0)*(291.0 - 919.0*eta + 97.0*e2);
        K.d[0][0][3] = -(1.0/56.0)*(834.0 - 1956.0*eta - 1743.0*e2);
        K.d[2][1][0] = -(1.0/252.0)*(37992.0 + 62832.0*eta + 9649.0*e2);
        K.d[2][0][1] =  (1.0/252.0)*(26703.0 + 21304.0*eta + 28486.0*e2);
        K.d[1][2][0] = -(1.0/252.0)*(99499.0 + 24002.0*eta + 33443.0*e2);
        K.d[1][1][1] =  (1.0/504.0)*(200244.0 + 65460.0*eta + 83501.0*e2);
        K.d[1][0][2] = -(1.0/504.0)*(16731.0 + 24785.0*eta + 41471.0*e2);
        K.d[0][2][1] = -(5.0/168.0)*(6889.0 - 21631.0*eta + 2380.0*e2);
        K.d[0][1][2] =  (1.0/168.0)*(21280.0 - 60733.0*eta - 11999.0*e2);
    }

    return K;
}


// Orbital kinematics
// e = sqrt(alpha^2 + beta^2) is computed from (alpha, beta)
static double normR(double p, double e, double alpha, double beta, double phi){
    (void)e;
    return p / (1.0 + alpha*cos(phi) + beta*sin(phi));
}
static double rDot(double p, double e, double alpha, double beta, double phi){
    (void)e;
    return sqrt(G*M*p)/p * (alpha*sin(phi) - beta*cos(phi));
}
static double normV2(double p, double e, double alpha, double beta, double phi){
    double r  = normR(p,e,alpha,beta,phi);
    double rd = rDot (p,e,alpha,beta,phi);
    return rd*rd + G*M*p/(r*r);
}

// Monomial sum over one coefficient table
static double sumTable(const PNCoeffs::Table& tab, int N,
                       double rd2, double v2, double gmr, double c_val)
{
    double s = 0.0, c2N = pow(c_val, 2*N);
    for (int l=0;l<=N;++l)
        for (int m=0;m<=N-l;++m){
            int n=N-l-m; if(n<0) continue;
            s += tab[l][m][n]*pow(rd2,m)*pow(v2,n)*pow(gmr,l)/c2N;
        }
    return s;
}

// QLT equations of motion
struct QLTrhs { double dp, dalpha, dbeta; };

QLTrhs computeQLT(const PNCoeffs& K, double p,
                  double alpha, double beta, double phi,
                  double c_val, int PNorder)
{
    double e   = sqrt(alpha*alpha + beta*beta);   // e = sqrt(alpha^2 + beta^2)
    double r   = normR (p,e,alpha,beta,phi);
    double rd  = rDot  (p,e,alpha,beta,phi);
    double v2  = normV2(p,e,alpha,beta,phi);
    double rd2 = rd*rd;
    double gmr = G*M/r;
    double eps = 1.0/c_val;

    double Atot=0, Btot=0;
    int conservativeOrder = std::min(PNorder, 2);
    for (int N=1;N<=conservativeOrder;++N){
        Atot += sumTable(K.a,N,rd2,v2,gmr,c_val);
        Btot += sumTable(K.b,N,rd2,v2,gmr,c_val);
    }
    double Crr=0, Drr=0;
    for (int N=1;N<=PNorder;++N){
        Crr += sumTable(K.c,N,rd2,v2,gmr,c_val);
        Drr += sumTable(K.d,N,rd2,v2,gmr,c_val);
    }
    double eps3 = eps*eps*eps;
    double rr_R_pre = (8.0/5.0)*eta*eps3*(G*M)*(G*M)/pow(r,3)*rd;
    double rr_S_pre = (8.0/5.0)*eta*eps3*pow(G*M/(r*r),2)*sqrt(G*M*p);

    double ScR = G*M/(r*r)*(Atot+Btot) + rr_R_pre*(Crr+Drr);
    double ScS = G*M/(r*r*r)*(sqrt(G*M*p)/rd)*Btot + rr_S_pre*Drr;

    // Definition of the QLT equations of motion (see eq. 2.8 in the paper) 
    double dp_dphi = 2.0*r*r*r/(G*M)*ScS;
    double dalpha  = r*r/(G*M)*(  ScR*sin(phi) + ScS*(alpha+cos(phi))*(1.0+r/p) - ScS*alpha );
    double dbeta   = r*r/(G*M)*( -ScR*cos(phi) + ScS*(beta +sin(phi))*(1.0+r/p) - ScS*beta  );

    return {dp_dphi, dalpha, dbeta};
}



//-----------------------------------------------
// JAN Feireisl RESULTS FOR EVOLUTION
//-----------------------------------------------





// Physical constants and parameters
struct PhysicalParams {
    double G;      // Gravitational constant
    double M;      // Total mass
    double eta;    // Symmetric mass ratio (eta = m1*m2/(m1+m2)^2)
    double phi;    // Phase angle
    double eps;    // PN expansion parameter: ε = 1/c
};

// State variables for binary evolution
struct BinaryState {
    double p;        // Orbital parameter (semi-latus rectum)
    double alpha;    // Spin parameter 1
    double beta;     // Spin parameter 2
};

// Container for RHS terms: {dp/dtheta, dalpha/dtheta, dbeta/dtheta}
using SecularRHS = std::array<double, 3>;

std::array<double, 3> oscillatory_1PN(const BinaryState& state, const PhysicalParams& params);
std::array<double, 3> oscillatory_2PN(const BinaryState& state, const PhysicalParams& params);

static constexpr int PHI_AVERAGE_SAMPLES = 512;
using QLTrhsFunc = std::function<QLTrhs(const BinaryState&)>;

static QLTrhs addQLT(const QLTrhs& a, const QLTrhs& b) {
    return {a.dp + b.dp, a.dalpha + b.dalpha, a.dbeta + b.dbeta};
}

static QLTrhs subQLT(const QLTrhs& a, const QLTrhs& b) {
    return {a.dp - b.dp, a.dalpha - b.dalpha, a.dbeta - b.dbeta};
}

static QLTrhs scaleQLT(const QLTrhs& a, double s) {
    return {a.dp * s, a.dalpha * s, a.dbeta * s};
}

static QLTrhs computeQLTDissipativeOrder(const PNCoeffs& K,
                                        const BinaryState& state,
                                        const PhysicalParams& params,
                                        int dissipativeOrder)
{
    double p = state.p;
    double alpha = state.alpha;
    double beta = state.beta;
    double phi = params.phi;
    double c_val = 1.0 / params.eps;
    double e = std::sqrt(alpha*alpha + beta*beta);
    double r = normR(p, e, alpha, beta, phi);
    double rd = rDot(p, e, alpha, beta, phi);
    double v2 = normV2(p, e, alpha, beta, phi);
    double rd2 = rd * rd;
    double GM = params.G * params.M;
    double eta = params.eta;
    double gmr = GM / r;

    double Crr = sumTable(K.c, dissipativeOrder, rd2, v2, gmr, c_val);
    double Drr = sumTable(K.d, dissipativeOrder, rd2, v2, gmr, c_val);

    double eps = params.eps;
    double eps3 = eps * eps * eps;
    double rr_R_pre = (8.0 / 5.0) * eta * eps3 * GM * GM / std::pow(r, 3) * rd;
    double rr_S_pre = (8.0 / 5.0) * eta * eps3 * std::pow(GM / (r * r), 2) * std::sqrt(GM * p);

    double ScR = rr_R_pre * (Crr + Drr);
    double ScS = rr_S_pre * Drr;

    double dp_dphi = 2.0 * r * r * r / GM * ScS;
    double dalpha  = r * r / GM * (ScR * std::sin(phi) + ScS * (alpha + std::cos(phi)) * (1.0 + r / p) - ScS * alpha);
    double dbeta   = r * r / GM * (-ScR * std::cos(phi) + ScS * (beta + std::sin(phi)) * (1.0 + r / p) - ScS * beta);

    return {dp_dphi, dalpha, dbeta};
}



SecularRHS secular_1PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    double eps2 = params.eps * params.eps;
    
    SecularRHS rhs;
    rhs[0] = 0.0;
    rhs[1] = eps2 * (-(3.0 * GM * beta) / p);
    rhs[2] = eps2 * ((3.0 * GM * alpha) / p);
    
    return rhs;
}

std::array<double, 3> oscillatory_1PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double& phi = params.phi;
    const double GM = params.G * params.M;
    
    std::array<double, 3> Y;
    
    // Y[p]
    Y[0] = 4.0 * GM * (-2.0 + eta) * (alpha * std::cos(phi) + beta * std::sin(phi));
    
    // Y[alpha]
    double cos_phi = std::cos(phi);
    double cos_2phi = std::cos(2.0 * phi);
    double cos_3phi = std::cos(3.0 * phi);
    double sin_phi = std::sin(phi);
    double sin_2phi = std::sin(2.0 * phi);
    double sin_3phi = std::sin(3.0 * phi);
    
    Y[1] = (1.0 / (8.0 * p)) * GM * (
        (8.0 * (-3.0 + eta) + beta * beta * (8.0 + 21.0 * eta) + 
         alpha * alpha * (-56.0 + 47.0 * eta)) * cos_phi +
        4.0 * alpha * (-5.0 + 4.0 * eta) * cos_2phi +
        (alpha * alpha - beta * beta) * eta * cos_3phi +
        2.0 * alpha * beta * (-32.0 + 13.0 * eta) * sin_phi +
        4.0 * beta * (-5.0 + 4.0 * eta) * sin_2phi +
        2.0 * alpha * beta * eta * sin_3phi
    );
    
    // Y[beta]
    Y[2] = (1.0 / (8.0 * p)) * GM * (
        2.0 * alpha * beta * (-32.0 + 13.0 * eta) * cos_phi +
        4.0 * beta * (5.0 - 4.0 * eta) * cos_2phi -
        2.0 * alpha * beta * eta * cos_3phi +
        (8.0 * (-3.0 + eta) + alpha * alpha * (8.0 + 21.0 * eta) + 
         beta * beta * (-56.0 + 47.0 * eta)) * sin_phi +
        4.0 * alpha * (-5.0 + 4.0 * eta) * sin_2phi +
        (alpha * alpha - beta * beta) * eta * sin_3phi
    );
    
    return Y;
}

// ============================================================================
// 2 PN ORDER TERMS
// ============================================================================

SecularRHS secular_2PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM2 = GM * GM;
    double eps2 = params.eps * params.eps;
    double eps4 = eps2 * eps2;
    
    SecularRHS rhs;
    
    double bracket = -10.0 + beta * beta - 4.0 * eta + 10.0 * beta * beta * eta + 
                     alpha * alpha * (1.0 + 10.0 * eta);
    
    rhs[0] = 0.0;
    rhs[1] = eps4 * (-(3.0 * GM2 * beta * bracket) / (4.0 * p * p));
    rhs[2] = eps4 * ((3.0 * GM2 * alpha * bracket) / (4.0 * p * p));
    
    return rhs;
}

std::array<double, 3> oscillatory_2PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double& phi = params.phi;
    const double GM = params.G * params.M;
    const double GM2 = GM * GM;
    
    std::array<double, 3> Y;
    
    double cos_phi = std::cos(phi);
    double cos_2phi = std::cos(2.0 * phi);
    double cos_3phi = std::cos(3.0 * phi);
    double cos_4phi = std::cos(4.0 * phi);
    double cos_5phi = std::cos(5.0 * phi);
    double sin_phi = std::sin(phi);
    double sin_2phi = std::sin(2.0 * phi);
    double sin_3phi = std::sin(3.0 * phi);
    double sin_4phi = std::sin(4.0 * phi);
    double sin_5phi = std::sin(5.0 * phi);
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    
    // Y[p]
    Y[0] = (1.0 / (4.0 * p)) * GM2 * (
        alpha * (-160.0 + eta * (256.0 - 33.0 * a2 - 33.0 * b2 + 
                 2.0 * (-8.0 + a2 + b2) * eta)) * cos_phi +
        (a2 - b2) * (68.0 + 3.0 * eta * (-15.0 + 4.0 * eta)) * cos_2phi -
        alpha * (a2 - 3.0 * b2) * eta * (3.0 + 2.0 * eta) * cos_3phi +
        beta * (-160.0 + eta * (256.0 - 33.0 * a2 - 33.0 * b2 + 
                2.0 * (-8.0 + a2 + b2) * eta)) * sin_phi +
        2.0 * alpha * beta * (68.0 + 3.0 * eta * (-15.0 + 4.0 * eta)) * sin_2phi +
        beta * (-3.0 * a2 + b2) * eta * (3.0 + 2.0 * eta) * sin_3phi
    );
    
    // Y[alpha]
    Y[1] = (1.0 / (128.0 * p * p)) * GM2 * (
        -2.0 * (5.0 * b4 * eta * (-27.0 + 41.0 * eta) + 
                a4 * eta * (477.0 + 161.0 * eta) + 
                b2 * (-944.0 + 92.0 * eta + 64.0 * eta * eta) - 
                16.0 * (60.0 + eta * (17.0 + 8.0 * eta)) + 
                a2 * (3056.0 + 2.0 * eta * (-1958.0 + 224.0 * eta + 
                      3.0 * b2 * (57.0 + 61.0 * eta)))) * cos_phi +
        16.0 * alpha * (28.0 * (-1.0 + 5.0 * eta) + 
                       b2 * (-162.0 + eta * (57.0 + 28.0 * eta)) + 
                       a2 * (106.0 + eta * (-163.0 + 60.0 * eta))) * cos_2phi +
        (-a4 * eta * (73.0 + 53.0 * eta) + 
         b2 * (-480.0 - 8.0 * eta * (9.0 + 16.0 * eta) + b2 * eta * (-33.0 + 19.0 * eta)) + 
         2.0 * a2 * (240.0 + eta * (36.0 + 64.0 * eta + 3.0 * b2 * (53.0 + 17.0 * eta)))) * cos_3phi -
        8.0 * alpha * (a2 - 3.0 * b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * cos_4phi -
        3.0 * (a4 - 6.0 * a2 * b2 + b4) * eta * (-1.0 + 3.0 * eta) * cos_5phi +
        8.0 * alpha * beta * (-1000.0 + eta * (1002.0 - 96.0 * eta + 
                            a2 * (-153.0 + 11.0 * eta) + 
                            b2 * (-153.0 + 11.0 * eta))) * sin_phi +
        16.0 * beta * (-28.0 * (1.0 + b2 - 5.0 * eta) + 
                      b2 * eta * (-53.0 + 44.0 * eta) + 
                      a2 * (240.0 + eta * (-273.0 + 76.0 * eta))) * sin_2phi -
        4.0 * alpha * beta * (-240.0 + eta * (b2 * (-43.0 + eta) + 
                            7.0 * a2 * (9.0 + 5.0 * eta) - 
                            4.0 * (9.0 + 16.0 * eta))) * sin_3phi +
        8.0 * beta * (-3.0 * a2 + b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * sin_4phi -
        12.0 * alpha * (alpha - beta) * beta * (alpha + beta) * eta * (-1.0 + 3.0 * eta) * sin_5phi
    );
    
    // Y[beta]
    Y[2] = (1.0 / (128.0 * p * p)) * GM2 * (
        8.0 * alpha * beta * (-1000.0 + eta * (1002.0 - 96.0 * eta + 
                            a2 * (-153.0 + 11.0 * eta) + 
                            b2 * (-153.0 + 11.0 * eta))) * cos_phi -
        16.0 * beta * (28.0 * (-1.0 + 5.0 * eta) + 
                      a2 * (-162.0 + eta * (57.0 + 28.0 * eta)) + 
                      b2 * (106.0 + eta * (-163.0 + 60.0 * eta))) * cos_2phi +
        4.0 * alpha * beta * (-240.0 + eta * (a2 * (-43.0 + eta) + 
                            7.0 * b2 * (9.0 + 5.0 * eta) - 
                            4.0 * (9.0 + 16.0 * eta))) * cos_3phi -
        8.0 * beta * (-3.0 * a2 + b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * cos_4phi +
        12.0 * alpha * (alpha - beta) * beta * (alpha + beta) * eta * (-1.0 + 3.0 * eta) * cos_5phi -
        2.0 * (5.0 * a4 * eta * (-27.0 + 41.0 * eta) + 
               b4 * eta * (477.0 + 161.0 * eta) - 
               16.0 * (60.0 + eta * (17.0 + 8.0 * eta)) + 
               4.0 * b2 * (764.0 + eta * (-979.0 + 112.0 * eta)) + 
               a2 * (-944.0 + 2.0 * eta * (46.0 + 32.0 * eta + 
                     3.0 * b2 * (57.0 + 61.0 * eta)))) * sin_phi +
        16.0 * alpha * (28.0 * (-1.0 + 5.0 * eta) + 
                       a2 * (-28.0 + eta * (-53.0 + 44.0 * eta)) + 
                       b2 * (240.0 + eta * (-273.0 + 76.0 * eta))) * sin_2phi +
        (a4 * (33.0 - 19.0 * eta) * eta + 
         b2 * (-480.0 - 8.0 * eta * (9.0 + 16.0 * eta) + b2 * eta * (73.0 + 53.0 * eta)) + 
         2.0 * a2 * (240.0 + eta * (36.0 + 64.0 * eta - 3.0 * b2 * (53.0 + 17.0 * eta)))) * sin_3phi -
        8.0 * alpha * (a2 - 3.0 * b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * sin_4phi -
        3.0 * (a4 - 6.0 * a2 * b2 + b4) * eta * (-1.0 + 3.0 * eta) * sin_5phi
    );
    
    return Y;
}

// ============================================================================
// 2.5 PN ORDER TERMS (RADIATION REACTION)
// ============================================================================

SecularRHS secular_2_5PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_5_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double eps = params.eps;
    double eps5 = eps * eps * eps * eps * eps;
    
    rhs[0] = -eps5 * (8.0 * GMp_5_2 * (8.0 + 7.0 * a2 + 7.0 * b2) * eta) / (5.0 * p * p * p * p);
    rhs[1] = -eps5 * (GMp_5_2 * alpha * (304.0 + 121.0 * a2 + 121.0 * b2) * eta) / (15.0 * p * p * p * p * p);
    rhs[2] = -eps5 * (GMp_5_2 * beta * (304.0 + 121.0 * a2 + 121.0 * b2) * eta) / (15.0 * p * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// 3.5 PN ORDER TERMS
// ============================================================================

SecularRHS secular_3_5PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM3 = GM * GM * GM;
    
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_3_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    
    // dp/dtheta
    double term_p = -8.0 * (2759.0 + 252.0 * eta) + 
                    8.0 * b2 * (758.0 + 889.0 * eta) + 
                    a4 * (1483.0 + 4424.0 * eta) + 
                    b4 * (1483.0 + 4424.0 * eta) + 
                    a2 * (8.0 * (758.0 + 889.0 * eta) + b2 * (2966.0 + 8848.0 * eta));
    
    double eps = params.eps;
    double eps2 = eps * eps;
    double eps5 = eps * eps * eps * eps * eps;
    double eps7 = eps5 * eps2;
    rhs[0] = -eps7 * (GM3 * GMp_3_2 * eta * term_p) / (210.0 * p * p * p);
    
    // dalpha/dtheta
    double term_alpha = -8.0 * (18049.0 + 4452.0 * eta) + 
                       4.0 * b2 * (8692.0 + 12803.0 * eta) + 
                       a4 * (2251.0 + 15064.0 * eta) + 
                       b4 * (2251.0 + 15064.0 * eta) + 
                       a2 * (34768.0 + 51212.0 * eta + b2 * (4502.0 + 30128.0 * eta));
    
    rhs[1] = -eps7 * (GM3 * GMp_3_2 * alpha * eta * term_alpha) / (840.0 * p * p * p * p);
    
    // dbeta/dtheta
    double term_beta = -8.0 * (18049.0 + 4452.0 * eta) + 
                      4.0 * b2 * (8692.0 + 12803.0 * eta) + 
                      a4 * (2251.0 + 15064.0 * eta) + 
                      b4 * (2251.0 + 15064.0 * eta) + 
                      a2 * (34768.0 + 51212.0 * eta + b2 * (4502.0 + 30128.0 * eta));
    
    rhs[2] = -eps7 * (GM3 * GMp_3_2 * beta * eta * term_beta) / (840.0 * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// 4.5 PN ORDER TERMS
// ============================================================================

SecularRHS secular_4_5PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM4 = GM * GM * GM * GM;
    

    // Tucker-Will Results beyond
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_5_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    double a6 = a4 * a2;
    double b6 = b4 * b2;
    
    double eta2 = eta * eta;
    
    // dp/dtheta
    double term_p = 8.0 * (-1034075.0 - 261369.0 * eta + 36288.0 * eta2) + 
                   9.0 * a6 * (527.0 - 6300.0 * eta + 53088.0 * eta2) + 
                   9.0 * b6 * (527.0 - 6300.0 * eta + 53088.0 * eta2) + 
                   12.0 * b2 * (-64831.0 + 1186209.0 * eta + 64575.0 * eta2) + 
                   b4 * (947991.0 - 5386365.0 * eta + 3988656.0 * eta2) + 
                   3.0 * a4 * (315997.0 - 1795455.0 * eta + 1329552.0 * eta2 + 
                              9.0 * b2 * (527.0 - 6300.0 * eta + 53088.0 * eta2)) + 
                   3.0 * a2 * (9.0 * b4 * (527.0 - 6300.0 * eta + 53088.0 * eta2) + 
                              4.0 * (-64831.0 + 1186209.0 * eta + 64575.0 * eta2) + 
                              b2 * (631994.0 - 3590910.0 * eta + 2659104.0 * eta2));
    
    double eps = params.eps;
    double eps2 = eps * eps;
    double eps4 = eps2 * eps2;
    double eps5 = eps * eps * eps * eps * eps;
    double eps9 = eps5 * eps4;
    rhs[0] = eps9 * (GM4 * GMp_5_2 * eta * term_p) / (11340.0 * p * p * p * p);
    
    // dalpha/dtheta
    double term_alpha = 16.0 * (-2739835.0 - 1394559.0 * eta + 145152.0 * eta2) + 
                       3.0 * a6 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                       3.0 * b6 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                       12.0 * b2 * (-354911.0 + 4848903.0 * eta + 511413.0 * eta2) + 
                       2.0 * b4 * (605645.0 - 8079297.0 * eta + 5758704.0 * eta2) + 
                       a4 * (9.0 * b2 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                             2.0 * (605645.0 - 8079297.0 * eta + 5758704.0 * eta2)) + 
                       a2 * (9.0 * b4 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                             12.0 * (-354911.0 + 4848903.0 * eta + 511413.0 * eta2) + 
                             4.0 * b2 * (605645.0 - 8079297.0 * eta + 5758704.0 * eta2));
    
    rhs[1] = eps9 * (GM4 * GMp_5_2 * alpha * eta * term_alpha) / (30240.0 * p * p * p * p * p);
    
    // dbeta/dtheta (same as dalpha/dtheta but multiplied by beta instead of alpha)
    rhs[2] = eps9 * (GM4 * GMp_5_2 * beta * eta * term_alpha) / (30240.0 * p * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// COMPOSITE FUNCTIONS - COMBINING ALL PN ORDERS
// ============================================================================

SecularRHS compute_secular_RHS(const BinaryState& state, const PhysicalParams& params, int max_PN_order) {
    SecularRHS total_rhs = {0.0, 0.0, 0.0};
    
    if (max_PN_order >= 1) {
        auto rhs = secular_1PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 2) {
        auto rhs = secular_2PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 3) {  // 2.5 PN
        auto rhs = secular_2_5PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 4) {  // 3.5 PN
        auto rhs = secular_3_5PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 5) {  // 4.5 PN
        auto rhs = secular_4_5PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    return total_rhs;
}




//----------------------------------------------------------------
// TUCKER-WILL RESULTS FOR EVOLUTION (ONLY DIFFER AT 4.5 PN ORDER)
//----------------------------------------------------------------

// Tucker-Will orders 1PN through 3.5PN are identical to Feireisl - use aliases
inline SecularRHS secular_1PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_1PN(state, params);
}

inline std::array<double, 3> oscillatory_1PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return oscillatory_1PN(state, params);
}

inline SecularRHS secular_2PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_2PN(state, params);
}

inline std::array<double, 3> oscillatory_2PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return oscillatory_2PN(state, params);
}

inline SecularRHS secular_2_5PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_2_5PN(state, params);
}

inline SecularRHS secular_3_5PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_3_5PN(state, params);
}

// ============================================================================
// 4.5 PN ORDER TERMS (Tucker-Will) - Only difference from Feireisl
// ============================================================================

SecularRHS secular_4_5PN_TW(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM4 = GM * GM * GM * GM;
    
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_5_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    double a6 = a4 * a2;
    double b6 = b4 * b2;
    
    double eta2 = eta * eta;
    
    // dp/dtheta (Tucker-Will 4.5PN)
    double term_p = -8272600.0 
                  + 72.0 * eta * (-29041.0 + 4032.0 * eta)
                  + 9.0 * a6 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta))
                  + 9.0 * b6 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta))
                  + 12.0 * b2 * (-64831.0 + 9.0 * eta * (131801.0 + 7175.0 * eta))
                  + b4 * (947991.0 + 27.0 * eta * (-199495.0 + 147728.0 * eta))
                  + a4 * (9.0 * b2 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta)) + 
                         (947991.0 + 27.0 * eta * (-199495.0 + 147728.0 * eta)))
                  + a2 * (12.0 * (-64831.0 + 9.0 * eta * (131801.0 + 7175.0 * eta)) + 
                         9.0 * b4 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta)) + 
                         b2 * (2.0 * (947991.0 + 27.0 * eta * (-199495.0 + 147728.0 * eta))));
    
    double eps = params.eps;
    double eps2 = eps * eps;
    double eps4 = eps2 * eps2;
    double eps5 = eps * eps * eps * eps * eps;
    double eps9 = eps5 * eps4;
    rhs[0] = eps9 * (GM4 * GMp_5_2 * eta * term_p) / (11340.0 * p * p * p * p);
    
    // dalpha/dtheta (Tucker-Will 4.5PN)
    double term_alpha = 43837360.0 
                      + 144.0 * eta * (154951.0 - 16128.0 * eta)
                      + 9.0 * a6 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta))
                      + 9.0 * b6 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta))
                      - 12.0 * b2 * (-354911.0 + 4848303.0 * eta + 511413.0 * eta * eta)
                      - 2.0 * b4 * (605645.0 + 9.0 * eta * (-898433.0 + 639856.0 * eta))
                      + a4 * (9.0 * b2 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta)) + 
                             (-2.0) * (605645.0 + 9.0 * eta * (-898433.0 + 639856.0 * eta)))
                      + a2 * (9.0 * b4 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta)) + 
                             (-12.0) * (-354911.0 + 4848303.0 * eta + 511413.0 * eta * eta) + 
                             (-4.0) * b2 * (605645.0 + 9.0 * eta * (-898433.0 + 639856.0 * eta)));
    
    rhs[1] = -eps9 * (alpha * eta * term_alpha * GM4 * GMp_5_2) / (30240.0 * p * p * p * p * p);
    
    // dbeta/dtheta (Tucker-Will 4.5PN - same as dalpha but with beta)
    rhs[2] = -eps9 * (beta * eta * term_alpha * GM4 * GMp_5_2) / (30240.0 * p * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// COMPOSITE FUNCTIONS - COMBINING ALL PN ORDERS (Tucker-Will)
// ============================================================================

SecularRHS compute_secular_RHS_TW(const BinaryState& state, const PhysicalParams& params, int max_PN_order) {
    // Tucker-Will equations are identical to Feireisl for orders 1-3.5PN
    SecularRHS total_rhs = compute_secular_RHS(state, params, std::min(max_PN_order, 4));
    
    // Replace 4.5PN order with Tucker-Will version
    if (max_PN_order >= 5) {
        auto rhs_TW_4_5 = secular_4_5PN_TW(state, params);
        total_rhs[0] += rhs_TW_4_5[0];
        total_rhs[1] += rhs_TW_4_5[1];
            total_rhs[2] += rhs_TW_4_5[2];
    }
    
    return total_rhs;
}

static double compute_dtheta_dphi(const BinaryState& state,
                                  const PhysicalParams& params,
                                  double phi) {
    (void)state;
    (void)phi;
    return 1.0; //params.eps;
}

SecularRHS compute_feireisl_physical_RHS(const BinaryState& state,
                                        const PhysicalParams& params,
                                        int max_PN_order,
                                        double phi) {
    auto rhs_theta = compute_secular_RHS(state, params, max_PN_order);
    double dtheta_dphi = compute_dtheta_dphi(state, params, phi);
    return {
        rhs_theta[0] * dtheta_dphi,
        rhs_theta[1] * dtheta_dphi,
        rhs_theta[2] * dtheta_dphi
    };
}

SecularRHS compute_tw_physical_RHS(const BinaryState& state,
                                   const PhysicalParams& params,
                                   int max_PN_order,
                                   double phi) {
    auto rhs_theta = compute_secular_RHS_TW(state, params, max_PN_order);
    double dtheta_dphi = compute_dtheta_dphi(state, params, phi);
    return {
        rhs_theta[0] * dtheta_dphi,
        rhs_theta[1] * dtheta_dphi,
        rhs_theta[2] * dtheta_dphi
    };
}

BinaryState transformTildeStateToActual(const BinaryState& tilde_state, const PhysicalParams& params, bool useTW = false) {
    double eps2 = params.eps * params.eps;
    double eps4 = eps2 * eps2;

    auto Y2 = useTW ? oscillatory_1PN_TW(tilde_state, params) : oscillatory_1PN(tilde_state, params);
    auto Y4 = useTW ? oscillatory_2PN_TW(tilde_state, params) : oscillatory_2PN(tilde_state, params);

    return {
        tilde_state.p + eps2 * Y2[0] + eps4 * Y4[0],
        tilde_state.alpha + eps2 * Y2[1] + eps4 * Y4[1],
        tilde_state.beta + eps2 * Y2[2] + eps4 * Y4[2]
    };
}

static BinaryState approximateTildeFromActual(const BinaryState& actual,
                                              const PhysicalParams& params,
                                              bool useTW = false,
                                              int iterations = 20) {
    BinaryState guess = actual;
    for (int i = 0; i < iterations; ++i) {
        BinaryState mapped = transformTildeStateToActual(guess, params, useTW);
        const double dp = actual.p - mapped.p;
        const double da = actual.alpha - mapped.alpha;
        const double db = actual.beta - mapped.beta;
        guess = {
            guess.p + dp,
            guess.alpha + da,
            guess.beta + db
        };
        if (std::abs(dp) < 1e-13 && std::abs(da) < 1e-13 && std::abs(db) < 1e-13) {
            break;
        }
    }
    return guess;
}

// ============================================================================
// ADAPTIVE RK4 INTEGRATOR WITH ERROR ESTIMATION
// ============================================================================

struct IntegrationResult {
    std::vector<double> theta;
    std::vector<BinaryState> states;
    std::vector<double> error_estimates;
    int num_steps;
};

static bool isFiniteState(const BinaryState& s) {
    return std::isfinite(s.p) && std::isfinite(s.alpha) && std::isfinite(s.beta);
}

static BinaryState reconstructPhysicalStateFromTilde(const BinaryState& tilde_state,
                                                     const PhysicalParams& params,
                                                     double phi,
                                                     bool useTW = false) {
    PhysicalParams params_phi = params;
    params_phi.phi = phi;
    return transformTildeStateToActual(tilde_state, params_phi, useTW);
}

static IntegrationResult reconstructPhysicalTrajectoryFromTilde(const IntegrationResult& tilde_result,
                                                                const PhysicalParams& params,
                                                                bool useTW = false) {
    IntegrationResult physical_result = tilde_result;
    physical_result.states.clear();
    physical_result.states.reserve(tilde_result.states.size());
    for (size_t i = 0; i < tilde_result.states.size() && i < tilde_result.theta.size(); ++i) {
        physical_result.states.push_back(
            reconstructPhysicalStateFromTilde(tilde_result.states[i], params, tilde_result.theta[i], useTW)
        );
    }
    return physical_result;
}

static BinaryState interpolateStateAtX(const IntegrationResult& result, double x) {
    if (result.theta.empty() || result.states.empty()) {
        return {0.0, 0.0, 0.0};
    }
    if (x <= result.theta.front()) {
        return result.states.front();
    }
    if (x >= result.theta.back()) {
        return result.states.back();
    }

    auto it = std::lower_bound(result.theta.begin(), result.theta.end(), x);
    size_t hi = static_cast<size_t>(it - result.theta.begin());
    if (hi == 0) {
        return result.states.front();
    }
    size_t lo = hi - 1;

    double x0 = result.theta[lo];
    double x1 = result.theta[hi];
    double t = (x1 != x0) ? (x - x0) / (x1 - x0) : 0.0;

    const auto& s0 = result.states[lo];
    const auto& s1 = result.states[hi];
    return {
        s0.p + t * (s1.p - s0.p),
        s0.alpha + t * (s1.alpha - s0.alpha),
        s0.beta + t * (s1.beta - s0.beta)
    };
}

static std::string g_output_base_dir;

static void initializeOutputBaseDir(const char* argv0) {
    char path_buf[PATH_MAX];
    if (argv0 != nullptr && realpath(argv0, path_buf) != nullptr) {
        std::string exe_path(path_buf);
        size_t slash = exe_path.find_last_of('/');
        if (slash != std::string::npos) {
            g_output_base_dir = exe_path.substr(0, slash);
            return;
        }
    }

    char cwd_buf[PATH_MAX];
    if (getcwd(cwd_buf, sizeof(cwd_buf)) != nullptr) {
        g_output_base_dir = cwd_buf;
    } else {
        g_output_base_dir.clear();
    }
}

static std::string resolveOutputPath(const std::string& output_file) {
    if (output_file.empty()) {
        return output_file;
    }

    if (!output_file.empty() && output_file[0] == '/') {
        return output_file;
    }

    if (!g_output_base_dir.empty()) {
        return g_output_base_dir + "/" + output_file;
    }

    char cwd_buf[PATH_MAX];
    if (getcwd(cwd_buf, sizeof(cwd_buf)) != nullptr) {
        return std::string(cwd_buf) + "/" + output_file;
    }

    return output_file;
}

static std::ofstream openOutputCsv(const std::string& output_file) {
    const std::string resolved_path = resolveOutputPath(output_file);
    std::ofstream outfile(resolved_path);

    if (!outfile.is_open()) {
        std::cerr << "ERROR: Failed to open output CSV: " << resolved_path << std::endl;
    } else {
        std::cout << "Writing CSV to: " << resolved_path << std::endl;
    }

    return outfile;
}

class AdaptiveRK4Integrator {
private:
    double tolerance;
    double min_step;
    double max_step;
    int max_steps;
    
    typedef std::function<SecularRHS(const BinaryState&, const PhysicalParams&, double)> RHSFunc;
    
public:
    AdaptiveRK4Integrator(double tol = 1e-4, double h_min = 1e-5, double h_max = 0.1, int max = 100)
        : tolerance(tol), min_step(h_min), max_step(h_max), max_steps(max) {}
    
    std::pair<BinaryState, double> adaptiveStep(
        const BinaryState& state, 
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double x,
        double h
    ) {
        BinaryState y1 = rk4Step(state, params, rhs_func, x, h);
        BinaryState y_half = rk4Step(state, params, rhs_func, x, h / 2.0);
        BinaryState y2 = rk4Step(y_half, params, rhs_func, x + 0.5 * h, h / 2.0);

        if (!isFiniteState(y1) || !isFiniteState(y_half) || !isFiniteState(y2)) {
            return {state, std::numeric_limits<double>::infinity()};
        }
        
        double max_rel_error = 0.0;
        if (std::abs(y2.p) > 1e-15) {
            max_rel_error = std::max(max_rel_error, std::abs(y1.p - y2.p) / std::abs(y2.p));
        }
        if (std::abs(y2.alpha) > 1e-15) {
            max_rel_error = std::max(max_rel_error, std::abs(y1.alpha - y2.alpha) / std::abs(y2.alpha));
        }
        if (std::abs(y2.beta) > 1e-15) {
            max_rel_error = std::max(max_rel_error, std::abs(y1.beta - y2.beta) / std::abs(y2.beta));
        }

        if (!std::isfinite(max_rel_error)) {
            return {state, std::numeric_limits<double>::infinity()};
        }
        
        return {y2, max_rel_error};
    }
    IntegrationResult integrate(
        const BinaryState& initial_state,
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double theta_start,
        double theta_end
    ) {
        IntegrationResult result;
        BinaryState current_state = initial_state;
        double current_theta = theta_start;
        double h = (theta_end > theta_start) ? max_step : -max_step;
        int step_count = 0;
        
        result.theta.push_back(current_theta);
        result.states.push_back(current_state);
        result.error_estimates.push_back(0.0);
        
        while ((h > 0 && current_theta < theta_end) || (h < 0 && current_theta > theta_end)) {
            if (step_count >= max_steps) break;
            
            double h_attempt = h;
            if (h > 0 && current_theta + h > theta_end) {
                h_attempt = theta_end - current_theta;
            } else if (h < 0 && current_theta + h < theta_end) {
                h_attempt = theta_end - current_theta;
            }
            
            auto step_result = adaptiveStep(current_state, params, rhs_func, current_theta, h_attempt);
            BinaryState next_state = step_result.first;
            double error = step_result.second;
            
            if (error < tolerance) {
                current_state = next_state;
                current_theta += h_attempt;
                result.theta.push_back(current_theta);
                result.states.push_back(current_state);
                result.error_estimates.push_back(error);
                step_count++;
                
                double scale = 0.9 * std::pow(tolerance / (error + 1e-16), 0.2);
                h *= std::max(0.1, std::min(2.0, scale));
                h = std::copysign(std::min(std::abs(h), max_step), h);
            } else {
                h *= 0.5;
                if (std::abs(h) < min_step) h = std::copysign(min_step, h);
            }
        }
        
        result.num_steps = step_count;
        return result;
    }

    // Integrate until p_physical (returned by get_p) drops to p_final.
    // Returns the interpolated phi at the crossing, or the last reached phi if not reached.
    // get_p(state, phi) : evaluates the physical p from the current integration state
    //   – for QLT (physical space): just returns state.p
    //   – for tilde-space methods: reconstructs p via the oscillatory transform
    double integrateUntilP(
        const BinaryState& initial_state,
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double phi_start,
        double p_final,
        std::function<double(const BinaryState&, const PhysicalParams&, double)> get_p,
        bool* reached_target = nullptr
    ) {
        BinaryState cur  = initial_state;
        double      phi  = phi_start;
        double      h    = max_step;
        int         steps = 0;

        double p_cur = get_p(cur, params, phi);

        while (true) {
            if (steps >= max_steps) {
                break;
            }

            double h_attempt = h;

            auto step_pair = adaptiveStep(cur, params, rhs_func, phi, h_attempt);
            BinaryState next = step_pair.first;
            double error = step_pair.second;

            if (std::isfinite(error) && error < tolerance) {
                double p_next = get_p(next, params, phi + h_attempt);
                if (!std::isfinite(p_next)) {
                    error = std::numeric_limits<double>::infinity();
                } else {
                    // Check for crossing
                    if (p_next <= p_final && p_cur > p_final) {
                        // Linear interpolation: find t in [0,1] where p = p_final
                        double t = (p_cur - p_final) / (p_cur - p_next);
                        if (reached_target != nullptr) {
                            *reached_target = true;
                        }
                        return phi + t * h_attempt;
                    }

                    cur   = next;
                    phi  += h_attempt;
                    p_cur = p_next;
                    steps++;

                    double scale = 0.9 * std::pow(tolerance / (error + 1e-16), 0.2);
                    h *= std::max(0.1, std::min(2.0, scale));
                    h  = std::min(h, max_step);
                    continue;
                }
            }

            h *= 0.5;
            if (h < min_step) {
                if (reached_target != nullptr) {
                    *reached_target = false;
                }
                return phi;
            }
        }
        if (reached_target != nullptr) {
            *reached_target = false;
        }
        return phi; // p_final not reached; return last physically integrated phi
    }

private:
    BinaryState rk4Step(
        const BinaryState& state,
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double x,
        double h
    ) {
        auto k1 = rhs_func(state, params, x);
        BinaryState state2{state.p + 0.5*h*k1[0], state.alpha + 0.5*h*k1[1], state.beta + 0.5*h*k1[2]};
        auto k2 = rhs_func(state2, params, x + 0.5 * h);
        BinaryState state3{state.p + 0.5*h*k2[0], state.alpha + 0.5*h*k2[1], state.beta + 0.5*h*k2[2]};
        auto k3 = rhs_func(state3, params, x + 0.5 * h);
        BinaryState state4{state.p + h*k3[0], state.alpha + h*k3[1], state.beta + h*k3[2]};
        auto k4 = rhs_func(state4, params, x + h);
        
        BinaryState next_state;
        next_state.p = state.p + (h/6.0) * (k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
        next_state.alpha = state.alpha + (h/6.0) * (k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
        next_state.beta = state.beta + (h/6.0) * (k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]);
        return next_state;
    }
};

// ============================================================================
// INTEGRATION COMPARISON
// ============================================================================

void compareEvolutionMethods(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double theta_start,
    double theta_end,
    double tolerance,
    const std::string& output_file
) {
    AdaptiveRK4Integrator integrator(tolerance);
    
    auto rhs_feireisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double /*theta*/) {
        return compute_secular_RHS(s, p, max_PN_order);
    };
    auto result_feireisl = integrator.integrate(initial_state, params, rhs_feireisl, theta_start, theta_end);
    
    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double /*theta*/) {
        return compute_secular_RHS_TW(s, p, max_PN_order);
    };
    auto result_tw = integrator.integrate(initial_state, params, rhs_tw, theta_start, theta_end);
    
    std::ofstream outfile = openOutputCsv(output_file);
    if (!outfile.is_open()) {
        return;
    }
    outfile << std::scientific << std::setprecision(10);
    outfile << "theta,p_feireisl,alpha_feireisl,beta_feireisl,p_tw,alpha_tw,beta_tw,"
            << "dp_diff,dalpha_diff,dbeta_diff,max_diff\n";
    
    size_t min_size = std::min(result_feireisl.states.size(), result_tw.states.size());
    for (size_t i = 0; i < min_size; ++i) {
        double theta = result_feireisl.theta[i];
        auto& state_f = result_feireisl.states[i];
        auto& state_tw_i = result_tw.states[i];
        
        double dp_diff = std::abs(state_f.p - state_tw_i.p);
        double da_diff = std::abs(state_f.alpha - state_tw_i.alpha);
        double db_diff = std::abs(state_f.beta - state_tw_i.beta);
        double max_diff = std::max({dp_diff, da_diff, db_diff});
        
        outfile << theta << "," 
                << state_f.p << "," << state_f.alpha << "," << state_f.beta << ","
                << state_tw_i.p << "," << state_tw_i.alpha << "," << state_tw_i.beta << ","
                << dp_diff << "," << da_diff << "," << db_diff << "," << max_diff << "\n";
    }
    outfile.close();
    
    std::cout << "\n=== CONVERGENCE COMPARISON ===" << std::endl;
    std::cout << "Feireisl steps: " << result_feireisl.num_steps << std::endl;
    std::cout << "Tucker-Will steps: " << result_tw.num_steps << std::endl;
    std::cout << "Results written to: " << output_file << std::endl;
}

void compareEvolutionMethodsPhi(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double phi_start,
    double phi_end,
    double tolerance,
    const std::string& output_file
) {
    AdaptiveRK4Integrator integrator(tolerance);

    auto rhs_feireisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
        return compute_feireisl_physical_RHS(s, p, max_PN_order, phi);
    };
    BinaryState initial_tilde_feireisl = approximateTildeFromActual(initial_state, params, false);
    auto result_feireisl_tilde = integrator.integrate(initial_tilde_feireisl, params, rhs_feireisl, phi_start, phi_end);
    auto result_feireisl = reconstructPhysicalTrajectoryFromTilde(result_feireisl_tilde, params, false);

    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
        return compute_tw_physical_RHS(s, p, max_PN_order, phi);
    };
    BinaryState initial_tilde_tw = approximateTildeFromActual(initial_state, params, true);
    auto result_tw_tilde = integrator.integrate(initial_tilde_tw, params, rhs_tw, phi_start, phi_end);
    auto result_tw = reconstructPhysicalTrajectoryFromTilde(result_tw_tilde, params, true);

    std::ofstream outfile = openOutputCsv(output_file);
    if (!outfile.is_open()) {
        return;
    }
    outfile << std::scientific << std::setprecision(10);
    outfile << "phi,p_feireisl_actual,alpha_feireisl_actual,beta_feireisl_actual,"
            << "p_tw_actual,alpha_tw_actual,beta_tw_actual,dp_diff,dalpha_diff,dbeta_diff,max_diff\n";

    size_t min_size = std::min(result_feireisl.states.size(), result_tw.states.size());
    for (size_t i = 0; i < min_size; ++i) {
        double phi = result_feireisl.theta[i];
        auto& state_f = result_feireisl.states[i];
        auto& state_tw_i = result_tw.states[i];

        auto actual_f = state_f;
        auto actual_tw = state_tw_i;

        double dp_diff = std::abs(actual_f.p - actual_tw.p);
        double da_diff = std::abs(actual_f.alpha - actual_tw.alpha);
        double db_diff = std::abs(actual_f.beta - actual_tw.beta);
        double max_diff = std::max({dp_diff, da_diff, db_diff});

        outfile << phi << ","
                << actual_f.p << "," << actual_f.alpha << "," << actual_f.beta << ","
                << actual_tw.p << "," << actual_tw.alpha << "," << actual_tw.beta << ","
                << dp_diff << "," << da_diff << "," << db_diff << "," << max_diff << "\n";
    }
    outfile.close();

    std::cout << "=== PHI-BASED TRANSFORMED COMPARISON ===" << std::endl;
    std::cout << "Feireisl steps: " << result_feireisl.num_steps << std::endl;
    std::cout << "Tucker-Will steps: " << result_tw.num_steps << std::endl;
    std::cout << "Phi-based actual comparison written to: " << output_file << std::endl;
}

struct EpsScanResult {
    double eps;
    double phi;
    double max_diff;
    int feireisl_steps;
    int tw_steps;
};

struct RHSFractionalDiffRow {
    double eps;
    double frac_feireisl;
    double frac_tw;
};

static int mapSecularOrderToQLTOrder(int max_PN_order) {
    if (max_PN_order <= 2) {
        return max_PN_order;
    }
    if (max_PN_order <= 4) {
        return 2;
    }
    return 3;
}

SecularRHS compute_QLT_RHS_phi(const BinaryState& state,
                               const PhysicalParams& params,
                               int max_PN_order,
                               double phi) {
    static const PNCoeffs K = buildCoefficients(true);

    // Evaluate the exact instantaneous QLT equations at the current phi.
    // PN mapping here is: 1->1PN, 2->2PN, 3/4->up to 3.5PN dissipative,
    // 5->up to 4.5PN dissipative.
    int qlt_order = mapSecularOrderToQLTOrder(max_PN_order);
    QLTrhs qlt = computeQLT(
        K,
        state.p,
        state.alpha,
        state.beta,
        phi,
        1.0 / params.eps,
        qlt_order
    );
    return {qlt.dp, qlt.dalpha, qlt.dbeta};
}

static double eccentricity(const BinaryState& state) {
    return std::sqrt(state.alpha * state.alpha + state.beta * state.beta);
}

std::vector<RHSFractionalDiffRow> computeRHSFractionalDifferences(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    const std::vector<double>& epsilon_values
) {
    std::cout << "\n=== Instantaneous RHS Convergence Analysis (de/dphi) ===" << std::endl;
    std::cout << "Initial state: p=" << initial_state.p
              << ", alpha=" << initial_state.alpha
              << ", beta=" << initial_state.beta << std::endl;
    std::cout << "Evaluating RHS at fixed initial state for each epsilon..." << std::endl;
    std::cout << "Observable: |de/dphi_QLT - de/dphi_method| / |de/dphi_QLT|" << std::endl;

    std::vector<RHSFractionalDiffRow> rows;
    rows.reserve(epsilon_values.size());
    const double probe_h = 1e-8;

    for (double eps : epsilon_values) {
        PhysicalParams scan_params = params;
        scan_params.eps = eps;

        auto rhs_qlt = compute_QLT_RHS_phi(initial_state, scan_params, max_PN_order, 0.0);
        BinaryState tilde0_f = approximateTildeFromActual(initial_state, scan_params, false);
        auto rhs_tilde_f = compute_feireisl_physical_RHS(tilde0_f, scan_params, max_PN_order, 0.0);
        BinaryState tilde0_tw = approximateTildeFromActual(initial_state, scan_params, true);
        auto rhs_tilde_tw = compute_tw_physical_RHS(tilde0_tw, scan_params, max_PN_order, 0.0);

        BinaryState qlt_next{
            initial_state.p + probe_h * rhs_qlt[0],
            initial_state.alpha + probe_h * rhs_qlt[1],
            initial_state.beta + probe_h * rhs_qlt[2]
        };
        double e_dot_qlt = (eccentricity(qlt_next) - eccentricity(initial_state)) / probe_h;

        BinaryState tilde1_f{
            tilde0_f.p + probe_h * rhs_tilde_f[0],
            tilde0_f.alpha + probe_h * rhs_tilde_f[1],
            tilde0_f.beta + probe_h * rhs_tilde_f[2]
        };
        BinaryState physical0_f = reconstructPhysicalStateFromTilde(tilde0_f, scan_params, 0.0, false);
        BinaryState physical1_f = reconstructPhysicalStateFromTilde(tilde1_f, scan_params, probe_h, false);
        double e_dot_f = (eccentricity(physical1_f) - eccentricity(physical0_f)) / probe_h;

        BinaryState tilde1_tw{
            tilde0_tw.p + probe_h * rhs_tilde_tw[0],
            tilde0_tw.alpha + probe_h * rhs_tilde_tw[1],
            tilde0_tw.beta + probe_h * rhs_tilde_tw[2]
        };
        BinaryState physical0_tw = reconstructPhysicalStateFromTilde(tilde0_tw, scan_params, 0.0, true);
        BinaryState physical1_tw = reconstructPhysicalStateFromTilde(tilde1_tw, scan_params, probe_h, true);
        double e_dot_tw = (eccentricity(physical1_tw) - eccentricity(physical0_tw)) / probe_h;

        double denom = std::abs(e_dot_qlt) + 1e-30;
        double frac_f = std::abs(e_dot_qlt - e_dot_f) / denom;
        double frac_tw = std::abs(e_dot_qlt - e_dot_tw) / denom;

        rows.push_back({eps, frac_f, frac_tw});

        std::cout << std::scientific << std::setprecision(6)
                  << "  eps=" << std::setprecision(4) << eps
                  << ": Feireisl frac=" << std::setprecision(6) << frac_f
                  << " | TW frac=" << frac_tw << std::endl;
        std::cout << std::defaultfloat;
    }

    return rows;
}

void plotRHSFractionalDifferences(const std::vector<RHSFractionalDiffRow>& rows, const std::string& output_image) {
    if (rows.empty()) {
        std::cout << "No RHS fractional-difference rows to plot." << std::endl;
        return;
    }

    std::string out_path = resolveOutputPath(output_image);
    auto plot_with_python = [&](const std::string& path) -> bool {
        const char* py_cmd =
            "python3 -c \"import sys,math;import matplotlib.pyplot as plt;"
            "rows=[tuple(map(float,l.split())) for l in sys.stdin if l.strip()];"
            "out=sys.argv[1];"
            "xf=[math.log10(r[0]) for r in rows if r[0]>0 and r[1]>0 and math.isfinite(r[1])];"
            "yf=[math.log10(r[1]) for r in rows if r[0]>0 and r[1]>0 and math.isfinite(r[1])];"
            "xt=[math.log10(r[0]) for r in rows if r[0]>0 and r[2]>0 and math.isfinite(r[2])];"
            "yt=[math.log10(r[2]) for r in rows if r[0]>0 and r[2]>0 and math.isfinite(r[2])];"
            "plt.figure(figsize=(8,5));"
            "plt.plot(xf,yf,'o-',label='Feireisl vs QLT') if xf else None;"
            "plt.plot(xt,yt,'s-',label='TW vs QLT') if xt else None;"
            "plt.xlabel('log10(epsilon)'); plt.ylabel('log10(fractional difference)');"
            "plt.title('Instantaneous RHS Fractional Differences');"
            "plt.grid(True, alpha=0.3); plt.legend(); plt.tight_layout(); plt.savefig(out, dpi=200)\"";

        std::string cmd = std::string(py_cmd) + " '" + path + "'";
        FILE* py = popen(cmd.c_str(), "w");
        if (py == nullptr) {
            return false;
        }
        for (const auto& r : rows) {
            std::fprintf(py, "%.17g %.17g %.17g\n", r.eps, r.frac_feireisl, r.frac_tw);
        }
        int py_status = pclose(py);
        return py_status != -1 && WIFEXITED(py_status) && WEXITSTATUS(py_status) == 0;
    };

    if (!plot_with_python(out_path)) {
        std::cout << "matplotlib plotting failed; skipping RHS fractional-difference plot generation." << std::endl;
        return;
    }
    std::cout << "RHS fractional-difference plot written to: " << out_path << std::endl;
}

void compareEvolutionMethodsPhiWithQLT(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double phi_start,
    double phi_end,
    double tolerance,
    const std::string& output_file
) {
    using Clock = std::chrono::steady_clock;
    auto section_start = Clock::now();

    AdaptiveRK4Integrator integrator(tolerance, 1e-13, 1e-3, 10000000);

    std::ofstream outfile = openOutputCsv(output_file);
    if (!outfile.is_open()) {
        return;
    }
    outfile << std::scientific << std::setprecision(10);
    outfile << "phi,p_qlt,alpha_qlt,beta_qlt,p_feireisl,alpha_feireisl,beta_feireisl,"
            << "p_tw,alpha_tw,beta_tw,"
            << "qlt_vs_feireisl_max,qlt_vs_tw_max,feireisl_vs_tw_max\n";
    outfile.flush();

    std::cout << "[progress] Starting QLT integration..." << std::endl;

    auto rhs_qlt = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
        return compute_QLT_RHS_phi(s, p, max_PN_order, phi);
    };
    auto result_qlt = integrator.integrate(initial_state, params, rhs_qlt, phi_start, phi_end);
    auto now = Clock::now();
    std::cout << "[progress] QLT integration complete in "
              << std::chrono::duration_cast<std::chrono::seconds>(now - section_start).count()
              << " s" << std::endl;
    section_start = now;

    std::cout << "[progress] Starting Feireisl integration..." << std::endl;
    auto rhs_feireisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
        return compute_feireisl_physical_RHS(s, p, max_PN_order, phi);
    };
    BinaryState initial_tilde_feireisl = approximateTildeFromActual(initial_state, params, false);
    auto result_feireisl_tilde = integrator.integrate(initial_tilde_feireisl, params, rhs_feireisl, phi_start, phi_end);
    auto result_feireisl = reconstructPhysicalTrajectoryFromTilde(result_feireisl_tilde, params, false);
    now = Clock::now();
    std::cout << "[progress] Feireisl integration complete in "
              << std::chrono::duration_cast<std::chrono::seconds>(now - section_start).count()
              << " s" << std::endl;
    section_start = now;

    std::cout << "[progress] Starting Tucker-Will integration..." << std::endl;
    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
        return compute_tw_physical_RHS(s, p, max_PN_order, phi);
    };
    BinaryState initial_tilde_tw = approximateTildeFromActual(initial_state, params, true);
    auto result_tw_tilde = integrator.integrate(initial_tilde_tw, params, rhs_tw, phi_start, phi_end);
    auto result_tw = reconstructPhysicalTrajectoryFromTilde(result_tw_tilde, params, true);
    now = Clock::now();
    std::cout << "[progress] Tucker-Will integration complete in "
              << std::chrono::duration_cast<std::chrono::seconds>(now - section_start).count()
              << " s" << std::endl;

    std::cout << "[progress] Writing dense comparison samples..." << std::endl;
    section_start = Clock::now();

    const int dense_samples = 20000;
    for (int i = 0; i <= dense_samples; ++i) {
        double frac = static_cast<double>(i) / static_cast<double>(dense_samples);
        double phi = phi_start + frac * (phi_end - phi_start);

        auto qlt = interpolateStateAtX(result_qlt, phi);
        auto fer = interpolateStateAtX(result_feireisl, phi);
        auto tw = interpolateStateAtX(result_tw, phi);

        double qf = std::max({std::abs(qlt.p - fer.p), std::abs(qlt.alpha - fer.alpha), std::abs(qlt.beta - fer.beta)});
        double qt = std::max({std::abs(qlt.p - tw.p), std::abs(qlt.alpha - tw.alpha), std::abs(qlt.beta - tw.beta)});
        double ft = std::max({std::abs(fer.p - tw.p), std::abs(fer.alpha - tw.alpha), std::abs(fer.beta - tw.beta)});

        outfile << phi << ","
                << qlt.p << "," << qlt.alpha << "," << qlt.beta << ","
                << fer.p << "," << fer.alpha << "," << fer.beta << ","
                << tw.p << "," << tw.alpha << "," << tw.beta << ","
                << qf << "," << qt << "," << ft << "\n";
    }
    outfile.close();
    now = Clock::now();
    std::cout << "[progress] Dense sampling complete in "
              << std::chrono::duration_cast<std::chrono::seconds>(now - section_start).count()
              << " s" << std::endl;

    std::cout << "=== PHI-BASED COMPARISON (QLT vs Feireisl vs Tucker-Will) ===" << std::endl;
    std::cout << "QLT steps: " << result_qlt.num_steps << std::endl;
    std::cout << "Feireisl steps: " << result_feireisl.num_steps << std::endl;
    std::cout << "Tucker-Will steps: " << result_tw.num_steps << std::endl;
    std::cout << "Comparison written to: " << output_file << std::endl;
}

EpsScanResult runEpsilonConvergenceScan(
    const BinaryState& initial_state,
    PhysicalParams params,
    int max_PN_order,
    double theta_start,
    double theta_end,
    double tolerance
) {
    AdaptiveRK4Integrator integrator(tolerance);
    auto rhs_feireisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double /*theta*/) {
        return compute_secular_RHS(s, p, max_PN_order);
    };
    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double /*theta*/) {
        return compute_secular_RHS_TW(s, p, max_PN_order);
    };

    auto result_feireisl = integrator.integrate(initial_state, params, rhs_feireisl, theta_start, theta_end);
    auto result_tw = integrator.integrate(initial_state, params, rhs_tw, theta_start, theta_end);

    double max_diff = 0.0;
    size_t min_size = std::min(result_feireisl.states.size(), result_tw.states.size());
    for (size_t i = 0; i < min_size; ++i) {
        auto& state_f = result_feireisl.states[i];
        auto& state_tw_i = result_tw.states[i];
        double dp_diff = std::abs(state_f.p - state_tw_i.p);
        double da_diff = std::abs(state_f.alpha - state_tw_i.alpha);
        double db_diff = std::abs(state_f.beta - state_tw_i.beta);
        double local_max = std::max(dp_diff, std::max(da_diff, db_diff));
        max_diff = std::max(max_diff, local_max);
    }

    return {params.eps, max_diff, max_diff, result_feireisl.num_steps, result_tw.num_steps};
}

struct EpsPhiRow {
    double eps;
    double phi_qlt;
    double phi_feireisl;
    double phi_tw;
    // end-state absolute errors
    double dp_qf, da_qf, db_qf;   // |QLT - Feireisl|
    double dp_qt, da_qt, db_qt;   // |QLT - TW|
    double dp_ft, da_ft, db_ft;   // |Feireisl - TW|
};

std::vector<EpsPhiRow> computeLogPhiVsLogEpsilonData(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    const std::vector<double>& epsilon_values,
    double phi_end,
    double tolerance
) {
    const double scan_min_step = 1e-13;
    const double scan_max_step = 1e-3;
    const int scan_max_steps = 10000000;
    AdaptiveRK4Integrator integrator(tolerance, scan_min_step, scan_max_step, scan_max_steps);

    std::cout << "[eps scan config] tol=" << std::scientific << std::setprecision(1) << tolerance
              << ", min_step=" << scan_min_step
              << ", max_step=" << scan_max_step
              << ", max_steps=" << scan_max_steps << std::endl;
    std::cout << std::defaultfloat;

    std::vector<EpsPhiRow> rows;
    rows.reserve(epsilon_values.size());

    for (double eps : epsilon_values) {
        PhysicalParams scan_params = params;
        scan_params.eps = eps;

        auto rhs_qlt = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
            return compute_QLT_RHS_phi(s, p, max_PN_order, phi);
        };
        auto rhs_feireisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
            return compute_feireisl_physical_RHS(s, p, max_PN_order, phi);
        };
        auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p, double phi) {
            return compute_tw_physical_RHS(s, p, max_PN_order, phi);
        };

        auto result_qlt = integrator.integrate(initial_state, scan_params, rhs_qlt, 0.0, phi_end);
        BinaryState initial_tilde_feireisl = approximateTildeFromActual(initial_state, scan_params, false);
        auto result_feireisl_tilde = integrator.integrate(initial_tilde_feireisl, scan_params, rhs_feireisl, 0.0, phi_end);
        auto result_feireisl = reconstructPhysicalTrajectoryFromTilde(result_feireisl_tilde, scan_params, false);
        BinaryState initial_tilde_tw = approximateTildeFromActual(initial_state, scan_params, true);
        auto result_tw_tilde = integrator.integrate(initial_tilde_tw, scan_params, rhs_tw, 0.0, phi_end);
        auto result_tw = reconstructPhysicalTrajectoryFromTilde(result_tw_tilde, scan_params, true);

        BinaryState s_qlt = result_qlt.states.empty() ? BinaryState{0.0, 0.0, 0.0} : result_qlt.states.back();
        BinaryState s_feireisl = result_feireisl.states.empty() ? BinaryState{0.0, 0.0, 0.0} : result_feireisl.states.back();
        BinaryState s_tw = result_tw.states.empty() ? BinaryState{0.0, 0.0, 0.0} : result_tw.states.back();

        double dp_qf = std::abs(s_qlt.p - s_feireisl.p);
        double da_qf = std::abs(s_qlt.alpha - s_feireisl.alpha);
        double db_qf = std::abs(s_qlt.beta - s_feireisl.beta);
        double dp_qt = std::abs(s_qlt.p - s_tw.p);
        double da_qt = std::abs(s_qlt.alpha - s_tw.alpha);
        double db_qt = std::abs(s_qlt.beta - s_tw.beta);
        double dp_ft = std::abs(s_feireisl.p - s_tw.p);
        double da_ft = std::abs(s_feireisl.alpha - s_tw.alpha);
        double db_ft = std::abs(s_feireisl.beta - s_tw.beta);

        double phi_qlt = result_qlt.theta.empty() ? phi_end : result_qlt.theta.back();
        double phi_feireisl = result_feireisl.theta.empty() ? phi_end : result_feireisl.theta.back();
        double phi_tw = result_tw.theta.empty() ? phi_end : result_tw.theta.back();

        rows.push_back({eps, phi_qlt, phi_feireisl, phi_tw,
                       dp_qf, da_qf, db_qf,
                       dp_qt, da_qt, db_qt,
                       dp_ft, da_ft, db_ft});

        std::cout << std::scientific << std::setprecision(6)
                  << "[eps scan] eps=" << std::setprecision(4) << eps
                  << " | phi_end: qlt=" << phi_qlt
                  << ", feireisl=" << phi_feireisl
                  << ", tw=" << phi_tw << std::endl;
        std::cout << std::setprecision(6)
                  << "           end-state: QLT vs Feireisl: "
                  << "dp=" << dp_qf << ", da=" << da_qf << ", db=" << db_qf << std::endl;
        std::cout << "           end-state: QLT vs TW: "
                  << "dp=" << dp_qt << ", da=" << da_qt << ", db=" << db_qt << std::endl;
        std::cout << "           end-state: Feireisl vs TW: "
                  << "dp=" << dp_ft << ", da=" << da_ft << ", db=" << db_ft << std::endl;
        std::cout << std::defaultfloat;
    }

    return rows;
}

void plotEndStateErrorsVsEpsilon(const std::vector<EpsPhiRow>& rows, const std::string& output_image) {
    if (rows.empty()) {
        std::cout << "No epsilon scan rows to plot error convergence." << std::endl;
        return;
    }

    std::string out_path = resolveOutputPath(output_image);

    // Write data to a temp CSV
    const std::string csv_path = "/tmp/eps_convergence_data.csv";
    {
        std::ofstream csv(csv_path);
        if (!csv.is_open()) {
            std::cout << "Failed to open temp CSV for error convergence plot." << std::endl;
            return;
        }
        csv << "eps,dp_qf,da_qf,db_qf,dp_qt,da_qt,db_qt,dp_ft,da_ft,db_ft\n";
        csv << std::scientific << std::setprecision(15);
        for (const auto& r : rows) {
            csv << r.eps   << ","
                << r.dp_qf << "," << r.da_qf << "," << r.db_qf << ","
                << r.dp_qt << "," << r.da_qt << "," << r.db_qt << ","
                << r.dp_ft << "," << r.da_ft << "," << r.db_ft << "\n";
        }
    }

    // Write Python plotting script to a temp file
    const std::string py_path = "/tmp/plot_eps_convergence.py";
    {
        std::ofstream py(py_path);
        if (!py.is_open()) {
            std::cout << "Failed to write Python script for error convergence plot." << std::endl;
            return;
        }
        py << R"PY(
import sys, csv, math
import matplotlib
import matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

data_file, out_file = sys.argv[1], sys.argv[2]

rows = []
with open(data_file) as f:
    for row in csv.DictReader(f):
        rows.append({k: float(v) for k, v in row.items()})

eps = np.array([r['eps'] for r in rows])
vars_cfg = [
    ('p',     'dp_qf', 'dp_qt'),
    ('alpha', 'da_qf', 'da_qt'),
    ('beta',  'db_qf', 'db_qt'),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (varname, fk, tk) in zip(axes, vars_cfg):
    ef = np.array([r[fk] for r in rows])
    et = np.array([r[tk] for r in rows])

    ax.loglog(eps, ef, 'o-', color='tab:blue',   label='|QLT - Feireisl|', linewidth=1.8, markersize=7)
    ax.loglog(eps, et, 's--', color='tab:orange', label='|QLT - TW|',       linewidth=1.8, markersize=7)

    # Power-law fits
    mask_f = (ef > 0) & np.isfinite(ef)
    mask_t = (et > 0) & np.isfinite(et)
    if mask_f.sum() >= 2:
        slope_f, ic_f = np.polyfit(np.log10(eps[mask_f]), np.log10(ef[mask_f]), 1)
        eps_fit = np.logspace(np.log10(eps.min()), np.log10(eps.max()), 80)
        ax.loglog(eps_fit, 10**(ic_f + slope_f*np.log10(eps_fit)),
                  ':', color='tab:blue', alpha=0.6, label=f'Feireisl fit (slope={slope_f:.2f})')
    if mask_t.sum() >= 2:
        slope_t, ic_t = np.polyfit(np.log10(eps[mask_t]), np.log10(et[mask_t]), 1)
        eps_fit = np.logspace(np.log10(eps.min()), np.log10(eps.max()), 80)
        ax.loglog(eps_fit, 10**(ic_t + slope_t*np.log10(eps_fit)),
                  ':', color='tab:orange', alpha=0.6, label=f'TW fit (slope={slope_t:.2f})')

    ax.set_xlabel(r'$\varepsilon$', fontsize=12)
    ax.set_ylabel(r'$|\Delta ' + varname + r'|$ at $\phi=100$', fontsize=12)
    ax.set_title(f'Convergence to QLT: {varname}', fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

    # Force every data eps value to appear as a tick
    ax.set_xticks(eps)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f'{x:g}')
    )
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())

fig.suptitle('End-state convergence to QLT (log-log): Feireisl vs Tucker-Will', fontsize=13)
plt.tight_layout()
plt.savefig(out_file, dpi=200)
print(f"Saved: {out_file}")
)PY";
    }

    std::string cmd = "python3 " + py_path + " " + csv_path + " " + out_path;
    int ret = std::system(cmd.c_str()); // NOLINT(cert-env33-c)
    if (ret != 0) {
        std::cout << "Warning: Python plot failed for end-state error convergence plot." << std::endl;
        return;
    }
    std::cout << "End-state convergence plot written to: " << out_path << std::endl;
}

void plotLogPhiVsLogEpsilon(const std::vector<EpsPhiRow>& rows, const std::string& output_image) {
    if (rows.empty()) {
        std::cout << "No epsilon scan rows to plot." << std::endl;
        return;
    }

    std::string out_path = resolveOutputPath(output_image);

    auto plot_with_python = [&](const std::string& path) -> bool {
        const char* py_cmd =
            "python3 -c \"import sys,math;import matplotlib.pyplot as plt;"
            "rows=[tuple(map(float,l.split())) for l in sys.stdin if l.strip()];"
            "out=sys.argv[1];"
            "xq=[math.log10(r[0]) for r in rows if r[0]>0.0 and r[1]>0.0 and math.isfinite(r[1])];"
            "yq=[math.log10(r[1]) for r in rows if r[0]>0.0 and r[1]>0.0 and math.isfinite(r[1])];"
            "xf=[math.log10(r[0]) for r in rows if r[0]>0.0 and r[2]>0.0 and math.isfinite(r[2])];"
            "yf=[math.log10(r[2]) for r in rows if r[0]>0.0 and r[2]>0.0 and math.isfinite(r[2])];"
            "xt=[math.log10(r[0]) for r in rows if r[0]>0.0 and r[3]>0.0 and math.isfinite(r[3])];"
            "yt=[math.log10(r[3]) for r in rows if r[0]>0.0 and r[3]>0.0 and math.isfinite(r[3])];"
            "plt.figure(figsize=(8,5));"
            "plt.plot(xq,yq,'o-',label='QLT') if xq else None;"
            "plt.plot(xf,yf,'s-',label='Feireisl') if xf else None;"
            "plt.plot(xt,yt,'^-',label='TW') if xt else None;"
            "plt.xlabel('log10(epsilon)'); plt.ylabel('log10(phi)');"
            "plt.title('log(phi) vs log(epsilon)'); plt.grid(True, alpha=0.3);"
            "plt.legend(); plt.tight_layout(); plt.savefig(out, dpi=200)\"";

        std::string cmd = std::string(py_cmd) + " '" + path + "'";
        FILE* py = popen(cmd.c_str(), "w");
        if (py == nullptr) {
            return false;
        }

        for (const auto& r : rows) {
            std::fprintf(py, "%.17g %.17g %.17g %.17g\n", r.eps, r.phi_qlt, r.phi_feireisl, r.phi_tw);
        }

        int py_status = pclose(py);
        return py_status != -1 && WIFEXITED(py_status) && WEXITSTATUS(py_status) == 0;
    };
    FILE* gp = popen("gnuplot", "w");
    if (gp == nullptr) {
        std::cout << "gnuplot is not installed; skipping log-log plot generation." << std::endl;
        return;
    }

    fprintf(gp, "set terminal pngcairo size 1000,600 enhanced\\n");
    fprintf(gp, "set output '%s'\\n", out_path.c_str());
    fprintf(gp, "set xlabel 'log10(epsilon)'\\n");
    fprintf(gp, "set ylabel 'log10(phi)'\\n");
    fprintf(gp, "set title 'log(phi) vs log(epsilon)'\\n");
    fprintf(gp, "set grid\\n");
    fprintf(gp, "plot '-' with linespoints title 'QLT', '-' with linespoints title 'Feireisl', '-' with linespoints title 'TW'\\n");

    for (const auto& r : rows) {
        if (r.eps > 0.0 && r.phi_qlt > 0.0 && std::isfinite(r.phi_qlt)) {
            fprintf(gp, "%.16e %.16e\\n", std::log10(r.eps), std::log10(r.phi_qlt));
        }
    }
    fprintf(gp, "e\\n");

    for (const auto& r : rows) {
        if (r.eps > 0.0 && r.phi_feireisl > 0.0 && std::isfinite(r.phi_feireisl)) {
            fprintf(gp, "%.16e %.16e\\n", std::log10(r.eps), std::log10(r.phi_feireisl));
        }
    }
    fprintf(gp, "e\\n");

    for (const auto& r : rows) {
        if (r.eps > 0.0 && r.phi_tw > 0.0 && std::isfinite(r.phi_tw)) {
            fprintf(gp, "%.16e %.16e\\n", std::log10(r.eps), std::log10(r.phi_tw));
        }
    }
    fprintf(gp, "e\\n");

    int close_status = pclose(gp);
    if (close_status == -1 || !WIFEXITED(close_status) || WEXITSTATUS(close_status) != 0) {
        std::cout << "gnuplot failed; trying matplotlib fallback..." << std::endl;
        if (!plot_with_python(out_path)) {
            std::cout << "matplotlib fallback also failed; skipping log-log plot generation." << std::endl;
            return;
        }
    }
    std::cout << "log(phi) vs log(epsilon) plot written to: " << out_path << std::endl;
}

void writeEpsilonConvergenceCSV(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double theta_start,
    double theta_end,
    double tolerance,
    const std::string& output_file
) {
    std::ofstream outfile = openOutputCsv(output_file);
    if (!outfile.is_open()) {
        return;
    }
    outfile << std::scientific << std::setprecision(10);
    outfile << "eps,phi,max_diff,feireisl_steps,tw_steps\n";

    std::vector<double> eps_values = {0.032, 0.016, 0.008, 0.004, 0.002};
    for (double eps : eps_values) {
        PhysicalParams scan_params = params;
        scan_params.eps = eps;
        auto result = runEpsilonConvergenceScan(initial_state, scan_params, max_PN_order,
                                                theta_start, theta_end, tolerance);
        outfile << result.eps << "," << result.phi << "," << result.max_diff << ","
                << result.feireisl_steps << "," << result.tw_steps << "\n";
    }
    outfile.close();
    std::cout << "Epsilon convergence scan written to: " << output_file << std::endl;
}

// ============================================================================
// DELTA-PHI VS EPSILON SCAN  (stop integration when p reaches p_final)
// ============================================================================

struct DeltaPhiRow {
    double eps;
    double phi_qlt;
    double phi_feireisl;
    double phi_tw;
    double delta_phi_feireisl; // phi_qlt - phi_feireisl
    double delta_phi_tw;       // phi_qlt - phi_tw
};

std::vector<DeltaPhiRow> computeDeltaPhiVsEpsilonScan(
    const PhysicalParams&      base_params,
    int                        max_PN_order,
    const BinaryState&         initial_physical,   // p_init=50, alpha=0.1, beta=0.1
    double                     p_final,            // 20
    const std::vector<double>& epsilon_values,
    double                     tolerance,
    bool                       fast_debug = false
) {
    // Keep QLT as the reference, but give it a dedicated budget because it must
    // resolve fast orbital oscillations. Secular methods remain cheap and can
    // use large macro-steps.
    const double qlt_min_step = 1e-6;
    const double qlt_max_step = 1.0;
    const int integration_max_steps = fast_debug ? 10000000 : 100000000;
    
    const int    qlt_retry_factor = 4;
    const int    qlt_retry_count = 3;

    const double secular_min_step = 1e-6;
    const double secular_max_step = 1.0;

    std::vector<DeltaPhiRow> rows;
    rows.reserve(epsilon_values.size());

    std::cout << "\n=== Delta-phi vs epsilon scan (p: " << initial_physical.p
              << " -> " << p_final << ") ==="  << std::endl;
    std::cout << "[scan config] tol=" << std::scientific << std::setprecision(1)
              << tolerance
              << ", qlt_max_step=" << qlt_max_step
              << ", secular_max_step=" << secular_max_step
              << ", max_steps=" << integration_max_steps
              << std::defaultfloat << std::endl;

    for (size_t idx = 0; idx < epsilon_values.size(); ++idx) {
        double eps = epsilon_values[idx];
        PhysicalParams p = base_params;
        p.eps = eps;

        std::cout << "[progress] eps index " << (idx + 1) << "/" << epsilon_values.size()
                  << " (eps=" << std::scientific << std::setprecision(4) << eps << ")"
                  << std::defaultfloat << std::endl;

        AdaptiveRK4Integrator qlt_integrator(tolerance, qlt_min_step, qlt_max_step, integration_max_steps);
        AdaptiveRK4Integrator secular_integrator(tolerance, secular_min_step, secular_max_step, integration_max_steps);

        // ── QLT: physical state, stop when state.p <= p_final ──────────────
        auto rhs_qlt = [max_PN_order](const BinaryState& s, const PhysicalParams& pp, double phi) {
            return compute_QLT_RHS_phi(s, pp, max_PN_order, phi);
        };
        auto get_p_physical = [](const BinaryState& s, const PhysicalParams&, double) {
            return s.p;
        };

        // ── Feireisl: tilde state, stop when reconstructed p <= p_final ────
        BinaryState tilde0_f = approximateTildeFromActual(initial_physical, p, false);
        BinaryState final_actual{p_final, initial_physical.alpha, initial_physical.beta};
        BinaryState tilde_final_f = approximateTildeFromActual(final_actual, p, false);
        std::cout << std::scientific << std::setprecision(16)
                  << "    [init tilde JF] p_tilde=" << tilde0_f.p
                  << ", alpha_tilde=" << tilde0_f.alpha
                  << ", beta_tilde=" << tilde0_f.beta
                  << std::defaultfloat << std::endl;
        std::cout << std::scientific << std::setprecision(16)
                  << "    [target tilde JF] p_tilde_final=" << tilde_final_f.p
                  << " (from p_final=" << p_final << ")"
                  << std::defaultfloat << std::endl;
        auto rhs_feireisl = [max_PN_order](const BinaryState& s, const PhysicalParams& pp, double phi) {
            return compute_feireisl_physical_RHS(s, pp, max_PN_order, phi);
        };
        auto get_p_feireisl = [](const BinaryState& s, const PhysicalParams&, double) {
            return s.p;
        };
        auto run_until_target_secular = [&](const std::string& label,
                                            const BinaryState& init_state,
                                            double p_stop,
                                            const auto& rhs,
                                            const auto& get_p) {
            bool reached = false;
            double phi_hit = secular_integrator.integrateUntilP(
                init_state, p, rhs, 0.0, p_stop, get_p, &reached);

            if (!reached) {
                std::cout << "    [warn] " << label
                          << " did not reach stop target within max_steps (phi="
                          << std::scientific << std::setprecision(16) << phi_hit
                          << ")." << std::defaultfloat << std::endl;
            }

            return std::make_pair(phi_hit, reached);
        };

        auto run_until_target_qlt = [&](const std::string& label,
                                        const BinaryState& init_state,
                                        const auto& rhs,
                                        const auto& get_p) {
            bool reached = false;
            int max_steps_cur = integration_max_steps;
            double phi_hit = 0.0;

            for (int attempt = 0; attempt <= qlt_retry_count; ++attempt) {
                AdaptiveRK4Integrator retry_integrator(tolerance, qlt_min_step, qlt_max_step, max_steps_cur);
                phi_hit = retry_integrator.integrateUntilP(
                    init_state, p, rhs, 0.0, p_final, get_p, &reached);
                if (reached) {
                    break;
                }
                if (attempt < qlt_retry_count) {
                    std::cout << "    [extend] " << label
                              << " did not reach p_final; increasing max_steps to "
                              << (max_steps_cur * qlt_retry_factor) << std::endl;
                    max_steps_cur *= qlt_retry_factor;
                }
            }

            if (!reached) {
                std::cout << "    [warn] " << label
                          << " did not reach p_final after retries (phi="
                          << std::scientific << std::setprecision(16) << phi_hit
                          << ")." << std::defaultfloat << std::endl;
            }

            return std::make_pair(phi_hit, reached);
        };

        auto qlt_run = run_until_target_qlt("QLT", initial_physical, rhs_qlt, get_p_physical);
        double phi_qlt = qlt_run.first;
        bool reached_qlt = qlt_run.second;

        auto feireisl_run = run_until_target_secular("Feireisl", tilde0_f, tilde_final_f.p, rhs_feireisl, get_p_feireisl);
        double phi_feireisl = feireisl_run.first;
        bool reached_feireisl = feireisl_run.second;

        // ── Tucker-Will: same pattern ────────────────────────────────────────
        BinaryState tilde0_tw = approximateTildeFromActual(initial_physical, p, true);
        BinaryState tilde_final_tw = approximateTildeFromActual(final_actual, p, true);
        std::cout << std::scientific << std::setprecision(6)
                  << "    [init tilde TW] p_tilde=" << tilde0_tw.p
                  << ", alpha_tilde=" << tilde0_tw.alpha
                  << ", beta_tilde=" << tilde0_tw.beta
                  << std::defaultfloat << std::endl;
        std::cout << std::scientific << std::setprecision(16)
                  << "    [target tilde TW] p_tilde_final=" << tilde_final_tw.p
                  << " (from p_final=" << p_final << ")"
                  << std::defaultfloat << std::endl;
        auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& pp, double phi) {
            return compute_tw_physical_RHS(s, pp, max_PN_order, phi);
        };

        // Diagnostic: compare effective initial physical dp/dphi across methods.
        {
            const double h_probe = 1e-6;
            auto qlt_rhs0 = rhs_qlt(initial_physical, p, 0.0);
            double dp0_qlt = qlt_rhs0[0];

            auto f_rhs0 = rhs_feireisl(tilde0_f, p, 0.0);
            double dp0_f_raw = f_rhs0[0];
            BinaryState f_step{
                tilde0_f.p + h_probe * f_rhs0[0],
                tilde0_f.alpha + h_probe * f_rhs0[1],
                tilde0_f.beta + h_probe * f_rhs0[2]
            };
            double p0_f = reconstructPhysicalStateFromTilde(tilde0_f, p, 0.0, false).p;
            double p1_f = reconstructPhysicalStateFromTilde(f_step, p, h_probe, false).p;
            double dp0_f = (p1_f - p0_f) / h_probe;

            auto tw_rhs0 = rhs_tw(tilde0_tw, p, 0.0);
            double dp0_tw_raw = tw_rhs0[0];
            BinaryState tw_step{
                tilde0_tw.p + h_probe * tw_rhs0[0],
                tilde0_tw.alpha + h_probe * tw_rhs0[1],
                tilde0_tw.beta + h_probe * tw_rhs0[2]
            };
            double p0_tw = reconstructPhysicalStateFromTilde(tilde0_tw, p, 0.0, true).p;
            double p1_tw = reconstructPhysicalStateFromTilde(tw_step, p, h_probe, true).p;
            double dp0_tw = (p1_tw - p0_tw) / h_probe;

            std::cout << std::scientific << std::setprecision(16)
                      << "    [diag raw dp/dphi @start] QLT=" << dp0_qlt
                      << ", JF=" << dp0_f_raw
                      << ", TW=" << dp0_tw_raw
                      << std::defaultfloat << std::endl;
            std::cout << std::scientific << std::setprecision(16)
                      << "    [diag eff dp/dphi @start] QLT=" << dp0_qlt
                      << ", JF=" << dp0_f
                      << ", TW=" << dp0_tw
                      << std::defaultfloat << std::endl;
        }

        auto get_p_tw = [](const BinaryState& s, const PhysicalParams&, double) {
            return s.p;
        };
        auto tw_run = run_until_target_secular("TW", tilde0_tw, tilde_final_tw.p, rhs_tw, get_p_tw);
        double phi_tw = tw_run.first;
        bool reached_tw = tw_run.second;

        double dphi_f  = std::numeric_limits<double>::quiet_NaN();
        double dphi_tw = std::numeric_limits<double>::quiet_NaN();
        if (reached_qlt && reached_feireisl) {
            dphi_f = phi_qlt - phi_feireisl;
        }
        if (reached_qlt && reached_tw) {
            dphi_tw = phi_qlt - phi_tw;
        }

        rows.push_back({eps, phi_qlt, phi_feireisl, phi_tw, dphi_f, dphi_tw});

        std::cout << std::scientific << std::setprecision(16)
                  << "  eps=" << eps
                  << " | phi_qlt=" << phi_qlt
                  << ", phi_F=" << phi_feireisl
                  << ", phi_TW=" << phi_tw
                  << " | delta_phi_F=" << dphi_f
                  << ", delta_phi_TW=" << dphi_tw
                  << " | reached(QLT/F/TW)=" << (reached_qlt ? "Y" : "N")
                  << "/" << (reached_feireisl ? "Y" : "N")
                  << "/" << (reached_tw ? "Y" : "N")
                  << std::defaultfloat << std::endl;
    }

    return rows;
}

void plotDeltaPhiVsEpsilon(const std::vector<DeltaPhiRow>& rows, const std::string& output_image) {
    if (rows.empty()) {
        std::cout << "No delta-phi rows to plot." << std::endl;
        return;
    }

    std::string out_path = resolveOutputPath(output_image);
    const std::string csv_path = "/tmp/delta_phi_data.csv";
    const std::string py_path  = "/tmp/plot_delta_phi.py";

    {
        std::ofstream csv(csv_path);
        if (!csv.is_open()) {
            std::cout << "Failed to open temp CSV for delta-phi plot." << std::endl;
            return;
        }
        csv << "eps,phi_qlt,phi_feireisl,phi_tw,delta_phi_F,delta_phi_TW\n";
        csv << std::scientific << std::setprecision(16);
        for (const auto& r : rows) {
            csv << r.eps            << ","
                << r.phi_qlt        << ","
                << r.phi_feireisl   << ","
                << r.phi_tw         << ","
                << r.delta_phi_feireisl << ","
                << r.delta_phi_tw   << "\n";
        }
    }

    {
        std::ofstream py(py_path);
        if (!py.is_open()) {
            std::cout << "Failed to write Python script for delta-phi plot." << std::endl;
            return;
        }
        py << R"PY(
import sys, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

data_file, out_file = sys.argv[1], sys.argv[2]

rows = []
with open(data_file) as f:
    for row in csv.DictReader(f):
        rows.append({k: float(v) for k, v in row.items()})

eps  = np.array([r['eps']           for r in rows])
df   = np.array([abs(r['delta_phi_F'])   for r in rows])
dt   = np.array([abs(r['delta_phi_TW'])  for r in rows])

fig, ax = plt.subplots(figsize=(8, 5))

mask_f = (df > 0) & np.isfinite(df)
mask_t = (dt > 0) & np.isfinite(dt)

xf = np.log10(eps[mask_f])
yf = np.log10(df[mask_f])
xt = np.log10(eps[mask_t])
yt = np.log10(dt[mask_t])

ax.plot(xf, yf, 'o-', color='tab:blue', linewidth=2.2, markersize=7,
    label=r'$\log_{10}|\phi_{\rm QLT}-\phi_{\rm Feireisl}|$')
ax.plot(xt, yt, 's--', color='tab:orange', linewidth=2.0, markersize=7,
    label=r'$\log_{10}|\phi_{\rm QLT}-\phi_{\rm TW}|$')

if mask_f.sum() >= 2:
    slope_f, ic_f = np.polyfit(np.log10(eps[mask_f]), np.log10(df[mask_f]), 1)
    x_line = np.linspace(np.log10(eps.min()), np.log10(eps.max()), 80)
    y_line = ic_f + slope_f * x_line
    ax.plot(x_line, y_line, ':', color='tab:blue', alpha=0.7,
        label=f'Feireisl fit slope={slope_f:.2f}')

if mask_t.sum() >= 2:
    slope_t, ic_t = np.polyfit(np.log10(eps[mask_t]), np.log10(dt[mask_t]), 1)
    x_line = np.linspace(np.log10(eps.min()), np.log10(eps.max()), 80)
    y_line = ic_t + slope_t * x_line
    ax.plot(x_line, y_line, ':', color='tab:orange', alpha=0.7,
        label=f'TW fit slope={slope_t:.2f}')

ax.set_xlabel(r'$\log_{10}(\varepsilon)$', fontsize=13)
ax.set_ylabel(r'$\log_{10}(|\Delta\phi|)$', fontsize=13)
ax.set_title(r'QLT-Referenced Phase Errors: $\log_{10}|\Delta\phi|$ vs $\log_{10}(\varepsilon)$', fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(out_file, dpi=200)
print(f"Saved: {out_file}")
)PY";
    }

    std::string cmd = "python3 " + py_path + " " + csv_path + " " + out_path;
    int ret = std::system(cmd.c_str()); // NOLINT(cert-env33-c)
    if (ret != 0) {
        std::cout << "Warning: Python plot failed for delta-phi vs epsilon." << std::endl;
        return;
    }
    std::cout << "Delta-phi vs epsilon plot written to: " << out_path << std::endl;
}

// ============================================================================
// MAIN
// ============================================================================

int main(int argc, char** argv) {
    (void)argc;
    initializeOutputBaseDir((argv != nullptr) ? argv[0] : nullptr);

    PhysicalParams params{1.0, 1.0, 0.25, 0.0, 1.0};
    const bool fast_debug = false;
    double tolerance = 1e-6;
    int max_PN_order = 5;   // 4.5PN order

    std::cout << "=== Adaptive RK4 Integration Comparison ===" << std::endl;
    std::cout << "PN order: " << max_PN_order << " (4.5PN)" << std::endl;
    std::cout << "Tolerance: " << tolerance << std::endl;
    std::cout << "Mode: " << (fast_debug ? "FAST DEBUG" : "FULL ACCURACY") << std::endl;

    // ── Delta-phi scan: p_init=50 -> p_final=20 ─────────────────────────────
    BinaryState initial_p50{50.0, 0.1, 0.1};
    double p_final = 20.0;
    std::vector<double> delta_phi_eps_values = fast_debug
        ? std::vector<double>{1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625}
        : std::vector<double>{1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625};

    auto delta_phi_rows = computeDeltaPhiVsEpsilonScan(
        params,
        max_PN_order,
        initial_p50,
        p_final,
        delta_phi_eps_values,
        tolerance,
        fast_debug
    );
    plotDeltaPhiVsEpsilon(delta_phi_rows, "delta_phi_vs_epsilon_p50top20_fromcpp_rk4.png");

    return 0;
}
