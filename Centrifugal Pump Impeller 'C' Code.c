#include <math.h>
#include <stdio.h>

#if !defined(M_E)
static long double M_E = 2.718281828459045235360287471352662498; // e

#endif // M_E

#if !defined(M_PI)
static long double M_PI = 3.141592653589793238462643383279502884; // pi

#endif // M_PI

double pump_head(double p1, double p2, double rho, double g) {
    return (p2 * 1000 - p1 * 1000) / (rho * g);
}

double specific_speed(int n, double Q, double H) {
    return 3.65 * n * (sqrt(Q) / pow(H, 3.0 / 4));
}

double normalized_diameter(double Q, int n) {
    return 4.25 * pow(Q / n, 1.0 / 3);
}

double volumetric_efficiency(double ns) {
    double a = 0.68;
    return 1.0 / (1 + a * pow(ns, -0.66));
}

double hydraulic_efficiency(double d1_adj) {
    return (1 - 0.42 / pow((log(d1_adj) - 0.172), 2));
}

double pump_efficiency(double nv, double nh, double nm) { return nv * nh * nm; }

double shaft_power(double rho, double Q, double H, double g, double pump_eff) {
    return (rho * Q * H * g) / (1000 * pump_eff);
}

double torque(double N, double n) { return 9600 * N / n; }

double pump_shaft_diameter(double M, double t_al) {
    return pow(M / (0.2 * t_al), 1 / 3.0);
}

double impeller_hub_diameter(double ds) { return 1.4 * ds * 10; }

double hub_length(double d_hub) { return 1.4 * d_hub; }

double peripheral_vel_impeller_inlet(int D1, int n) {
    return M_PI * D1 * n / 60;
}

double vel_impeller_eye(double Q, double nv, int D0, double d_hub) {
    return (4 * Q) / (nv * M_PI * (pow(D0, 2) - pow(d_hub, 2)));
}

double inlet_vane_width(double Q, double c1, double D1, double mu) {
    return Q / (c1 * M_PI * (D1 / 1000) * mu);
}

double peripheral_vel(double c_2r, double beta_2, double g, double H,
                      double nh) {
    double x = c_2r / (2 * tan(beta_2));
    return x + sqrt(x + (g * H / nh));
}

double no_impeller_vanes(double m, double beta_1v, double beta_2v) {
    return 1 + 6.5 * ((m + 1) / (m - 1)) *
                   sin((beta_1v + beta_2v) / 2 * M_PI /
                       180); // 1 added for proper estimates
}

int main(void) {
    system("color 70");
    double rho = 1000;
    double g = 9.81;
    double Q = 0;  // m^3 / s
    double p1 = 0; // kPA
    double p2 = 0; // kPA

    printf("\nVolumetric Flow Rate Q(m^3 / s) = ");
    scanf("%lf", &Q);

    printf("\nInlet pressure p1(kPa) = ");
    scanf("%lf", &p1);

    printf("\nOutlet pressure p2(kPa) = ");
    scanf("%lf", &p2);
    int n = 0; // rpm

    printf("\nRevolution per minute n (rpm) = ");
    scanf("%d", &n);

    double H = pump_head(p1, p2, rho, g); // m
    double ns = specific_speed(n, Q, H);
    double d1_adj = normalized_diameter(Q, n); // m

    double nv = volumetric_efficiency(ns);
    double nh = hydraulic_efficiency(d1_adj);
    double nm = 0.93; // Mechanical Efficiency

    double pump_eff = pump_efficiency(nv, nh, nm);
    double N = shaft_power(rho, Q, H, g, pump_eff); // kW
    double M = torque(N, n) * 100;                  // N.cm
    double t_al = 1500;                             // N / cm^2
    double ds = pump_shaft_diameter(M, t_al);       // cm
    double d_hub = impeller_hub_diameter(ds);       // mm

    double D0 = 140;     // mm
    double D1 = D0 + 20; // mm (Impeller Inlet Diameter)

    double l_hub = hub_length(d_hub); // mm

    double u1 = peripheral_vel_impeller_inlet(D1, n) / 1000; // m/s

    // Velocity at impeller eye with D_O = D_1
    D0 = D1;
    double c0 = vel_impeller_eye(Q, nv, D0, d_hub) * pow(10, 6); // m/s

    double c1 = c0;
    double beta_1 = atan(c1 / u1) * 180 / M_PI; // degree

    // Assuming i = 4, vane angle at inlet is
    int I = 4;
    double beta_1v = beta_1 + I;

    double mu = 0.9; // Blockage factor

    double b1 = inlet_vane_width(Q, c1, D1, mu) / 2 * 1000; // mm

    double beta_2 = 32; // degree
    double c_2r = 5.6;

    double u2 = peripheral_vel(c_2r, beta_2 * M_PI / 180, g, H, nh); // m/s
    double D2 = 60 * u2 / (M_PI * n);                                // m
    double m = D2 / (D1 / 1000); // outlet-to-inlet diameter ratio

    // Assuming c_1r = c_2r,
    double b2 = b1 * (D1 / 1000) / D2; // mm

    int z = no_impeller_vanes(m, beta_1v * M_PI / 180, beta_2);

    printf("\n\n");

    printf("\t******Impeller Dimensions******\n\n");
    printf("\tPump Head = %0.2lf m\n", H);
    printf("\tSpecific speed = %0.2lf\n", ns);
    printf("\tNormalized Diameter = %0.2lf m\n\n", d1_adj);
    printf("\tVolumetric Efficiency = %0.2lf\n", nv);
    printf("\tHydraulic Efficiency = %0.2lf\n", nh);
    printf("\tMechanical Efficiency = %0.2lf\n", nm);
    printf("\tPump Efficiency = %0.2lf\n\n", pump_eff);
    printf("\tShaft Power = %0.2lf kW\n", N);
    printf("\tTorque = %0.2lf N.cm\n", M);
    printf("\tAllowable torsional stress = %0.2lf cm\n", ds);
    printf("\tImpeller Hub Diameter = %0.2lf mm\n", d_hub);
    printf("\tImpeller Inlet Diameter = %0.2lf mm\n\n", D1);
    printf("\tPeripheral velocity at impeller inlet = %0.2lf m/s\n", u1);
    printf("\tVelocity at impeller eye = %0.2lf m/s\n", c0);
    printf("\tVane angle at inlet = %0.2lf degree\n", beta_1v);
    printf("\tVane angle at outlet = 32 degree\n\n");
    printf("\tInlet vane width = %0.2lf mm\n", b1);
    printf("\tPeripheral velocity = %0.2lf m/s\n", u2);
    printf("\tOutlet-to-inlet diameter ratio = %0.2lf\n", m);
    printf("\tDischarge vane width = %0.2lf mm\n", b2);
    printf("\tNumber of impeller vanes = %d\n", z);

    return 0;
}
