#include <dde_solver.hpp>
#include <dde_options.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

struct dde_solver_ex441_test : public testing::Test {
  static double tau, omega;
  int sol_npts, sol_flag, sol_ne;
  double* sol_t_ptr;
  double* sol_y_ptr;
  double* sol_te_ptr;
  double* sol_ye_ptr;
  std::vector<int> nvar{4, 2};
  int n_nvar = nvar.size();
  std::vector<double> tspan{0.0, 350.0};
  int n_tspan = tspan.size();
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{default_rel_err()};
  std::vector<double> ae_vector{default_abs_err()};
  int nre = re_vector.size();
  int nae = ae_vector.size();

  /* we have to pass arrays to fortran, but n_xxx = 0
     indicates these arrays are dummies */
  int n_jumps = 0;
  int n_thit = 0;
  int n_direction = 0;
  int n_isterminal = 0;
  double jumps[1];
  double thit_exactly[1];
  int direction[1];
  bool isterminal[1];

  std::vector<double> delay{ 42.0, 0.15 };
  std::vector<double> history{15.0, 0.0, 2.0, 3.0};

  dde_solver_ex441_test() {}

  static void ddes_cc(double* t, int* n, int* nlags, double y[], double z[], double dy[]) {
    // physical parameters
    const double a = 0.330, d = 0.006, lambda = 0.308,
      gamma = 0.040 , epsilon=0.060;

    // local solution variables and delayed solution values:
    double s = y[0];
    double e = y[1];
    double i = y[2];
    double r = y[3];
    double itau   = z[2];
    double somega = z[4];
    double eomega = z[5];
    double iomega = z[6];
    double romega = z[7];

    double noft = s + e + i + r;
    double nomega = somega + eomega + iomega + romega;

    // local derivatives:
    double dsdt = a - d*s - lambda*((s*i)/noft) + gamma*itau*exp(-d*tau);
    double dedt = lambda*((s*i)/noft) -
      lambda*((somega*iomega)/nomega)* exp(-d*omega) - d*e;
    double didt = lambda*((somega*iomega)/nomega)*exp(-d*omega) -
      (gamma+epsilon+d)*i;
    double drdt = gamma*i - gamma*itau*exp(-d*tau) - d*r;

    // derivatives for the integrator:
    dy[0] = dsdt;
    dy[1] = dedt;
    dy[2] = didt;
    dy[3] = drdt;
  }
};

double dde_solver_ex441_test::tau = 42.0;
double dde_solver_ex441_test::omega = 0.15;

TEST_F(dde_solver_ex441_test, sol) {

  integrate_dde_4(&n_nvar, nvar.data(),
                  ddes_cc, &nvar[1], delay.data(), &nvar[0], history.data(),
                  &n_tspan, tspan.data(),
                  &sol_npts, &sol_flag, &sol_ne,
                  &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr,
                  &nre, re_vector.data(),
                  &nae, ae_vector.data(),
                  &n_jumps, jumps,
                  &n_thit, thit_exactly,
                  &n_direction, direction,
                  &n_isterminal, isterminal,
                  &opts_cc,
                  NULL, NULL, NULL, NULL);

  EXPECT_EQ(sol_flag, 0);       // solution successful

  // Fortran solutions at the last point
  std::cout.precision(10);
  EXPECT_FLOAT_EQ(sol_t_ptr[sol_npts - 1], tspan[1]);

  EXPECT_FLOAT_EQ(sol_y_ptr[sol_npts - 1]   , 5.231271201   );
  EXPECT_FLOAT_EQ(sol_y_ptr[2*sol_npts - 1] , 0.05490837936 );
  EXPECT_FLOAT_EQ(sol_y_ptr[3*sol_npts - 1] , 3.985113815   );
  EXPECT_FLOAT_EQ(sol_y_ptr[4*sol_npts - 1] , 5.915635409   );
}
