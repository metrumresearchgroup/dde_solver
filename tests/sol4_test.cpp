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

  // output
  int sol_npts, sol_flag, sol_ne;
  double* sol_t_ptr;
  double* sol_y_ptr;
  double* sol_te_ptr;
  double* sol_ye_ptr;
  double* sol_queue_ptr;
  double* sol_yoft_ptr;
  double* sol_tqueue_ptr;
  int* sol_stats_ptr;
  int* sol_ie_ptr;
  int* sol_ipoint_ptr;
  bool sol_shift;
  double sol_tshift;

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

  static void ddes_cc(const double* t, const int* n, const int* nlags, const double y[], const double z[], double dy[]) {
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
                  &sol_npts, &sol_flag, &sol_ne,                    // output
                  &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr, // output
                  &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,   // output
                  &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,     // output
                  &sol_shift, &sol_tshift,                          // output
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

struct dde_solver_ex442_test : public testing::Test {
  // output
  int sol_npts, sol_flag, sol_ne;
  double* sol_t_ptr;
  double* sol_y_ptr;
  double* sol_te_ptr;
  double* sol_ye_ptr;
  double* sol_queue_ptr;
  double* sol_yoft_ptr;
  double* sol_tqueue_ptr;
  int* sol_stats_ptr;
  int* sol_ie_ptr;
  int* sol_ipoint_ptr;
  bool sol_shift;
  double sol_tshift;

  std::vector<int> nvar{1, 1};
  int n_nvar = nvar.size();
  std::vector<double> tspan{0.0, 100.0};
  int n_tspan = tspan.size();
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{1.0E-6};
  std::vector<double> ae_vector{1.0E-6};
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

  std::vector<double> delay{ 2.0 };
  std::vector<double> history{0.5};

  dde_solver_ex442_test() {
    opts_cc.interpolation = true;
  }

  static void ddes_cc(const double* t, const int* n, const int* nlags, const double y[], const double z[], double dy[]) {
    dy[0] = 2.0 * z[0] / (1.0 + pow(z[0], 9.65)) - y[0];
  }
};

TEST_F(dde_solver_ex442_test, sol) {

  integrate_dde_4(&n_nvar, nvar.data(),
                  ddes_cc, &nvar[1], delay.data(), &nvar[0], history.data(),
                  &n_tspan, tspan.data(),
                  &sol_npts, &sol_flag, &sol_ne,                    // output
                  &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr, // output
                  &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,   // output
                  &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,     // output
                  &sol_shift, &sol_tshift,                          // output
                  &nre, re_vector.data(),
                  &nae, ae_vector.data(),
                  &n_jumps, jumps,
                  &n_thit, thit_exactly,
                  &n_direction, direction,
                  &n_isterminal, isterminal,
                  &opts_cc,
                  NULL, NULL, NULL, NULL);

  EXPECT_EQ(sol_npts, 547);
  EXPECT_EQ(sol_flag, 0);       // solution successful

  // TODO: test interpolated results
}

struct dde_solver_ex443_test : public testing::Test {
  static bool useODEmodel;
  static constexpr double d=5.0, c=0.5, g=1.0, n=100.0;

  // output
  int sol_npts, sol_flag, sol_ne;
  double* sol_t_ptr;
  double* sol_y_ptr;
  double* sol_te_ptr;
  double* sol_ye_ptr;
  double* sol_queue_ptr;
  double* sol_yoft_ptr;
  double* sol_tqueue_ptr;
  int* sol_stats_ptr;
  int* sol_ie_ptr;
  int* sol_ipoint_ptr;
  bool sol_shift;
  double sol_tshift;

  std::vector<int> nvar{6, 1, 1};
  int n_nvar = nvar.size();
  std::vector<double> tspan{0.0, 4 * d};
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

  std::vector<double> delay{ d };
  std::vector<double> history{ 0.8*n, 0.2*n, 0.0, 0.0, 0.0, 0.0 };

  dde_solver_ex443_test() {}

  static void ddes_cc(const double* t, const int* nn, const int* nlags, const double y[], const double z[], double dy[]) {
    if (useODEmodel) {
      // integrate the ode model for t <= d:
      dy[0] = -y[0] * y[2]  +  g * y[1];
      dy[1] = -dy[0];
      dy[3] = exp(g * (*t)) * y[1];
      dy[4] = (*t) * exp(g * (*t)) * y[0] * y[2];
      dy[5] = exp(g * (*t)) * y[0] * y[2];
      dy[2] = (c/n) * exp(-g * (*t)) * ((dy[3] + dy[4])-g * (y[3] + y[4]));
    } else {
      // integrate the dde model for t >= d:
      dy[0] = -y[0] * y[2]  +  g * y[1];
      dy[1] = -dy[0];
      dy[3] = exp(g * (*t)) * y[1] - exp(g * ((*t)-d)) * z[1];
      dy[4] = d * exp(g * (*t)) * y[0] * y[2] - y[5];
      dy[5] = exp(g * (*t)) * y[0] * y[2] - exp(g * ((*t)-d)) * z[0] * z[2];
      dy[2] = (c/n) * exp(-g * (*t)) * ((dy[3] + dy[4])-g * (y[3] + y[4]));
    }
  }

    // event function
  static void ef_cc(const double* t, const int* n, const int* nlag,
                    const double y[], const double dy[], const double z[], double g[]) {
    double ylag = z[0];
    g[0] = (*t) - d;
  }

  static void chng_cc(const int* nevent, const double* tevent,
                      double yevent[], double dyevent[], double* hinit,
                      const int* n_direction, int direction[],
                      const int* n_isterminal, bool isterminal[], bool* quit) {
    if (*nevent == 1) {
      useODEmodel = false;
    }
  }
};

bool dde_solver_ex443_test::useODEmodel = false;

TEST_F(dde_solver_ex443_test, sol) {

  useODEmodel = true;

  integrate_dde_4(&n_nvar, nvar.data(),
                  ddes_cc, &nvar[1], delay.data(), &nvar[0], history.data(),
                  &n_tspan, tspan.data(),
                  &sol_npts, &sol_flag, &sol_ne,                    // output
                  &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr, // output
                  &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,   // output
                  &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,     // output
                  &sol_shift, &sol_tshift,                          // output
                  &nre, re_vector.data(),
                  &nae, ae_vector.data(),
                  &n_jumps, jumps,
                  &n_thit, thit_exactly,
                  &n_direction, direction,
                  &n_isterminal, isterminal,
                  &opts_cc,
                  ef_cc, chng_cc, NULL, NULL);

  EXPECT_EQ(sol_npts, 46);

  EXPECT_FLOAT_EQ(sol_t_ptr[1]              , 0.005682091767);
  EXPECT_FLOAT_EQ(sol_y_ptr[1]              , 80.11319105   );
  EXPECT_FLOAT_EQ(sol_y_ptr[sol_npts + 1]   , 19.88680895   );
  EXPECT_FLOAT_EQ(sol_y_ptr[2*sol_npts + 1] , 0.000564993367);

  EXPECT_FLOAT_EQ(sol_t_ptr[sol_npts - 1]   , 20.0);
  EXPECT_FLOAT_EQ(sol_y_ptr[sol_npts - 1]   , 39.999859);
  EXPECT_FLOAT_EQ(sol_y_ptr[2*sol_npts - 1] , 60.000141);
  EXPECT_FLOAT_EQ(sol_y_ptr[3*sol_npts - 1] , 1.5000278);
}

struct dde_solver_ex445_test : public testing::Test {
  static int state;
  static constexpr double t0 = 0.0, tfinal = 12.0;
  static constexpr int nout = 1000;

  // output
  int sol_npts, sol_flag, sol_ne;
  double* sol_t_ptr;
  double* sol_y_ptr;
  double* sol_te_ptr;
  double* sol_ye_ptr;
  double* sol_queue_ptr;
  double* sol_yoft_ptr;
  double* sol_tqueue_ptr;
  int* sol_stats_ptr;
  int* sol_ie_ptr;
  int* sol_ipoint_ptr;
  bool sol_shift;
  double sol_tshift;

  std::vector<int> nvar{2, 1, 2};
  int n_nvar = nvar.size();
  std::vector<double> tspan;
  int n_tspan = nout;
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{1.E-5};
  std::vector<double> ae_vector{default_abs_err()};
  int nre = re_vector.size();
  int nae = ae_vector.size();

  /* we have to pass arrays to fortran, but n_xxx = 0
     indicates these arrays are dummies */
  int n_jumps = 0;
  int n_thit = 0;
  double jumps[1];
  double thit_exactly[1];
  int n_direction = 2;
  int direction[2] = { -1, 0 };
  int n_isterminal = 2;
  bool isterminal[2] = { false, true };

  std::vector<double> delay{ 0.1 };
  std::vector<double> history{ 0.0, 0.0 };

  dde_solver_ex445_test() {
    const double h = (tfinal - t0) / (nout - 1);
    for (int i = 0; i < nout; ++i) {
      tspan.push_back(t0 + i * h);
    }
  }

  static void ddes_cc(const double* t, const int* n, const int* nlags, const double y[], const double z[], double dy[]) {
    double ylag;
    static const double gamma=0.248, beta=1.0, A=0.75, omega=1.37;
    ylag = z[0];
    dy[0] = y[1];
    dy[1] = sin(y[0]) - state * gamma * cos(y[0]) - beta * ylag
      + A * sin(omega * (*t) + asin(gamma / A));
  }

    // event function
  static void ef_cc(const double* t, const int* n, const int* nlag,
                    const double y[], const double dy[], const double z[], double g[]) {
    double ylag = z[0];
    g[0] = y[0];
    g[1] = fabs(y[0]) - asin(1.0);
  }

  static void chng_cc(const int* nevent, const double* tevent,
                      double yevent[], double dyevent[], double* hinit,
                      const int* n_direction, int direction[],
                      const int* n_isterminal, bool isterminal[], bool* quit) {
    // Restart the integration with initial values
    // that correspond to a bounce of the suitcase.
    if (*nevent == 1) {
      state = -state;
      yevent[0] = 0.0;
      yevent[1] = 0.913 * yevent[1];
      direction[0] = - direction[0];
    }
  }
};

int dde_solver_ex445_test::state = 1;

TEST_F(dde_solver_ex445_test, sol) {

  state = 1;

  integrate_dde_4(&n_nvar, nvar.data(),
                  ddes_cc, &nvar[1], delay.data(), &nvar[0], history.data(),
                  &n_tspan, tspan.data(),
                  &sol_npts, &sol_flag, &sol_ne,                    // output
                  &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr, // output
                  &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,   // output
                  &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,     // output
                  &sol_shift, &sol_tshift,                          // output
                  &nre, re_vector.data(),
                  &nae, ae_vector.data(),
                  &n_jumps, jumps,
                  &n_thit, thit_exactly,
                  &n_direction, direction,
                  &n_isterminal, isterminal,
                  &opts_cc,
                  ef_cc, chng_cc, NULL, NULL);

  EXPECT_EQ(sol_npts, 972);
  EXPECT_EQ(sol_ne, 3);
  
  // event time
  EXPECT_FLOAT_EQ(sol_te_ptr[0], 4.5167570648757893);
  EXPECT_FLOAT_EQ(sol_te_ptr[1], 9.7510095674153696);
  EXPECT_FLOAT_EQ(sol_te_ptr[2], 11.670403419456539);
}
