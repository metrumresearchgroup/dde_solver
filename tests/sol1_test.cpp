#include <dde_solver.hpp>
#include <dde_options.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

struct dde_solver_b2_test : public testing::Test {
  static int state;
  int n_nvar = 2, n_tspan = 2;

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
  std::vector<double> tspan{0.0, 2.0 * log(66)};
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{1.0E-5};
  std::vector<double> ae_vector{1.0E-5};
  int nre = 1;
  int nae = 1;

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

  dde_solver_b2_test() {}

  static void ddes_cc(double* t, int* n, int* nlags, double y[], double z[], double dy[]) {
    double ylag = z[0];
    double f;
    if (ylag < 0.0) {
      f = 1.0;
    } else {
      f = -1.0;
    }
    dy[0] = f - y[0];
  }

  static void ddes2_cc(double* t, int* n, int* nlags, double y[], double z[], double dy[]) {
    double f;
    if (state < 0.0) {
      f = 1.0;
    } else {
      f = -1.0;
    }
    dy[0] = f - y[0];
  }

  static void history_cc(double* t, int* n, double y[]) {
    y[0] = 1.0;
  }

  static void beta_cc(double* t, int* n, int* nlags, double y[], double bval[]) {
    bval[0] = *t/2.0;
  }

  // event function
  static void ef_cc(double* t, int* n, int* nlag, double y[], double dy[], double z[], double g[]) {
    double ylag = z[0];
    g[0] = ylag;
  }

  // change function when an event occurs
  static void chng_cc(int* nevent, double* tevent,
                      double yevent[], double dyevent[], double* hinit,
                      int* n_direction, int direction[],
                      int* n_isterminal, bool isterminal[], bool* quit) {
    if (*nevent == 1) {
      state = -state;
    }
  }
};

int dde_solver_b2_test::state = 1;

TEST_F(dde_solver_b2_test, b2) {

  integrate_dde_1(&n_nvar, nvar.data(),
                  ddes_cc, beta_cc, history_cc, 
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

  EXPECT_EQ(sol_npts, 75);
  EXPECT_EQ(sol_flag, 0);       // solution successful

  // // Fortran solutions at the first/last 3 points
  std::vector<double> sol_t_head {0.0, 0.08379309484, 0.2513792845, 0.5865516639};
  std::vector<double> sol_t_tail {8.379279328, 8.37930049, 8.379309484};

  std::cout.precision(10);
  for (size_t i = 0; i < sol_t_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[i], sol_t_head[i]);
  }

  for (size_t i = 0; i < sol_t_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[sol_npts - sol_t_tail.size() + i], sol_t_tail[i]);
  }

  // // Fortran solutions at the first/last 3 points
  std::vector<double> sol_y_head {1.0, 0.8392430217, 0.5554546718, 0.112484232};
  std::vector<double> sol_y_tail {-0.9848479861, -0.9848483068, -0.9848355217};
  for (size_t i = 0; i < sol_y_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[i], sol_y_head[i]);
  }

  for (size_t i = 0; i < sol_y_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[sol_npts - sol_y_tail.size() + i], sol_y_tail[i]);
  }  
}

TEST_F(dde_solver_b2_test, b2g) {

  state = 1;
  n_nvar = 3;
  int nef = 1;
  nvar.push_back(nef);
  isterminal[0] = false;
  direction[0] = 0;
  n_direction = 1;
  n_isterminal = 1;

  integrate_dde_1(&n_nvar, nvar.data(),
                  ddes2_cc, beta_cc, history_cc, 
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

  EXPECT_EQ(sol_npts, 26);
  EXPECT_EQ(sol_flag, 0);       // solution successful

  // // Fortran solutions at the first/last 3 points
  std::vector<double> sol_t_head {0.0, 0.08379309484, 0.2513792845, 0.5865516639};
  std::vector<double> sol_t_tail {7.732039825, 8.356846033, 8.379309484};
  std::cout.precision(10);
  for (size_t i = 0; i < sol_t_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[i], sol_t_head[i]);
  }
  for (size_t i = 0; i < sol_t_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[sol_npts - sol_t_tail.size() + i], sol_t_tail[i]);
  }

  // Fortran solutions at the first/last 3 points
  std::vector<double> sol_y_head {1.0, 0.8392430217, 0.5554546718, 0.112484232};
  std::vector<double> sol_y_tail {-0.9710556888, -0.9845041444, -0.9848387599};
  for (size_t i = 0; i < sol_y_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[i], sol_y_head[i]);
  }

  for (size_t i = 0; i < sol_y_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[sol_npts - sol_y_tail.size() + i], sol_y_tail[i]);
  }

  // there are events occurred
  EXPECT_EQ(sol_ne, 3);
  std::vector<double> sol_te {1.38629378, 3.58351673, 8.379304687};
  std::vector<double> sol_ye {-0.4999996322, 0.8333329703, -0.9848482816};
  for (size_t i = 0; i < sol_ne; ++i) {
    EXPECT_FLOAT_EQ(sol_te_ptr[i], sol_te[i]);
    EXPECT_FLOAT_EQ(sol_ye_ptr[i], sol_ye[i]);
  }
}

struct dde_solver_c2_test : public testing::Test {
  static int state1, state2;
  int n_nvar = 3, n_tspan = 2;

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
  std::vector<double> tspan{0.0, 40.0};
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{1.0E-5};
  std::vector<double> ae_vector{1.0E-5};
  int nre = 1;
  int nae = 1;

  /* we have to pass arrays to fortran, but n_xxx = 0
     indicates these arrays are dummies */
  int n_jumps = 0;
  int n_thit = 0;
  int n_direction = 2;
  int n_isterminal = 2;
  double jumps[1];
  double thit_exactly[1];
  int direction[2] = { 0, 0 };
  bool isterminal[2] = { false, false };

  dde_solver_c2_test() {}

  static void ddes_cc(double* t, int* n, int* nlags, double y[], double z[], double dy[]) {
    double ylag1 = z[0], y1 = y[0];
    dy[0] = -2.0 * ylag1;
    if (state2 == -1) {
      y1 = -y1;      
    }
    if (state1 == -1) {
      ylag1 = -ylag1;      
    }
    dy[1] = (ylag1 - y1) / (1.0 + ylag1);
  }

  static void history_cc(double* t, int* n, double y[]) {
    y[0] = 1.0;
    y[1] = 0.5;
  }

  static void beta_cc(double* t, int* n, int* nlags, double y[], double bval[]) {
    bval[0] = *t - y[1];
  }

  // event function
  static void ef_cc(double* t, int* n, int* nlag, double y[], double dy[], double z[], double g[]) {
    double ylag1 = z[0];
    g[0] = ylag1;
    g[1] = y[0];
  }

  // change function when an event occurs
  static void chng_cc(int* nevent, double* tevent,
                      double yevent[], double dyevent[], double* hinit,
                      int* n_direction, int direction[],
                      int* n_isterminal, bool isterminal[], bool* quit) {
    if (*nevent == 1) {
      state1 = -state1;
    }
    if (*nevent == 2) {
      state2 = -state2;
    }    
  }
};

// derivative switches
int dde_solver_c2_test::state1 = 1;
int dde_solver_c2_test::state2 = 1;

TEST_F(dde_solver_c2_test, c2) {

  integrate_dde_1(&n_nvar, nvar.data(),
                  ddes_cc, beta_cc, history_cc, 
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

  EXPECT_EQ(sol_npts, 271);
  EXPECT_EQ(sol_flag, 0);       // solution successful

  // // Fortran solutions at the first/last 3 points
  std::vector<double> sol_t_head {0.0, 0.05, 0.15, 0.35};
  std::vector<double> sol_t_tail {39.70777588, 39.87330109, 40};
  for (size_t i = 0; i < sol_t_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[i], sol_t_head[i]);
  }
  for (size_t i = 0; i < sol_t_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[sol_npts - sol_t_tail.size() + i], sol_t_tail[i]);
  }

  // Fortran solutions at the first/last 3 points: y1
  std::vector<double> sol_y1_head {1.0, 0.9, 0.7, 0.3};
  std::vector<double> sol_y1_tail {0.595325, 0.470919, 0.33905858};
  for (size_t i = 0; i < sol_y1_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[i], sol_y1_head[i]);
  }
  for (size_t i = 0; i < sol_y1_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[sol_npts - sol_y1_tail.size() + i], sol_y1_tail[i]);
  }

  // Fortran solutions at the first/last 3 points: y2
  std::vector<double> sol_y2_head {0.5, 0.50125, 0.51125, 0.56125};
  std::vector<double> sol_y2_tail {0.6536570292, 0.6335469724, 0.6428966006};
  for (size_t i = 0; i < sol_y2_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[i + sol_npts], sol_y2_head[i]);
  }
  for (size_t i = 0; i < sol_y2_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol_y_ptr[2*sol_npts - sol_y2_tail.size() + i], sol_y2_tail[i]);
  }

  // events
  EXPECT_EQ(sol_ne, 56);
}
