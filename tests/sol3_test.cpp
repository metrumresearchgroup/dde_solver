#include <dde_solver.hpp>
#include <dde_options.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

struct dde_solver_exsd1_test : public testing::Test {
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
  std::vector<double> tspan{1.0, 10.0};
  int n_tspan = tspan.size();
  double history[1] = { 1.0 };
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{1.E-5};
  std::vector<double> ae_vector{1.E-5};
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

  dde_solver_exsd1_test() {}

  static void ddes_cc(double* t, int* n, int* nlags, double y[], double z[], double dy[]) {
    dy[0] = y[0] * z[0] / (*t);
  }

  static void beta_cc(double* t, int* n, int* nlags, double y[], double bval[]) {
    bval[0] = log(y[0]);
  }
};

TEST_F(dde_solver_exsd1_test, sol) {

  integrate_dde_3(&n_nvar, nvar.data(),
                  ddes_cc, beta_cc, &nvar[0], history,
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

  EXPECT_EQ(sol_npts, 21);
  EXPECT_EQ(sol_flag, 0);       // solution successful
  
  std::vector<double> sol_t {
      1          ,  
      1.09       ,  
      1.27       ,  
      1.63       ,  
      2.35       ,  
      2.575      ,  
      2.6875     ,  
      2.808290335,  
      2.988290335,  
      3.348290335,  
      4.039252256,  
      4.790665749,  
      5.579446823,  
      6.398737789,  
      7.243951735,  
      7.479083337,  
      7.659083337,  
      8.019083337,  
      8.739083337,  
      9.62198824 ,  
      10         };

  std::vector<double> sol_y {
    1          ,
    1.09       ,
    1.27       ,
    1.63       ,
    2.35       ,
    2.575      ,
    2.6875     ,
    2.809797119,
    3.00215565 ,
    3.427280805,
    4.419196104,
    5.826329976,
    7.787839984,
    10.52721425,
    14.36651937,
    15.6645679 ,
    16.73731622,
    19.11227919,
    24.96810772,
    34.87223264,
    40.36169238 };

  for (int i = 0; i < sol_npts; ++i) {
    EXPECT_FLOAT_EQ(sol_t_ptr[i], sol_t[i]);
    EXPECT_FLOAT_EQ(sol_y_ptr[i], sol_y[i]);
  }
}
