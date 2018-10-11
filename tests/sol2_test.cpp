#include <dde_solver.hpp>
#include <dde_options.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

struct dde_solver_ex444_test : public testing::Test {
  static constexpr double R=3.5, M=19.0;

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

  std::vector<int> nvar{1, 1, 2};
  int n_nvar = nvar.size();
  std::vector<double> tspan{0.0, 40.0};
  int n_tspan = tspan.size();
  struct dde_opts_cc opts_cc = default_dde_opts();
  std::vector<double> re_vector{1.E-4};
  std::vector<double> ae_vector{1.E-7};
  int nre = re_vector.size();
  int nae = ae_vector.size();

  /* we have to pass arrays to fortran, but n_xxx = 0
     indicates these arrays are dummies */
  int n_jumps = 1;
  int n_thit = 0;
  int n_direction = 2;
  int n_isterminal = 0;
  double jumps[1] = {0.0};
  double thit_exactly[1];
  int direction[2] = {1, -1};
  bool isterminal[1];

  std::vector<double> delay{ 0.74 };

  dde_solver_ex444_test() {
    opts_cc.interpolation = true;
  }

  static void ddes_cc(double* t, int* n, int* nlags, double y[], double z[], double dy[]) {
    dy[0] = R * y[0] * (1.0 - z[0] / M);
  }

  static void history_cc(double* t, int* n, double y[]) {
    // Use the JUMPS option to tell solver about
    // the discontinuity at t = 0.
    if (*t == 0.0) {
      y[0] = 19.001;
    } else {
      y[0] = 19.0;      
    }
  }

  // event function
  static void ef_cc(double* t, int* n, int* nlag, double y[], double dy[], double z[], double g[]) {
    double ylag = z[0];
    g[0] = dy[0];
    g[1] = dy[0];
  }
};

TEST_F(dde_solver_ex444_test, sol) {

  integrate_dde_2(&n_nvar, nvar.data(),
                  ddes_cc, &nvar[1], delay.data(), history_cc,
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
                  ef_cc, NULL, NULL, NULL);

  EXPECT_EQ(sol_npts, 247);
  EXPECT_EQ(sol_ne, 29);
  EXPECT_EQ(sol_flag, 0);       // solution successful
  
  std::vector<double> sol_te{
      0.222      ,
      0.444      ,        
      0.666      ,        
      0.74       ,         
      1.769775904,  
      3.096061136,  
      4.408009079,  
      5.720544093,  
      7.032395192,  
      8.345671622,  
      9.657853266,  
      10.97088795,  
      12.28251841,  
      13.59696625,  
      14.90608014,  
      16.22553137,  
      17.52578094,  
      18.86458542,  
      20.13748323,  
      21.56558453,  
      22.77891694,  
      24.81088036,  
      25.99040395,  
      28.95463428,  
      30.13403615,  
      33.10305783,  
      34.28245981,  
      37.25148153,   
      38.43088351};

  std::vector<double> sol_ye {
      19.001       ,
      19.001       ,
      19.001       ,
      19.001       ,
      18.99791092  ,
      19.00394652  ,
      18.9925348   ,
      19.01413737  ,
      18.97325699  ,
      19.05069113  ,
      18.90428332  ,
      19.18205183  ,
      18.65847038  ,
      19.65756904  ,
      17.794963    ,
      21.42262739  ,
      14.92442709  ,
      28.48921772  ,
      7.391971634  ,
      58.96076945  ,
      0.3866781655 ,
      96.73055858  ,
      0.01190277783,
      96.92295087  ,
      0.0116982637 ,
      96.92295011  ,
      0.01169826332,
      96.92295025  ,
      0.01169826334 };

  for (int i = 0; i < sol_ne; ++i) {
    EXPECT_FLOAT_EQ(sol_te_ptr[i], sol_te[i]);
    EXPECT_FLOAT_EQ(sol_ye_ptr[i], sol_ye[i]);
  }
}
