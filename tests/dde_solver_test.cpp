#include <dde_solver.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

struct TestFn_b2 {
  static int state;

  struct FnDdes {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z) const {
      double ylag = z[0];
      double f;
      if (ylag < 0.0) {
        f = 1.0;
      } else {
        f = -1.0;
      }
      return { f - y[0] };
    }
  };

  struct FnDdes2 {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z) const {
      double f;
      if (state < 0) {
        f = 1.0;
      } else {
        f = -1.0;
      }
      return { f - y[0] };
    }
  };

  struct FnBeta {
    std::vector<double> operator()(const double& t,
                                   const int& nlags,
                                   const std::vector<double>& y) const {
      return { 0.50 * t };
    }
  };
  
  struct FnHistory {
    std::vector<double> operator()(const double& t) const {
      return { 1.0 };
    }
  };

  struct FnEvent
  {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& dydt,
                                   const std::vector<double>& z) const {
      return { z[0] };
    }  
  };

  struct FnChange
  {
    void operator()(const double& t, const double& h,
                    const std::vector<double>& y,
                    const std::vector<double>& dy,
                    const std::vector<int>& v,
                    const std::vector<bool>& i,
                    const bool& quit) const {
      if (y.size() == 1) state = -state;
    }
  };

  FnDdes ddes;
  FnDdes2 ddes2;
  FnBeta beta;
  FnHistory his;
  FnEvent ef;
  FnChange chng;
};

int TestFn_b2::state = 1;

TEST(dde_solver_test, b2) {
  using dde_solver::DdeUserOption;
  using dde_solver::DdeSol;
  using dde_solver::DdeIntegrator;

  TestFn_b2 test_fn;
  using Dde = DdeUserOption<TestFn_b2::FnDdes, TestFn_b2::FnBeta, TestFn_b2::FnHistory,
                            nullptr_t, nullptr_t, nullptr_t, nullptr_t>;
  Dde o1(test_fn.ddes, test_fn.beta, test_fn.his, nullptr, nullptr, nullptr, nullptr);

  EXPECT_EQ(!o1.f_event(), true);
  EXPECT_EQ(!o1.f_change(), true);
  EXPECT_EQ(!o1.f_out(), true);
  EXPECT_EQ(!o1.f_user(), true);

  o1.nvar[0] = 1;               // n
  o1.nvar[1] = 1;               // nlags

  // tspan
  o1.tspan[0] = 0.0;
  o1.tspan[1] = 2.0 * log(66);

  // set error tol
  o1.rerr[0] = 1.E-5;
  o1.aerr[0] = 1.E-5;

  DdeSol sol = DdeIntegrator<Dde>(o1)();
                
  EXPECT_EQ(sol.npts, 75);
  EXPECT_EQ(sol.flag, 0);       // solution successful

  // Fortran solutions at the first/last 3 points
  std::vector<double> sol_t_head {0.0, 0.08379309484, 0.2513792845, 0.5865516639};
  std::vector<double> sol_t_tail {8.379279328, 8.37930049, 8.379309484};
  std::cout.precision(10);
  for (size_t i = 0; i < sol_t_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.t[i], sol_t_head[i]);
  }
  for (size_t i = 0; i < sol_t_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.t[sol.npts - sol_t_tail.size() + i], sol_t_tail[i]);
  }

  // // Fortran solutions at the first/last 3 points
  std::vector<double> sol_y_head {1.0, 0.8392430217, 0.5554546718, 0.112484232};
  std::vector<double> sol_y_tail {-0.9848479861, -0.9848483068, -0.9848355217};
  for (size_t i = 0; i < sol_y_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.y[i], sol_y_head[i]);
  }
  for (size_t i = 0; i < sol_y_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.y[sol.npts - sol_y_tail.size() + i], sol_y_tail[i]);
  }
}

TEST(dde_solver_test, b2g) {
  using dde_solver::DdeUserOption;
  using dde_solver::DdeSol;
  using dde_solver::DdeIntegrator;
  
  TestFn_b2 test_fn;
  using Dde = DdeUserOption<TestFn_b2::FnDdes2, TestFn_b2::FnBeta, TestFn_b2::FnHistory,
                            TestFn_b2::FnEvent, TestFn_b2::FnChange>;
  Dde o1(test_fn.ddes2, test_fn.beta, test_fn.his, test_fn.ef, test_fn.chng, nullptr, nullptr);

  EXPECT_EQ(!o1.f_event(), false);
  EXPECT_EQ(!o1.f_change(), false);
  EXPECT_EQ(!o1.f_out(), true);
  EXPECT_EQ(!o1.f_user(), true);

  TestFn_b2::state = 1;
  o1.nvar.resize(3);
  o1.nvar[0] = 1;               // n
  o1.nvar[1] = 1;               // nlags
  o1.nvar[2] = 1;               // nef

  // tspan
  o1.tspan[0] = 0.0;
  o1.tspan[1] = 2.0 * log(66);

  o1.isterminal.push_back(0);
  o1.direction.push_back(0);

  // set error tol
  o1.rerr[0] = 1.E-5;
  o1.aerr[0] = 1.E-5;

  DdeSol sol = DdeIntegrator<Dde>(o1)();

  EXPECT_EQ(sol.npts, 26);
  EXPECT_EQ(sol.flag, 0);       // solution successful

  // Fortran solutions at the first/last 3 points
  std::vector<double> sol_t_head {0.0, 0.08379309484, 0.2513792845, 0.5865516639};
  std::vector<double> sol_t_tail {7.732039825, 8.356846033, 8.379309484};
  std::cout.precision(10);
  for (size_t i = 0; i < sol_t_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.t[i], sol_t_head[i]);
  }
  for (size_t i = 0; i < sol_t_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.t[sol.npts - sol_t_tail.size() + i], sol_t_tail[i]);
  }

  // Fortran solutions at the first/last 3 points
  std::vector<double> sol_y_head {1.0, 0.8392430217, 0.5554546718, 0.112484232};
  std::vector<double> sol_y_tail {-0.9710556888, -0.9845041444, -0.9848387599};
  for (size_t i = 0; i < sol_y_head.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.y[i], sol_y_head[i]);
  }

  for (size_t i = 0; i < sol_y_tail.size(); ++i) {
    EXPECT_FLOAT_EQ(sol.y[sol.npts - sol_y_tail.size() + i], sol_y_tail[i]);
  }

  // there are events occurred
  EXPECT_EQ(sol.ne, 3);
  std::vector<double> sol_te {1.38629378, 3.58351673, 8.379304687};
  std::vector<double> sol_ye {-0.4999996322, 0.8333329703, -0.9848482816};
  for (size_t i = 0; i < sol.ne; ++i) {
    EXPECT_FLOAT_EQ(sol.te[i], sol_te[i]);
    EXPECT_FLOAT_EQ(sol.ye[i], sol_ye[i]);
  }
}

TEST(dde_solver_test, ex444) {
  using dde_solver::DdeUserOption;
  using dde_solver::DdeSol;
  using dde_solver::DdeIntegrator;
  
  static constexpr double R=3.5, M=19.0;

  struct FnDdes {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z) const {
      return { R * y[0] * (1.0 - z[0] / M) };
    }
  };

  struct FnHistory {
    std::vector<double> operator()(const double& t) const {
      if (t == 0.0) {
        return {19.001};
      } else {
        return {19.0};
      }
    }
  };

  struct FnEvent
  {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& dydt,
                                   const std::vector<double>& z) const {
      return { dydt[0], dydt[0] };
    }
  };

  FnDdes ddes;
  FnHistory his;
  FnEvent ef;
  std::vector<double> delay{ 0.74 };
  using Dde = DdeUserOption<FnDdes, double, FnHistory,
                            FnEvent>;
  Dde o1(ddes, delay, his, ef, nullptr, nullptr, nullptr);

  o1.nvar.resize(3);
  o1.nvar[0] = 1;               // n
  o1.nvar[1] = 1;               // nlags
  o1.nvar[2] = 2;               // nef

  // tspan
  o1.tspan[0] = 0.0;
  o1.tspan[1] = 40.0;

  o1.jumps.push_back(0.0);
  o1.direction.push_back(1);
  o1.direction.push_back(-1);

  // set error tol
  o1.rerr[0] = 1.E-4;
  o1.aerr[0] = 1.E-7;

  o1.opts.interpolation = true;

  DdeSol sol = DdeIntegrator<Dde>(o1)();

  EXPECT_EQ(sol.npts, 247);
  EXPECT_EQ(sol.ne, 29);
  EXPECT_EQ(sol.flag, 0);       // solution successful
  EXPECT_EQ(sol.ye.size(), sol.ne);
  EXPECT_EQ(sol.te.size(), sol.ne);

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

  for (int i = 0; i < sol.ne; ++i) {
    EXPECT_FLOAT_EQ(sol.te[i], sol_te[i]);
    EXPECT_FLOAT_EQ(sol.ye[i], sol_ye[i]);
  }
}

TEST(dde_solver_test, exsd1) {
  using dde_solver::DdeUserOption;
  using dde_solver::DdeSol;
  using dde_solver::DdeIntegrator;
  
  struct FnDdes {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z) const {
      return { y[0] * z[0] / t };
    }
  };

  struct FnBeta {
    std::vector<double> operator()(const double& t,
                                   const int& nlags,
                                   const std::vector<double>& y) const {
      return { log(y[0]) };
    }
  };

  FnDdes ddes;
  FnBeta beta;
  std::vector<double> history{ 1.0 };
  using Dde = DdeUserOption<FnDdes, FnBeta, double>;
  Dde o1(ddes, beta, history, nullptr, nullptr, nullptr, nullptr);

  o1.nvar[0] = 1;               // n
  o1.nvar[1] = 1;               // nlags

  // tspan
  o1.tspan[0] = 1.0;
  o1.tspan[1] = 10.0;

  o1.jumps.push_back(0.0);
  o1.direction.push_back(1);
  o1.direction.push_back(-1);

  // set error tol
  o1.rerr[0] = 1.E-5;
  o1.aerr[0] = 1.E-5;

  DdeSol sol = DdeIntegrator<Dde>(o1)();

  EXPECT_EQ(sol.npts, 21);
  EXPECT_EQ(sol.flag, 0);       // solution successful

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

  for (int i = 0; i < sol.npts; ++i) {
    EXPECT_FLOAT_EQ(sol.t[i], sol_t[i]);
    EXPECT_FLOAT_EQ(sol.y[i], sol_y[i]);
  }
}

TEST(dde_solver_test, ex441) {
  using dde_solver::DdeUserOption;
  using dde_solver::DdeSol;
  using dde_solver::DdeIntegrator;
  
  struct FnDdes {
    std::vector<double> operator()(const double& t,
                                   const std::vector<double>& y,
                                   const std::vector<double>& z) const {
      const double tau = 42.0, omega = 0.15;
      const double a = 0.330, d = 0.006, lambda = 0.308, gamma = 0.040 , epsilon=0.060;

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
      return { dsdt, dedt, didt, drdt };
    }
  };

  FnDdes ddes;
  std::vector<double> history{15.0, 0.0, 2.0, 3.0};
  std::vector<double> delay{ 42.0, 0.15 };
  using Dde = DdeUserOption<FnDdes, double, double>;
  Dde o1(ddes, delay, history, nullptr, nullptr, nullptr, nullptr);

  o1.nvar[0] = 4;               // n
  o1.nvar[1] = 2;               // nlags

  // tspan
  o1.tspan[0] = 0.0;
  o1.tspan[1] = 350.0;

  // set error tol
  o1.rerr[0] = default_rel_err();
  o1.aerr[0] = default_abs_err();

  DdeSol sol = DdeIntegrator<Dde>(o1)();

  EXPECT_EQ(sol.flag, 0);       // solution successful
  EXPECT_FLOAT_EQ(sol.t[sol.npts - 1], o1.tspan[1]);

  EXPECT_FLOAT_EQ(sol.y[sol.npts - 1]   , 5.231271201   );
  EXPECT_FLOAT_EQ(sol.y[2*sol.npts - 1] , 0.05490837936 );
  EXPECT_FLOAT_EQ(sol.y[3*sol.npts - 1] , 3.985113815   );
  EXPECT_FLOAT_EQ(sol.y[4*sol.npts - 1] , 5.915635409   );
}
