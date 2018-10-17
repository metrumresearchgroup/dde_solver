#ifndef DDE_SOLVER_HPP
#define DDE_SOLVER_HPP

#include <iostream>
#include <sstream>
#include <vector>
#include <string>

extern "C" {

#include "dde_solver.h"

  // fuction pointers for C++ interface, only difference
  // being the attached @c user_data for callbacks
  typedef void (*f_ddes_fcn_cc) (const double*, const int*,
                                 const int*, const double[],
                                 const double[], double[], void* user_data);

  typedef void (*f_beta_fcn_cc) (const double*, const int*, const int*, const double[],
                                double[], void* user_data);

  typedef void (*f_history_fcn_cc) (const double*, const int*, double[], void* user_data);

  typedef void (*f_event_fcn_cc) (const double*, const int*,
                                 const int*, const double[], const double[],
                                 const double[], double[], void* user_data);

  typedef void (*f_change_fcn_cc) (const int*, const double*, 
                                  double[], double[], double*, const int*, int[],
                                  const int*, bool[], bool*, void* user_data);

  typedef void (*f_out_fcn_cc) (double* t, double[], double[], int*, int*, void* user_data);

  typedef void (*f_user_trim_get_cc) (void* user_data);

  void integrate_dde_1_cc(void* user_data,
                          const int* n_nvar, const int nvar[],
                          f_ddes_fcn_cc,
                          f_beta_fcn_cc,
                          f_history_fcn_cc,
                          const int* n_tspan, double tspan[],
                          // output
                          int* sol_npts, int* sol_nlags, int* sol_ne,
                          double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                          double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                          int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                          bool* sol_shift, double* sol_tshift,
                          // options
                          const int* n_re_vector, double re_vector[],
                          const int* n_ae_vector, double ae_vector[],
                          const int* n_jumps, const double jumps[],
                          const int* n_thit, const double thit_exactly[],
                          const int* n_direction, const int direction[],
                          const int* n_isterminal, int isterminal[],
                          const struct dde_opts_cc* opts_cc,
                          f_event_fcn_cc, f_change_fcn_cc, f_out_fcn_cc, f_user_trim_get_cc);

  void integrate_dde_2_cc(void* user_data,
                          const int* n_nvar, const int nvar[],
                          f_ddes_fcn_cc,
                          const int* n_delay, const double delay[],
                          f_history_fcn_cc,
                          const int* n_tspan, double tspan[],
                          // output
                          int* sol_npts, int* sol_nlags, int* sol_ne,
                          double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                          double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                          int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                          bool* sol_shift, double* sol_tshift,
                          // options
                          const int* n_re_vector, double re_vector[],
                          const int* n_ae_vector, double ae_vector[],
                          const int* n_jumps, const double jumps[],
                          const int* n_thit, const double thit_exactly[],
                          const int* n_direction, const int direction[],
                          const int* n_isterminal, int isterminal[],
                          const struct dde_opts_cc* opts_cc,
                          f_event_fcn_cc, f_change_fcn_cc, f_out_fcn_cc, f_user_trim_get_cc);

  void integrate_dde_3_cc(void* user_data,
                          const int* n_nvar, const int nvar[],
                          f_ddes_fcn_cc,
                          f_beta_fcn_cc,
                          const int* n_his, const double history[],
                          const int* n_tspan, double tspan[],
                          // output
                          int* sol_npts, int* sol_nlags, int* sol_ne,
                          double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                          double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                          int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                          bool* sol_shift, double* sol_tshift,
                          // options
                          const int* n_re_vector, double re_vector[],
                          const int* n_ae_vector, double ae_vector[],
                          const int* n_jumps, const double jumps[],
                          const int* n_thit, const double thit_exactly[],
                          const int* n_direction, const int direction[],
                          const int* n_isterminal, int isterminal[],
                          const struct dde_opts_cc* opts_cc,
                          f_event_fcn_cc, f_change_fcn_cc, f_out_fcn_cc, f_user_trim_get_cc);

  void integrate_dde_4_cc(void* user_data,
                          const int* n_nvar, const int nvar[],
                          f_ddes_fcn_cc,
                          const int* n_delay, const double delay[],
                          const int* n_his, const double history[],
                          const int* n_tspan, double tspan[],
                          // output
                          int* sol_npts, int* sol_nlags, int* sol_ne,
                          double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                          double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                          int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                          bool* sol_shift, double* sol_tshift,
                          // options
                          const int* n_re_vector, double re_vector[],
                          const int* n_ae_vector, double ae_vector[],
                          const int* n_jumps, const double jumps[],
                          const int* n_thit, const double thit_exactly[],
                          const int* n_direction, const int direction[],
                          const int* n_isterminal, int isterminal[],
                          const struct dde_opts_cc* opts_cc,
                          f_event_fcn_cc, f_change_fcn_cc, f_out_fcn_cc, f_user_trim_get_cc);
}

namespace dde_solver {
  /*
   * map input @c beta or @c history type. When input is a functor, return itself.
   * @tparam F type of @c F_beta or @c F_his
   */
  template<typename F>
  struct DdeUserData {
    using type = F;
  };

  /*
   * Instead of functor, @c F_beta and @c F_his can also be arranged data
   */
  template<>
  struct DdeUserData<double> {
    using type = std::vector<double>;
  };

  /*
   * In @c DdeUserOption, depend on the type of @c F_event,
   * we need to either return a function pointer/lambda, or
   * return @c NULL, as Fortran's optional argument. This
   * helper function generates the return.
   */
  template<typename Dde, typename F>
  struct DdeEventFunc {
    f_event_fcn_cc operator()() const {
      return [](const double* t, const int* n, const int* nlags,
                const double y[], const double dydt[], const double z[],
                double g[], void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        std::vector<double> yv(y, y + *n);
        std::vector<double> dydtv(dydt, dydt + *n);
        std::vector<double> zv(z, z + (*n) * (*nlags));
        std::vector<double> gv = dde -> ef_cc(*t, yv, dydtv, zv);
        for (int i = 0; i < gv.size(); ++i) g[i] = gv[i];
      };
    }
  };
  
  template<typename Dde>
  struct DdeEventFunc<Dde, std::nullptr_t> {
    f_event_fcn_cc operator()() {
      return NULL;
    }
  };

  /*
   * Generate a function pointer for @c CHANGE_FCN according to the
   * template parameter type. @c user_data in the argument
   * list of the generated function will be used for callback,
   * and points to an object of @c DdeUserOption.
   *
   * @return a function pointer that can be passed to Fortran
   * function call @c integrate_dde_x_cc()
   */
  template<typename Dde, typename F>
  struct DdeChangeFunc {
    f_change_fcn_cc operator()() {
      return [](const int* nevent, const double* tevent,
                double yevent[], double dyevent[],
                double* hinit,
                const int* n_direction, int direction[],
                const int* n_isterminal, bool isterminal[],
                bool* quit, void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        double t = *tevent;
        double h = *hinit;
        std::vector<double> yv(yevent, yevent + *nevent);
        std::vector<double> dy(dyevent, dyevent + *nevent);
        std::vector<int> dv(direction, direction + *n_direction);
        std::vector<bool> iv(isterminal, isterminal + *n_isterminal);
        dde -> chng_cc(t, h, yv, dy, dv, iv, *quit);
      };
    }
  };

  template<typename Dde>
  struct DdeChangeFunc<Dde, std::nullptr_t> {
    f_change_fcn_cc operator()() {
      return NULL;
    }
  };

  /*
   * Generate a function pointer for @c OUT_FNC according to the
   * template parameter type. @c user_data in the argument
   * list of the generated function will be used for callback,
   * and points to an object of @c DdeUserOption.
   *
   * @return a function pointer that can be passed to Fortran
   * function call @c integrate_dde_x_cc()
   */
  template<typename Dde, typename F>
  struct DdeOutFunc {
    f_out_fcn_cc operator()() {
      return [](const double* t, const double y[], const double dy[],
                const int* n, const int* nevent, void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        std::vector<double> yv(y, y + *n);
        std::vector<double> dyv(dy, dy + *n);
        dde -> out_cc(*t, yv, dyv);
      };
    }
  };

  template<typename Dde>
  struct DdeOutFunc<Dde, std::nullptr_t> {
    f_out_fcn_cc operator()() {
      return NULL;
    }
  };

  /*
   * Generate a function pointer for @c USER_TRIM_GET_FNC according to the
   * template parameter type. @c user_data in the argument
   * list of the generated function will be used for callback,
   * and points to an object of @c DdeUserOption.
   *
   * @return a function pointer that can be passed to Fortran
   * function call @c integrate_dde_x_cc()
   */
  template<typename Dde, typename F>
  struct DdeUserTrimFunc {
    f_user_trim_get_cc operator()() {
      return [](void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        dde -> user_cc();
      };
    }
  };

  template<typename Dde>
  struct DdeUserTrimFunc<Dde, std::nullptr_t> {
    f_user_trim_get_cc operator()() {
      return NULL;
    }
  };

  /*
   * User options for DDE solver. Note that @c istermial is
   * of type @c int instead of @c bool to match Fortran's @c
   * logical. This is because C++'s @c vector<bool> is of
   * proxy type and requires special treatment. We avoid
   * this by using @c int and implicit cast between @c int
   * and @c bool.
   * @tparam F_ddes C++ functor with signature f(real, vec, vec) -> vec
   * @tparam F_beta C++ functor with signature f(real, int, vec) -> vec
   * @tparam F_his C++ functor with signature f(real) -> vec
   * @tparam F_event C++ functor with signature f(real, vec, vec, vec) -> vec
   * @tparam F_change C++ functor with signature f(real, real, vec, vec, vec_i, vec_bool, bool);
   * @tparam F_out C++ functor with signature f(real);
   * @tparam F_user C++ functor with signature f();
   */
  template<typename F_ddes, typename F_beta, typename F_his,
           typename F_event = std::nullptr_t,
           typename F_change = std::nullptr_t,
           typename F_out = std::nullptr_t,
           typename F_user = std::nullptr_t>
  struct DdeUserOption {
    static constexpr double default_rel_err = 1.0E-3;
    static constexpr double default_abs_err = 1.0E-6;

    dde_opts_cc opts;
    std::vector<int> nvar;
    std::vector<double> tspan;
    const F_ddes& ddes_cc;
    const F_beta& beta_cc;
    const F_his& history_cc;
    // const typename DdeUserData<F_beta>::type& beta_cc;
    // const typename DdeUserData<F_his>::type& history_cc;
    std::vector<double> rerr;
    std::vector<double> aerr;
    std::vector<double> jumps;
    std::vector<double> thit_exactly;
    std::vector<int> direction;
    std::vector<int> isterminal;
    const F_event& ef_cc;
    const F_change& chng_cc;
    const F_out& out_cc;
    const F_user& user_cc;
    DdeUserOption(const F_ddes& ddes,
                  const typename DdeUserData<F_beta>::type& beta,
                  const typename DdeUserData<F_his>::type& his,
                  const F_event& ef,
                  const F_change& chng,
                  const F_out& out,
                  const F_user& user) :
      opts       (default_dde_opts()),
      nvar       {0, 0},
      tspan      {0, 1},
      ddes_cc    (ddes ),
      beta_cc    (beta ),
      history_cc (his  ),
      rerr       {default_rel_err},
      aerr       {default_abs_err},
      ef_cc      (ef   ),
      chng_cc    (chng ),
      out_cc     (out  ),
      user_cc    (user )
    {}

    using Dde = DdeUserOption<F_ddes, F_beta, F_his, F_event, F_change, F_out, F_user>;

    /*
     * Generate a function pointer for @c ddes according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_ddes_fcn_cc f_ddes() {
      return [](const double* t, const int* n, const int* nlags,
                const double y[], const double z[], double dy[], void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        std::vector<double> yv(y, y + *n);
        std::vector<double> zv(z, z + (*n) * (*nlags));
        std::vector<double> dyv = dde -> ddes_cc(*t, yv, zv);
        for (int i = 0; i < *n; ++i) {
          dy[i] = dyv[i];
        }
      };
    }

    /*
     * Generate a function pointer for @c beta according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_beta_fcn_cc f_beta() {
      return [](const double* t, const int* n, const int* nlags,
                const double y[], double bval[], void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        std::vector<double> yv(y, y + *n);
        std::vector<double> bv = dde -> beta_cc(*t, *nlags, yv);
        for (int i = 0; i < *nlags; ++i) bval[i] = bv[i];
      };
    }

    /*
     * Generate a function pointer for @c history according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_history_fcn_cc f_history() {
      return [](const double* t, const int* n, double y[], void* user_data) {
        Dde* dde = static_cast<Dde*>(user_data);
        std::vector<double> yv = dde -> history_cc(*t);
        for (int i = 0; i < *n; ++i) y[i] = yv[i];
      };
    }

    /*
     * Generate a function pointer for @c EVENT_FNC according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_event_fcn_cc f_event() {
      return DdeEventFunc<Dde, F_event>()();
    }

    /*
     * Generate a function pointer for @c CHANGE_FCN according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_change_fcn_cc f_change() {
      return DdeChangeFunc<Dde, F_change>()();
    }

    /*
     * Generate a function pointer for @c OUT_FNC according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_out_fcn_cc f_out() {
      return DdeOutFunc<Dde, F_out>()();
    }

    /*
     * Generate a function pointer for @c USER_TRIM_GET_FNC according to the
     * template parameter type. @c user_data in the argument
     * list of the generated function will be used for callback,
     * and points to an object of @c DdeUserOption.
     *
     * @return a function pointer that can be passed to Fortran
     * function call @c integrate_dde_x_cc()
     */
    static f_user_trim_get_cc f_user() {
      return DdeUserTrimFunc<Dde, F_user>()();
    }
  };

  struct DdeSol {
    static constexpr int lipoint = 27; // size of ipoint
    static constexpr int lstats = 6;   // size of stats
    const int npts;                      // nb. of points
    const int flag;                      // nb. of delays
    const int ne;                        // nb. of events
    std::vector<double> t;               // time at points
    std::vector<double> y;               // solution at points
    std::vector<double> te;              // time at events
    std::vector<double> ye;              // solution at events
    // std::vector<double> queue;
    std::vector<double> yoft;
    // std::vector<double> tqueue;
    std::vector<int> stats;
    std::vector<int> ie;
    std::vector<int> ipoint;
    bool shift;
    double tshift;

    DdeSol(const int neqn, const bool& interpolation,
           const int& sol_npts, const int& sol_flag, const int& sol_ne,
           double* sol_t_ptr, double* sol_y_ptr, double* sol_te_ptr, double* sol_ye_ptr,
           double* sol_queue_ptr, double* sol_yoft_ptr, double* sol_tqueue_ptr,
           int* sol_stats_ptr, int* sol_ie_ptr, int* sol_ipoint_ptr,
           const bool& sol_shift, const double& sol_tshift) :
      npts(sol_npts), flag(sol_flag), ne(sol_ne),
      t(npts), y(npts * neqn), te(ne), ye(ne * neqn),
      yoft(neqn), stats(lstats), ipoint(lipoint),
      shift(sol_shift), tshift(sol_tshift)
    {
      for (int i = 0; i < npts; ++i)        t[i] = sol_t_ptr[i];
      for (int i = 0; i < npts * neqn; ++i) y[i] = sol_y_ptr[i];
      for (int i = 0; i < ne; ++i)          te[i] = sol_te_ptr[i];
      for (int i = 0; i < ne * neqn; ++i)   ye[i] = sol_ye_ptr[i];
      if(interpolation) {
        for (int i = 0; i < neqn; ++i)      yoft[i] = sol_yoft_ptr[i];
        for (int i = 0; i < lipoint; ++i)   ipoint[i] = sol_ipoint_ptr[i];
      }
      for (int i = 0; i < lstats; ++i)      stats[i] = sol_stats_ptr[i];
      for (int i = 0; i < lipoint; ++i)     ipoint[i] = sol_ipoint_ptr[i];
    }
  };

  /*
   * Integrator for DDE.
   *
   * @tparam Dde type of @c DdeUserOptions.
   */
  template<typename Dde>
  struct DdeIntegrator;

  /*
   * When @c F_ddes, @c F_beta, and @c F_his are all
   * functors, we call @c integrate_dde_1_cc()
   */
  template<typename F_ddes, typename F_beta, typename F_his,
           typename F_event, typename F_change, typename F_out, typename F_user>
  struct DdeIntegrator<DdeUserOption<F_ddes, F_beta, F_his,
                                     F_event, F_change, F_out, F_user> >
  {
    using Dde = DdeUserOption<F_ddes, F_beta, F_his, F_event, F_change, F_out, F_user>;
    Dde& dde;

    DdeIntegrator(Dde& dde_in) : dde(dde_in) {}

    DdeSol operator()() {

      void* user_data = static_cast<void*>(&dde);
      const int n_nvar = dde.nvar.size();
      const int n_tspan = dde.tspan.size();
      const int nre = dde.rerr.size();
      const int nae = dde.aerr.size();

      // some fortran compilers don't like empty arrays, so if
      // an array is of size 0, we push in a dummy value.
      const int n_jumps = dde.jumps.size();
      if (dde.jumps.empty()) dde.jumps.push_back(0);
      const int n_thit = dde.thit_exactly.size();
      if (dde.thit_exactly.empty()) dde.thit_exactly.push_back(0);
      const int n_direction = dde.direction.size();
      if (dde.direction.empty()) dde.direction.push_back(0);
      const int n_isterminal = dde.isterminal.size();
      if (dde.isterminal.empty()) dde.isterminal.push_back(0);

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

      integrate_dde_1_cc(user_data,
                         &n_nvar, dde.nvar.data(),
                         Dde::f_ddes(), Dde::f_beta(), Dde::f_history(),
                         &n_tspan, dde.tspan.data(),
                         // output
                         &sol_npts, &sol_flag, &sol_ne,
                         &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr,
                         &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,
                         &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,
                         &sol_shift, &sol_tshift,
                         // options
                         &nre, dde.rerr.data(),
                         &nae, dde.aerr.data(),
                         &n_jumps, dde.jumps.data(),
                         &n_thit, dde.thit_exactly.data(),
                         &n_direction, dde.direction.data(),
                         &n_isterminal, dde.isterminal.data(),
                         &dde.opts,
                         Dde::f_event(), Dde::f_change(),
                         Dde::f_out(), Dde::f_user());

      DdeSol sol(dde.nvar[0], dde.opts.interpolation,
                 sol_npts, sol_flag, sol_ne,
                 sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,
                 sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,
                 sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,
                 sol_shift, sol_tshift);

      return sol;
    }    
  };

  /*
   * For constant delays, @c F_beta will be the arrary of delays.
   * We call @c integrate_dde_2_cc() in this case.
   */
  template<typename F_ddes, typename F_his,
           typename F_event, typename F_change, typename F_out, typename F_user>
  struct DdeIntegrator<DdeUserOption<F_ddes, std::vector<double>, F_his,
                                     F_event, F_change, F_out, F_user> >
  {
    using Dde = DdeUserOption<F_ddes, std::vector<double>, F_his,
                              F_event, F_change, F_out, F_user>;

    Dde& dde;

    DdeIntegrator(Dde& dde_in) : dde(dde_in) {}

    DdeSol operator()() {

      void* user_data = static_cast<void*>(&dde);
      const int n_nvar = dde.nvar.size();
      const int n_tspan = dde.tspan.size();
      const int nre = dde.rerr.size();
      const int nae = dde.aerr.size();

      // some fortran compilers don't like empty arrays, so if
      // an array is of size 0, we push in a dummy value.
      const int n_jumps = dde.jumps.size();
      if (dde.jumps.empty()) dde.jumps.push_back(0);
      const int n_thit = dde.thit_exactly.size();
      if (dde.thit_exactly.empty()) dde.thit_exactly.push_back(0);
      const int n_direction = dde.direction.size();
      if (dde.direction.empty()) dde.direction.push_back(0);
      const int n_isterminal = dde.isterminal.size();
      if (dde.isterminal.empty()) dde.isterminal.push_back(0);

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

      const int n_beta = dde.beta_cc.size();

      integrate_dde_2_cc(user_data,
                         &n_nvar, dde.nvar.data(),
                         Dde::f_ddes(), &n_beta, dde.beta_cc.data(), Dde::f_history(),
                         &n_tspan, dde.tspan.data(),
                         // output
                         &sol_npts, &sol_flag, &sol_ne,
                         &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr,
                         &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,
                         &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,
                         &sol_shift, &sol_tshift,
                         // options
                         &nre, dde.rerr.data(),
                         &nae, dde.aerr.data(),
                         &n_jumps, dde.jumps.data(),
                         &n_thit, dde.thit_exactly.data(),
                         &n_direction, dde.direction.data(),
                         &n_isterminal, dde.isterminal.data(),
                         &dde.opts,
                         Dde::f_event(), Dde::f_change(),
                         Dde::f_out(), Dde::f_user());

      DdeSol sol(dde.nvar[0], dde.opts.interpolation,
                 sol_npts, sol_flag, sol_ne,
                 sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,
                 sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,
                 sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,
                 sol_shift, sol_tshift);

      return sol;
    }
  };

  /*
   * For constant history, @c F_his will be the arrary of history.
   * We call @c integrate_dde_3_cc() in this case.
   */
  template<typename F_ddes, typename F_beta,
           typename F_event, typename F_change, typename F_out, typename F_user>
  struct DdeIntegrator<DdeUserOption<F_ddes, F_beta, std::vector<double>,
                                     F_event, F_change, F_out, F_user> >
  {
    using Dde = DdeUserOption<F_ddes, F_beta, std::vector<double>,
                              F_event, F_change, F_out, F_user>;

    Dde& dde;

    DdeIntegrator(Dde& dde_in) : dde(dde_in) {}

    DdeSol operator()() {

      void* user_data = static_cast<void*>(&dde);
      const int n_nvar = dde.nvar.size();
      const int n_tspan = dde.tspan.size();
      const int nre = dde.rerr.size();
      const int nae = dde.aerr.size();

      // some fortran compilers don't like empty arrays, so if
      // an array is of size 0, we push in a dummy value.
      const int n_jumps = dde.jumps.size();
      if (dde.jumps.empty()) dde.jumps.push_back(0);
      const int n_thit = dde.thit_exactly.size();
      if (dde.thit_exactly.empty()) dde.thit_exactly.push_back(0);
      const int n_direction = dde.direction.size();
      if (dde.direction.empty()) dde.direction.push_back(0);
      const int n_isterminal = dde.isterminal.size();
      if (dde.isterminal.empty()) dde.isterminal.push_back(0);

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

      const int n_his = dde.history_cc.size();

      integrate_dde_3_cc(user_data,
                         &n_nvar, dde.nvar.data(),
                         Dde::f_ddes(), Dde::f_beta(), &n_his, dde.history_cc.data(),
                         &n_tspan, dde.tspan.data(),
                         // output
                         &sol_npts, &sol_flag, &sol_ne,
                         &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr,
                         &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,
                         &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,
                         &sol_shift, &sol_tshift,
                         // options
                         &nre, dde.rerr.data(),
                         &nae, dde.aerr.data(),
                         &n_jumps, dde.jumps.data(),
                         &n_thit, dde.thit_exactly.data(),
                         &n_direction, dde.direction.data(),
                         &n_isterminal, dde.isterminal.data(),
                         &dde.opts,
                         Dde::f_event(), Dde::f_change(),
                         Dde::f_out(), Dde::f_user());

      DdeSol sol(dde.nvar[0], dde.opts.interpolation,
                 sol_npts, sol_flag, sol_ne,
                 sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,
                 sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,
                 sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,
                 sol_shift, sol_tshift);

      return sol;
    }
  };

  /*
   * Use @c integrate_dde_4_cc() when both delay and history
   * are arrays.
   */
  template<typename F_ddes,
           typename F_event, typename F_change, typename F_out, typename F_user>
  struct DdeIntegrator<DdeUserOption<F_ddes, std::vector<double>, std::vector<double>,
                                     F_event, F_change, F_out, F_user> >
  {
    using Dde = DdeUserOption<F_ddes, std::vector<double>, std::vector<double>,
                              F_event, F_change, F_out, F_user>;

    Dde& dde;

    DdeIntegrator(Dde& dde_in) : dde(dde_in) {}

    DdeSol operator()() {

      void* user_data = static_cast<void*>(&dde);
      const int n_nvar = dde.nvar.size();
      const int n_tspan = dde.tspan.size();
      const int nre = dde.rerr.size();
      const int nae = dde.aerr.size();

      // some fortran compilers don't like empty arrays, so if
      // an array is of size 0, we push in a dummy value.
      const int n_jumps = dde.jumps.size();
      if (dde.jumps.empty()) dde.jumps.push_back(0);
      const int n_thit = dde.thit_exactly.size();
      if (dde.thit_exactly.empty()) dde.thit_exactly.push_back(0);
      const int n_direction = dde.direction.size();
      if (dde.direction.empty()) dde.direction.push_back(0);
      const int n_isterminal = dde.isterminal.size();
      if (dde.isterminal.empty()) dde.isterminal.push_back(0);

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

      const int n_beta = dde.beta_cc.size();
      const int n_his = dde.history_cc.size();

      integrate_dde_4_cc(user_data,
                         &n_nvar, dde.nvar.data(),
                         Dde::f_ddes(), &n_beta, dde.beta_cc.data(), &n_his, dde.history_cc.data(),
                         &n_tspan, dde.tspan.data(),
                         // output
                         &sol_npts, &sol_flag, &sol_ne,
                         &sol_t_ptr, &sol_y_ptr, &sol_te_ptr, &sol_ye_ptr,
                         &sol_queue_ptr, &sol_yoft_ptr, &sol_tqueue_ptr,
                         &sol_stats_ptr, &sol_ie_ptr, &sol_ipoint_ptr,
                         &sol_shift, &sol_tshift,
                         // options
                         &nre, dde.rerr.data(),
                         &nae, dde.aerr.data(),
                         &n_jumps, dde.jumps.data(),
                         &n_thit, dde.thit_exactly.data(),
                         &n_direction, dde.direction.data(),
                         &n_isterminal, dde.isterminal.data(),
                         &dde.opts,
                         Dde::f_event(), Dde::f_change(),
                         Dde::f_out(), Dde::f_user());

      DdeSol sol(dde.nvar[0], dde.opts.interpolation,
                 sol_npts, sol_flag, sol_ne,
                 sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,
                 sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,
                 sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,
                 sol_shift, sol_tshift);

      return sol;
    }
  };
}
#endif
