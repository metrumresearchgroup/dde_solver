#ifndef DDE_SOLVER_HPP
#define DDE_SOLVER_HPP

#include <vector>

extern "C" {

  void integrate_dde_1(int* n_nvar, int nvar[],
                       void (*f_ddes_cc) (double*, int*, int*, double[], double[], double[]),
                       void (*f_beta_cc) (double*, int*, int*, double[], double[]),
                       void (*f_history_cc) (double*, int*, double[]),
                       int* n_tspan, double tspan[],
                       int* sol_npts,
                       double** sol_t_ptr,
                       double** sol_y_ptr,
                       double** sol_te_ptr,
                       double** sol_ye_ptr,
                       int* n_re_vector, double re_vector[],
                       int* n_ae_vector, double ae_vector[],
                       int* n_jumps, double jumps[],
                       int* n_thit, double thit_exactly[],
                       int* n_direction, int direction[],
                       int* n_isterminal, bool isterminal[],
                       struct dde_opts_cc* opts_cc,
                       void (*f_event_fcn_cc) (double*, int*, int*, double[], double[], double[], double[]),
                       void (*f_change_fcn_cc) (int*, int*, double*, double[], double[], double*, int*, int[], int*, bool[], bool*),
                       void (*f_out_fcn_cc) (double* t, double[], double[], int*, int*),
                       void (*f_user_trim_get_cc) ());
}

#endif
