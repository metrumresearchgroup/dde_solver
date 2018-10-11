#ifndef DDE_SOLVER_HPP
#define DDE_SOLVER_HPP

#include <vector>

extern "C" {

  void integrate_dde_1(int* n_nvar, int nvar[],
                       void (*f_ddes_cc) (double*, int*, int*, double[], double[], double[]),
                       void (*f_beta_cc) (double*, int*, int*, double[], double[]),
                       void (*f_history_cc) (double*, int*, double[]),
                       int* n_tspan, double tspan[],
                       // output
                       int* sol_npts, int* sol_flag, int* sol_ne,
                       double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                       double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                       int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                       bool* sol_shift, double* sol_tshift,
                       // options
                       int* n_re_vector, double re_vector[],
                       int* n_ae_vector, double ae_vector[],
                       int* n_jumps, double jumps[],
                       int* n_thit, double thit_exactly[],
                       int* n_direction, int direction[],
                       int* n_isterminal, const bool isterminal[],
                       struct dde_opts_cc* opts_cc,
                       void (*f_event_fcn_cc) (double*, int*, int*, double[], double[], double[], double[]),
                       void (*f_change_fcn_cc) (int*, double*, double[], double[], double*, int*, int[], int*, bool[], bool*),
                       void (*f_out_fcn_cc) (double* t, double[], double[], int*, int*),
                       void (*f_user_trim_get_cc) ());

  void integrate_dde_2(int* n_nvar, int nvar[],
                       void (*f_ddes_cc) (double*, int*, int*, double[], double[], double[]),
                       int* n_delay, double delay[],
                       void (*f_history_cc) (double*, int*, double[]),
                       int* n_tspan, double tspan[],
                       // output
                       int* sol_npts, int* sol_flag, int* sol_ne,
                       double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                       double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                       int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                       bool* sol_shift, double* sol_tshift,
                       // options
                       int* n_re_vector, double re_vector[],
                       int* n_ae_vector, double ae_vector[],
                       int* n_jumps, double jumps[],
                       int* n_thit, double thit_exactly[],
                       int* n_direction, int direction[],
                       int* n_isterminal, const bool isterminal[],
                       struct dde_opts_cc* opts_cc,
                       void (*f_event_fcn_cc) (double*, int*, int*, double[], double[], double[], double[]),
                       void (*f_change_fcn_cc) (int*, double*, double[], double[], double*, int*, int[], int*, bool[], bool*),
                       void (*f_out_fcn_cc) (double* t, double[], double[], int*, int*),
                       void (*f_user_trim_get_cc) ());

  void integrate_dde_3(int* n_nvar, int nvar[],
                       void (*f_ddes_cc) (double*, int*, int*, double[], double[], double[]),
                       void (*f_beta_cc) (double*, int*, int*, double[], double[]),
                       int* n_history, double history[],
                       int* n_tspan, double tspan[],
                       // output
                       int* sol_npts, int* sol_flag, int* sol_ne,
                       double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                       double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                       int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                       bool* sol_shift, double* sol_tshift,
                       // options
                       int* n_re_vector, double re_vector[],
                       int* n_ae_vector, double ae_vector[],
                       int* n_jumps, double jumps[],
                       int* n_thit, double thit_exactly[],
                       int* n_direction, int direction[],
                       int* n_isterminal, const bool isterminal[],
                       struct dde_opts_cc* opts_cc,
                       void (*f_event_fcn_cc) (double*, int*, int*, double[], double[], double[], double[]),
                       void (*f_change_fcn_cc) (int*, double*, double[], double[], double*, int*, int[], int*, bool[], bool*),
                       void (*f_out_fcn_cc) (double* t, double[], double[], int*, int*),
                       void (*f_user_trim_get_cc) ());

  void integrate_dde_4(int* n_nvar, int nvar[],
                       void (*f_ddes_cc) (double*, int*, int*, double[], double[], double[]),
                       int* n_delay, double delay[],
                       int* n_history, double history[],
                       int* n_tspan, double tspan[],
                       // output
                       int* sol_npts, int* sol_flag, int* sol_ne,
                       double** sol_t_ptr, double** sol_y_ptr, double** sol_te_ptr, double** sol_ye_ptr,
                       double** sol_queue_ptr, double** sol_yoft_ptr, double** sol_tqueue_ptr,
                       int** sol_stats_ptr, int** sol_ie_ptr, int** sol_ipoint_ptr,
                       bool* sol_shift, double* sol_tshift,
                       // options
                       int* n_re_vector, double re_vector[],
                       int* n_ae_vector, double ae_vector[],
                       int* n_jumps, double jumps[],
                       int* n_thit, double thit_exactly[],
                       int* n_direction, int direction[],
                       int* n_isterminal, const bool isterminal[],
                       struct dde_opts_cc* opts_cc,
                       void (*f_event_fcn_cc) (double*, int*, int*, double[], double[], double[], double[]),
                       void (*f_change_fcn_cc) (int*, double*, double[], double[], double*, int*, int[], int*, bool[], bool*),
                       void (*f_out_fcn_cc) (double* t, double[], double[], int*, int*),
                       void (*f_user_trim_get_cc) ());
}

#endif
