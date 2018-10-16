#ifndef DDE_SOLVER_H
#define DDE_SOLVER_H

// fuction pointers for C interface
typedef void (*f_ddes_fcn_c) (const double*, const int*,
                              const int*, const double[], const double[], double[]);
typedef void (*f_beta_fcn_c) (const double*, const int*, const int*, const double[],
                              double[]);
typedef void (*f_history_fcn_c) (const double*, const int*, double[]);

typedef void (*f_event_fcn_c) (const double*, const int*,
                               const int*, const double[], const double[],
                               const double[], double[]);
typedef void (*f_change_fcn_c) (const int*, const double*, 
                                double[], double[], double*, const int*, int[],
                                const int*, bool[], bool*);
typedef void (*f_out_fcn_c) (double* t, double[], double[], int*, int*);
typedef void (*f_user_trim_get_c) ();

void integrate_dde_1(int* n_nvar, int nvar[],
                     f_ddes_fcn_c,
                     f_beta_fcn_c,
                     f_history_fcn_c,
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
                     f_event_fcn_c,
                     f_change_fcn_c,
                     f_out_fcn_c,
                     f_user_trim_get_c);

void integrate_dde_2(int* n_nvar, int nvar[],
                     f_ddes_fcn_c,
                     int* n_delay, double delay[],
                     f_history_fcn_c,
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
                     f_event_fcn_c,
                     f_change_fcn_c,
                     f_out_fcn_c,
                     f_user_trim_get_c);

void integrate_dde_3(int* n_nvar, int nvar[],
                     f_ddes_fcn_c,
                     f_beta_fcn_c,
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
                     f_event_fcn_c,
                     f_change_fcn_c,
                     f_out_fcn_c,
                     f_user_trim_get_c);

void integrate_dde_4(int* n_nvar, int nvar[],
                     f_ddes_fcn_c,
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
                     f_event_fcn_c,
                     f_change_fcn_c,
                     f_out_fcn_c,
                     f_user_trim_get_c);

#endif
