#ifndef DDE_SOLVER_HPP
#define DDE_SOLVER_HPP

#include <vector>

extern "C" {

  void integrate_dde_1(double* re, double[] re_vector, int* n_re_vector, 
                       double* ae, double[] ae_vector, int* n_ae_vector, 
                       double* hinit, double* hmax,
                       int* n_isterminal, double[] isterminal, 
                       int* n_jumps, int* n_thit, int* n_direction, 
                       int[] direction, double[] jumps,
                       bool* neutral, bool* track_discontinuities,
                       int* tracking_level, 
                       bool* interpolation, double* thit_exactly,
                       int* max_events, int* max_steps, 
                       int* moving_average, int* max_delay, int* trim_frequency,
                       int* n_nvar, int[] nvar,
                       void (*ddes) (),
                       void (*beta) (),
                       void (*history) (),
                       int* n_tspan, double[] tspan,
                       void (*event_fcn) (), 
                       void (*change_fcn) (),
                       void (*out_fcn) (),
                       void (*user_trim_get) (), 
                       int* sol_npts,
                       double** sol_t_ptr,
                       double** sol_y_ptr,
                       double** sol_te_ptr,
                       double** sol_ye_ptr);

}

#endif
