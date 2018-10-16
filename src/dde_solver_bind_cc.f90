module dde_solver_bind
  use iso_c_binding

  implicit none

  integer, parameter:: dp=kind(0.d0)

  ! C++ interface callback functions
  interface
     subroutine f_ddes_cc(t, n, nlags, y, z, dy, user_data) bind (c)
       use iso_c_binding
       double precision :: t
       integer :: n, nlags
       double precision, dimension(n) :: y, dy
       double precision, dimension(n,nlags) :: z
       type(c_ptr), intent(inout) :: user_data;
       intent(in):: t,n,nlags,y,z
       intent(out) :: dy
     end subroutine f_ddes_cc
  end interface
  interface
     subroutine f_history_cc(t,n,y, user_data) bind (c)
       use iso_c_binding
       double precision :: t
       integer :: n
       double precision, dimension(n) :: y
       type(c_ptr), intent(inout) :: user_data;
       intent(in):: t,n
       intent(out) :: y
     end subroutine f_history_cc
  end interface
  interface
     subroutine f_beta_cc(t, n, nlags, y, bval, user_data) bind (c)
       use iso_c_binding
       double precision :: t
       integer, intent(in) :: n, nlags
       double precision, dimension(n) :: y
       double precision, dimension(nlags) :: bval
       type(c_ptr), intent(inout) :: user_data;
       intent(in):: t,y
       intent(out) :: bval
     end subroutine f_beta_cc
  end interface
  interface
     subroutine f_change_fcn_cc(nevent,tevent,yevent,dyevent,hinit, &
          n_direction, direction,&
          n_isterminal, isterminal, quit, user_data) bind (c)
       use iso_c_binding
       integer :: nevent, n_direction, n_isterminal
       integer, dimension(n_direction) :: direction
       double precision :: tevent,hinit
       double precision, dimension(nevent) :: yevent,dyevent
       logical :: quit
       logical, dimension(n_isterminal) :: isterminal
       type(c_ptr), intent(inout) :: user_data;
       intent(in) :: nevent,tevent,n_direction,n_isterminal
       intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
     end subroutine f_change_fcn_cc
  end interface
  interface
     subroutine f_event_fcn_cc(t,n, nlags, y,dydt,z,g,user_data) bind (c)
       use iso_c_binding
       integer, intent(in) :: n, nlags
       double precision :: t
       double precision, dimension(n) :: y,dydt
       double precision, dimension(n,nlags) :: z
       double precision, dimension(n) :: g
       type(c_ptr), intent(inout) :: user_data;
       intent(in):: t,y,dydt,z
       intent(out) :: g
     end subroutine f_event_fcn_cc
  end interface
  interface
     subroutine f_out_fcn_cc(t,y,dy,n,nevent,user_data) bind (c)
       use iso_c_binding
       integer :: n,nevent
       double precision :: t
       double precision, dimension(n) :: y, dy
       type(c_ptr), intent(inout) :: user_data;
     end subroutine f_out_fcn_cc
  end interface
  interface
     subroutine f_user_trim_get_cc(user_data) bind (c)
       use iso_c_binding
       type(c_ptr), intent(inout) :: user_data;
     end subroutine f_user_trim_get_cc
  end interface

contains

  ! C++ binding for DKL_1
  subroutine integrate_dde_1_cc(user_data,&
       n_nvar, nvar, ddes_cc, beta_cc, history_cc,&
       n_tspan, tspan, &
       ! output
       sol_npts, sol_flag, sol_ne, &
       sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,&
       sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,&
       sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,&
       sol_shift, sol_tshift, &
       ! input parameters
       n_re_vector, re_vector, &
       n_ae_vector, ae_vector, &
       n_jumps, jumps, &
       n_thit, thit_exactly, &
       n_direction, direction, &
       n_isterminal, isterminal,&
       opts_cc, &
       event_fcn_cc, change_fcn_cc, out_fcn_cc, user_trim_get_cc) bind(c)

    use dde_solver_m
    implicit none

    type(c_ptr), intent(inout) :: user_data

    integer(c_int), intent(in) :: n_re_vector, n_ae_vector, n_jumps, n_thit
    real(c_double), intent(in) :: re_vector(n_re_vector)
    real(c_double), intent(in) :: ae_vector(n_ae_vector)
    real(c_double), intent(in) :: jumps(n_jumps), thit_exactly(n_thit)
    integer, intent(in) :: n_direction
    integer(c_int), intent(in) :: direction(n_direction)
    integer, intent(in) :: n_isterminal
    logical(c_bool), intent(in) :: isterminal(n_isterminal)
    type (dde_opts_cc), intent(in) :: opts_cc

    integer, intent(in) :: n_tspan, n_nvar
    real(c_double) :: tspan(n_tspan)
    integer :: nvar(n_nvar)

    procedure(f_ddes_cc)          :: ddes_cc
    procedure(f_beta_cc)          :: beta_cc
    procedure(f_history_cc)       :: history_cc
    procedure(f_event_fcn_cc)    , optional :: event_fcn_cc
    procedure(f_change_fcn_cc)   , optional :: change_fcn_cc
    procedure(f_out_fcn_cc)      , optional :: out_fcn_cc
    procedure(f_user_trim_get_cc), optional :: user_trim_get_cc

    ! output
    integer(c_int), intent(out) :: sol_npts, sol_flag, sol_ne
    type(c_ptr) :: sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr
    type(c_ptr) :: sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr
    type(c_ptr) :: sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr
    logical(c_bool), intent(out) :: sol_shift
    real(c_double), intent(out) :: sol_tshift

    ! local
    type (dde_opts) :: options
    type (dde_sol) :: sol
    logical :: have_event
    logical :: have_change
    logical :: have_out
    logical :: have_trim
    logical :: isterminal_ff(n_isterminal)

    integer :: i

    have_event = present(event_fcn_cc)
    have_change = present(change_fcn_cc)
    have_out = present(out_fcn_cc)
    have_trim = present(user_trim_get_cc)

    ! fortran's logical is logical(4), while logical(c_bool)
    ! is logical(1)
    isterminal_ff = isterminal

    options = dde_set_cc(n_re_vector, re_vector, &
         n_ae_vector, ae_vector, &
         n_jumps, jumps, &
         n_thit, thit_exactly, &
         n_direction, direction, &
         n_isterminal, isterminal_ff, &
         opts_cc)

    if(have_event) then
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c)
             endif
          endif
       endif
    else
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_1(nvar,ddes_c,beta_c,history_c,tspan,options)
             endif
          endif
       endif
    endif

    ! copy to C output
    sol_npts       = sol%npts
    sol_flag       = sol%flag
    sol_ne         = sol%ne
    sol_t_ptr      = c_loc(sol%t)
    sol_y_ptr      = c_loc(sol%y)
    sol_te_ptr     = c_loc(sol%te)
    sol_ye_ptr     = c_loc(sol%ye)
    sol_queue_ptr  = c_loc(sol%queue)
    sol_yoft_ptr   = c_loc(sol%yoft)
    sol_tqueue_ptr = c_loc(sol%tqueue)
    sol_stats_ptr  = c_loc(sol%stats)
    sol_ie_ptr     = c_loc(sol%ie)
    sol_ipoint_ptr = c_loc(sol%ipoint)
    sol_shift      = sol%shift
    sol_tshift     = sol%tshift

    ! call print_stats(sol)
    ! call release_arrays(sol, options)

  contains
    subroutine ddes_c(t, y, z, dy) bind (c)
       double precision :: t
       double precision, dimension(:) :: y, dy
       double precision, dimension(:,:) :: z
       intent(in):: t,y,z
       intent(out) :: dy
      call ddes_cc(t, size(y), size(z, 2), y, z, dy, user_data)
    end subroutine ddes_c

    subroutine history_c(t,y) bind (c)
      double precision :: t
      double precision, dimension(:) :: y
      intent(in):: t
      intent(out) :: y
      call history_cc(t, size(y), y, user_data)
    end subroutine history_c

    subroutine beta_c(t, y, bval) bind (c)
      double precision :: t
      double precision, dimension(:) :: y
      double precision, dimension(:) :: bval
      intent(in):: t,y
      intent(out) :: bval
      call beta_cc(t, size(y), size(bval), y, bval, user_data)
    end subroutine beta_c

    subroutine change_fcn_c(nevent,tevent,yevent,dyevent,hinit, &
         direction, isterminal, quit) bind (c)
      integer, intent(in) :: nevent
      integer, dimension(:) :: direction
      double precision :: tevent,hinit
      double precision, dimension(:) :: yevent,dyevent
      logical :: quit
      logical, dimension(:) :: isterminal
      intent(in) :: tevent
      intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
      call change_fcn_cc(nevent,tevent,yevent,dyevent,hinit, &
         size(direction), direction,&
         size(isterminal), isterminal, quit, user_data)
    end subroutine change_fcn_c

    subroutine event_fcn_c(t,y,dydt,z,g) bind (c)
      double precision :: t
      double precision, dimension(:) :: y,dydt
      double precision, dimension(:,:) :: z
      double precision, dimension(:) :: g
      intent(in):: t,y,dydt,z
      intent(out) :: g
      call event_fcn_cc(t, size(y), size(z, 2), y, dydt, z, g, user_data)
    end subroutine event_fcn_c

    subroutine out_fcn_c(t,y,dy,n,nevent) bind (c)
      integer :: n,nevent
      double precision :: t
      double precision, dimension(:) :: y, dy
      call out_fcn_cc(t,y,dy,n,nevent,user_data)
    end subroutine out_fcn_c

    subroutine user_trim_get_c() bind (c)
      call user_trim_get_cc(user_data)
    end subroutine user_trim_get_c

  end subroutine integrate_dde_1_cc

  ! C++ binding for DKL_2
  subroutine integrate_dde_2_cc(user_data,&
       n_nvar, nvar, ddes_cc, n_delay, delay, history_cc,&
       n_tspan, tspan, &
       ! output
       sol_npts, sol_flag, sol_ne, &
       sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,&
       sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,&
       sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,&
       sol_shift, sol_tshift, &
       ! input parameters
       n_re_vector, re_vector, &
       n_ae_vector, ae_vector, &
       n_jumps, jumps, &
       n_thit, thit_exactly, &
       n_direction, direction, &
       n_isterminal, isterminal,&
       opts_cc, &
       event_fcn_cc, change_fcn_cc, out_fcn_cc, user_trim_get_cc) bind(c)

    use dde_solver_m
    implicit none

    type(c_ptr), intent(inout) :: user_data

    integer(c_int), intent(in) :: n_re_vector, n_ae_vector, n_jumps, n_thit
    real(c_double), intent(in) :: re_vector(n_re_vector)
    real(c_double), intent(in) :: ae_vector(n_ae_vector)
    real(c_double), intent(in) :: jumps(n_jumps), thit_exactly(n_thit)
    integer, intent(in) :: n_direction
    integer(c_int), intent(in) :: direction(n_direction)
    integer, intent(in) :: n_isterminal
    logical(c_bool), intent(in) :: isterminal(n_isterminal)
    type (dde_opts_cc), intent(in) :: opts_cc

    integer, intent(in) :: n_tspan, n_nvar
    real(c_double) :: tspan(n_tspan)
    integer :: nvar(n_nvar)

    procedure(f_ddes_cc)          :: ddes_cc
    integer(c_int), intent(in)    :: n_delay
    real(c_double)                :: delay(n_delay)
    procedure(f_history_cc)       :: history_cc
    procedure(f_event_fcn_cc)    , optional :: event_fcn_cc
    procedure(f_change_fcn_cc)   , optional :: change_fcn_cc
    procedure(f_out_fcn_cc)      , optional :: out_fcn_cc
    procedure(f_user_trim_get_cc), optional :: user_trim_get_cc

    ! output
    integer(c_int), intent(out) :: sol_npts, sol_flag, sol_ne
    type(c_ptr) :: sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr
    type(c_ptr) :: sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr
    type(c_ptr) :: sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr
    logical(c_bool), intent(out) :: sol_shift
    real(c_double), intent(out) :: sol_tshift

    ! local
    type (dde_opts) :: options
    type (dde_sol) :: sol
    logical :: have_event
    logical :: have_change
    logical :: have_out
    logical :: have_trim
    logical :: isterminal_ff(n_isterminal)

    integer :: i

    have_event = present(event_fcn_cc)
    have_change = present(change_fcn_cc)
    have_out = present(out_fcn_cc)
    have_trim = present(user_trim_get_cc)

    ! fortran's logical is logical(4), while logical(c_bool)
    ! is logical(1)
    isterminal_ff = isterminal

    options = dde_set_cc(n_re_vector, re_vector, &
         n_ae_vector, ae_vector, &
         n_jumps, jumps, &
         n_thit, thit_exactly, &
         n_direction, direction, &
         n_isterminal, isterminal_ff, &
         opts_cc)

    if(have_event) then
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     EVENT_FCN = event_fcn_c)
             endif
          endif
       endif
    else
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_2(nvar,ddes_c,delay,history_c,tspan,options)
             endif
          endif
       endif
    endif

    ! copy to C output
    sol_npts       = sol%npts
    sol_flag       = sol%flag
    sol_ne         = sol%ne
    sol_t_ptr      = c_loc(sol%t)
    sol_y_ptr      = c_loc(sol%y)
    sol_te_ptr     = c_loc(sol%te)
    sol_ye_ptr     = c_loc(sol%ye)
    sol_queue_ptr  = c_loc(sol%queue)
    sol_yoft_ptr   = c_loc(sol%yoft)
    sol_tqueue_ptr = c_loc(sol%tqueue)
    sol_stats_ptr  = c_loc(sol%stats)
    sol_ie_ptr     = c_loc(sol%ie)
    sol_ipoint_ptr = c_loc(sol%ipoint)
    sol_shift      = sol%shift
    sol_tshift     = sol%tshift

    ! call print_stats(sol)
    ! call release_arrays(sol, options)

  contains
    subroutine ddes_c(t, y, z, dy) bind (c)
       double precision :: t
       double precision, dimension(:) :: y, dy
       double precision, dimension(:,:) :: z
       intent(in):: t,y,z
       intent(out) :: dy
      call ddes_cc(t, size(y), size(z, 2), y, z, dy, user_data)
    end subroutine ddes_c

    subroutine history_c(t,y) bind (c)
      double precision :: t
      double precision, dimension(:) :: y
      intent(in):: t
      intent(out) :: y
      call history_cc(t, size(y), y, user_data)
    end subroutine history_c

    subroutine change_fcn_c(nevent,tevent,yevent,dyevent,hinit, &
         direction, isterminal, quit) bind (c)
      integer, intent(in) :: nevent
      integer, dimension(:) :: direction
      double precision :: tevent,hinit
      double precision, dimension(:) :: yevent,dyevent
      logical :: quit
      logical, dimension(:) :: isterminal
      intent(in) :: tevent
      intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
      call change_fcn_cc(nevent,tevent,yevent,dyevent,hinit, &
         size(direction), direction,&
         size(isterminal), isterminal, quit, user_data)
    end subroutine change_fcn_c

    subroutine event_fcn_c(t,y,dydt,z,g) bind (c)
      double precision :: t
      double precision, dimension(:) :: y,dydt
      double precision, dimension(:,:) :: z
      double precision, dimension(:) :: g
      intent(in):: t,y,dydt,z
      intent(out) :: g
      call event_fcn_cc(t, size(y), size(z, 2), y, dydt, z, g, user_data)
    end subroutine event_fcn_c

    subroutine out_fcn_c(t,y,dy,n,nevent) bind (c)
      integer :: n,nevent
      double precision :: t
      double precision, dimension(:) :: y, dy
      call out_fcn_cc(t,y,dy,n,nevent,user_data)
    end subroutine out_fcn_c

    subroutine user_trim_get_c() bind (c)
      call user_trim_get_cc(user_data)
    end subroutine user_trim_get_c

  end subroutine integrate_dde_2_cc

  ! C++ binding for DKL_3
  subroutine integrate_dde_3_cc(user_data,&
       n_nvar, nvar, ddes_cc, beta_cc, n_his, history,&
       n_tspan, tspan, &
       ! output
       sol_npts, sol_flag, sol_ne, &
       sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,&
       sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,&
       sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,&
       sol_shift, sol_tshift, &
       ! input parameters
       n_re_vector, re_vector, &
       n_ae_vector, ae_vector, &
       n_jumps, jumps, &
       n_thit, thit_exactly, &
       n_direction, direction, &
       n_isterminal, isterminal,&
       opts_cc, &
       event_fcn_cc, change_fcn_cc, out_fcn_cc, user_trim_get_cc) bind(c)

    use dde_solver_m
    implicit none

    type(c_ptr), intent(inout) :: user_data

    integer(c_int), intent(in) :: n_re_vector, n_ae_vector, n_jumps, n_thit
    real(c_double), intent(in) :: re_vector(n_re_vector)
    real(c_double), intent(in) :: ae_vector(n_ae_vector)
    real(c_double), intent(in) :: jumps(n_jumps), thit_exactly(n_thit)
    integer, intent(in) :: n_direction
    integer(c_int), intent(in) :: direction(n_direction)
    integer, intent(in) :: n_isterminal
    logical(c_bool), intent(in) :: isterminal(n_isterminal)
    type (dde_opts_cc), intent(in) :: opts_cc

    integer, intent(in) :: n_tspan, n_nvar
    real(c_double) :: tspan(n_tspan)
    integer :: nvar(n_nvar)

    procedure(f_ddes_cc)          :: ddes_cc
    procedure(f_beta_cc)          :: beta_cc
    integer(c_int)                :: n_his
    real(c_double)                :: history(n_his)
    procedure(f_event_fcn_cc)    , optional :: event_fcn_cc
    procedure(f_change_fcn_cc)   , optional :: change_fcn_cc
    procedure(f_out_fcn_cc)      , optional :: out_fcn_cc
    procedure(f_user_trim_get_cc), optional :: user_trim_get_cc

    ! output
    integer(c_int), intent(out) :: sol_npts, sol_flag, sol_ne
    type(c_ptr) :: sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr
    type(c_ptr) :: sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr
    type(c_ptr) :: sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr
    logical(c_bool), intent(out) :: sol_shift
    real(c_double), intent(out) :: sol_tshift

    ! local
    type (dde_opts) :: options
    type (dde_sol) :: sol
    logical :: have_event
    logical :: have_change
    logical :: have_out
    logical :: have_trim
    logical :: isterminal_ff(n_isterminal)

    integer :: i

    have_event = present(event_fcn_cc)
    have_change = present(change_fcn_cc)
    have_out = present(out_fcn_cc)
    have_trim = present(user_trim_get_cc)

    ! fortran's logical is logical(4), while logical(c_bool)
    ! is logical(1)
    isterminal_ff = isterminal

    options = dde_set_cc(n_re_vector, re_vector, &
         n_ae_vector, ae_vector, &
         n_jumps, jumps, &
         n_thit, thit_exactly, &
         n_direction, direction, &
         n_isterminal, isterminal_ff, &
         opts_cc)

    if(have_event) then
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     EVENT_FCN = event_fcn_c)
             endif
          endif
       endif
    else
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_3(nvar,ddes_c,beta_c,history,tspan,options)
             endif
          endif
       endif
    endif

    ! copy to C output
    sol_npts       = sol%npts
    sol_flag       = sol%flag
    sol_ne         = sol%ne
    sol_t_ptr      = c_loc(sol%t)
    sol_y_ptr      = c_loc(sol%y)
    sol_te_ptr     = c_loc(sol%te)
    sol_ye_ptr     = c_loc(sol%ye)
    sol_queue_ptr  = c_loc(sol%queue)
    sol_yoft_ptr   = c_loc(sol%yoft)
    sol_tqueue_ptr = c_loc(sol%tqueue)
    sol_stats_ptr  = c_loc(sol%stats)
    sol_ie_ptr     = c_loc(sol%ie)
    sol_ipoint_ptr = c_loc(sol%ipoint)
    sol_shift      = sol%shift
    sol_tshift     = sol%tshift

    ! call print_stats(sol)
    ! call release_arrays(sol, options)

  contains
    subroutine ddes_c(t, y, z, dy) bind (c)
       double precision :: t
       double precision, dimension(:) :: y, dy
       double precision, dimension(:,:) :: z
       intent(in):: t,y,z
       intent(out) :: dy
      call ddes_cc(t, size(y), size(z, 2), y, z, dy, user_data)
    end subroutine ddes_c

    subroutine beta_c(t, y, bval) bind (c)
      double precision :: t
      double precision, dimension(:) :: y
      double precision, dimension(:) :: bval
      intent(in):: t,y
      intent(out) :: bval
      call beta_cc(t, size(y), size(bval), y, bval, user_data)
    end subroutine beta_c

    subroutine change_fcn_c(nevent,tevent,yevent,dyevent,hinit, &
         direction, isterminal, quit) bind (c)
      integer, intent(in) :: nevent
      integer, dimension(:) :: direction
      double precision :: tevent,hinit
      double precision, dimension(:) :: yevent,dyevent
      logical :: quit
      logical, dimension(:) :: isterminal
      intent(in) :: tevent
      intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
      call change_fcn_cc(nevent,tevent,yevent,dyevent,hinit, &
         size(direction), direction,&
         size(isterminal), isterminal, quit, user_data)
    end subroutine change_fcn_c

    subroutine event_fcn_c(t,y,dydt,z,g) bind (c)
      double precision :: t
      double precision, dimension(:) :: y,dydt
      double precision, dimension(:,:) :: z
      double precision, dimension(:) :: g
      intent(in):: t,y,dydt,z
      intent(out) :: g
      call event_fcn_cc(t, size(y), size(z, 2), y, dydt, z, g, user_data)
    end subroutine event_fcn_c

    subroutine out_fcn_c(t,y,dy,n,nevent) bind (c)
      integer :: n,nevent
      double precision :: t
      double precision, dimension(:) :: y, dy
      call out_fcn_cc(t,y,dy,n,nevent,user_data)
    end subroutine out_fcn_c

    subroutine user_trim_get_c() bind (c)
      call user_trim_get_cc(user_data)
    end subroutine user_trim_get_c

  end subroutine integrate_dde_3_cc

  ! C++ binding for DKL_4
  subroutine integrate_dde_4_cc(user_data,&
       n_nvar, nvar, ddes_cc, n_delay, delay, n_his, history,&
       n_tspan, tspan, &
       ! output
       sol_npts, sol_flag, sol_ne, &
       sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,&
       sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr,&
       sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr,&
       sol_shift, sol_tshift, &
       ! input parameters
       n_re_vector, re_vector, &
       n_ae_vector, ae_vector, &
       n_jumps, jumps, &
       n_thit, thit_exactly, &
       n_direction, direction, &
       n_isterminal, isterminal,&
       opts_cc, &
       event_fcn_cc, change_fcn_cc, out_fcn_cc, user_trim_get_cc) bind(c)

    use dde_solver_m
    implicit none

    type(c_ptr), intent(inout) :: user_data

    integer(c_int), intent(in) :: n_re_vector, n_ae_vector, n_jumps, n_thit
    real(c_double), intent(in) :: re_vector(n_re_vector)
    real(c_double), intent(in) :: ae_vector(n_ae_vector)
    real(c_double), intent(in) :: jumps(n_jumps), thit_exactly(n_thit)
    integer, intent(in) :: n_direction
    integer(c_int), intent(in) :: direction(n_direction)
    integer, intent(in) :: n_isterminal
    logical(c_bool), intent(in) :: isterminal(n_isterminal)
    type (dde_opts_cc), intent(in) :: opts_cc

    integer, intent(in) :: n_tspan, n_nvar
    real(c_double) :: tspan(n_tspan)
    integer :: nvar(n_nvar)

    procedure(f_ddes_cc)          :: ddes_cc
    integer(c_int)                :: n_delay, n_his
    real(c_double)                :: delay(n_delay), history(n_his)
    procedure(f_event_fcn_cc)    , optional :: event_fcn_cc
    procedure(f_change_fcn_cc)   , optional :: change_fcn_cc
    procedure(f_out_fcn_cc)      , optional :: out_fcn_cc
    procedure(f_user_trim_get_cc), optional :: user_trim_get_cc

    ! output
    integer(c_int), intent(out) :: sol_npts, sol_flag, sol_ne
    type(c_ptr) :: sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr
    type(c_ptr) :: sol_queue_ptr, sol_yoft_ptr, sol_tqueue_ptr
    type(c_ptr) :: sol_stats_ptr, sol_ie_ptr, sol_ipoint_ptr
    logical(c_bool), intent(out) :: sol_shift
    real(c_double), intent(out) :: sol_tshift

    ! local
    type (dde_opts) :: options
    type (dde_sol) :: sol
    logical :: have_event
    logical :: have_change
    logical :: have_out
    logical :: have_trim
    logical :: isterminal_ff(n_isterminal)

    integer :: i

    have_event = present(event_fcn_cc)
    have_change = present(change_fcn_cc)
    have_out = present(out_fcn_cc)
    have_trim = present(user_trim_get_cc)

    ! fortran's logical is logical(4), while logical(c_bool)
    ! is logical(1)
    isterminal_ff = isterminal

    options = dde_set_cc(n_re_vector, re_vector, &
         n_ae_vector, ae_vector, &
         n_jumps, jumps, &
         n_thit, thit_exactly, &
         n_direction, direction, &
         n_isterminal, isterminal_ff, &
         opts_cc)

    if(have_event) then
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     EVENT_FCN = event_fcn_c)
             endif
          endif
       endif
    else
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     CHANGE_FCN = change_fcn_c)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     OUT_FCN = out_fcn_c, &
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     OUT_FCN = out_fcn_c)
             endif
          else
             if(have_trim) then
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options,&
                     USER_TRIM_GET = user_trim_get_c)
             else
                sol = dkl_4(nvar,ddes_c,delay,history,tspan,options)
             endif
          endif
       endif
    endif

    ! copy to C output
    sol_npts       = sol%npts
    sol_flag       = sol%flag
    sol_ne         = sol%ne
    sol_t_ptr      = c_loc(sol%t)
    sol_y_ptr      = c_loc(sol%y)
    sol_te_ptr     = c_loc(sol%te)
    sol_ye_ptr     = c_loc(sol%ye)
    sol_queue_ptr  = c_loc(sol%queue)
    sol_yoft_ptr   = c_loc(sol%yoft)
    sol_tqueue_ptr = c_loc(sol%tqueue)
    sol_stats_ptr  = c_loc(sol%stats)
    sol_ie_ptr     = c_loc(sol%ie)
    sol_ipoint_ptr = c_loc(sol%ipoint)
    sol_shift      = sol%shift
    sol_tshift     = sol%tshift

    ! call print_stats(sol)
    ! call release_arrays(sol, options)

  contains
    subroutine ddes_c(t, y, z, dy) bind (c)
       double precision :: t
       double precision, dimension(:) :: y, dy
       double precision, dimension(:,:) :: z
       intent(in):: t,y,z
       intent(out) :: dy
      call ddes_cc(t, size(y), size(z, 2), y, z, dy, user_data)
    end subroutine ddes_c

    subroutine change_fcn_c(nevent,tevent,yevent,dyevent,hinit, &
         direction, isterminal, quit) bind (c)
      integer, intent(in) :: nevent
      integer, dimension(:) :: direction
      double precision :: tevent,hinit
      double precision, dimension(:) :: yevent,dyevent
      logical :: quit
      logical, dimension(:) :: isterminal
      intent(in) :: tevent
      intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
      call change_fcn_cc(nevent,tevent,yevent,dyevent,hinit, &
         size(direction), direction,&
         size(isterminal), isterminal, quit, user_data)
    end subroutine change_fcn_c

    subroutine event_fcn_c(t,y,dydt,z,g) bind (c)
      double precision :: t
      double precision, dimension(:) :: y,dydt
      double precision, dimension(:,:) :: z
      double precision, dimension(:) :: g
      intent(in):: t,y,dydt,z
      intent(out) :: g
      call event_fcn_cc(t, size(y), size(z, 2), y, dydt, z, g, user_data)
    end subroutine event_fcn_c

    subroutine out_fcn_c(t,y,dy,n,nevent) bind (c)
      integer :: n,nevent
      double precision :: t
      double precision, dimension(:) :: y, dy
      call out_fcn_cc(t,y,dy,n,nevent,user_data)
    end subroutine out_fcn_c

    subroutine user_trim_get_c() bind (c)
      call user_trim_get_cc(user_data)
    end subroutine user_trim_get_c

  end subroutine integrate_dde_4_cc

end module dde_solver_bind
