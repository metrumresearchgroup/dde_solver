module dde_solver_bind
  use iso_c_binding

  implicit none

  integer, parameter:: dp=kind(0.d0)

  interface
     subroutine f_ddes_cc(t, n, nlags, y, z, dy) bind (c)
       double precision :: t
       integer :: n, nlags
       double precision, dimension(n) :: y, dy
       double precision, dimension(n,nlags) :: z
       intent(in):: t,y,z
       intent(out) :: dy
     end subroutine f_ddes_cc
  end interface
  interface
     subroutine f_history_cc(t,n,y) bind (c)
       double precision :: t
       integer :: n
       double precision, dimension(n) :: y
       intent(in):: t
       intent(out) :: y
     end subroutine f_history_cc
  end interface
  interface
     subroutine f_beta_cc(t, n, nlags, y, bval) bind (c)
       double precision :: t
       integer, intent(in) :: n, nlags
       double precision, dimension(n) :: y
       double precision, dimension(nlags) :: bval
       intent(in):: t,y
       intent(out) :: bval
     end subroutine f_beta_cc
  end interface
  interface
     subroutine f_change_fcn_cc(nevent,tevent,yevent,dyevent,hinit, &
          n_direction, direction,&
          n_isterminal, isterminal, quit) bind (c)
       integer :: nevent, n_direction, n_isterminal
       integer, dimension(n_direction) :: direction
       double precision :: tevent,hinit
       double precision, dimension(nevent) :: yevent,dyevent
       logical :: quit
       logical, dimension(n_isterminal) :: isterminal
       intent(in) :: nevent,tevent
       intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
     end subroutine f_change_fcn_cc
  end interface
  interface
     subroutine f_event_fcn_cc(t,n, nlags, y,dydt,z,g) bind (c)
       integer, intent(in) :: n, nlags
       double precision :: t
       double precision, dimension(n) :: y,dydt
       double precision, dimension(n,nlags) :: z
       double precision, dimension(n) :: g
       intent(in):: t,y,dydt,z
       intent(out) :: g
     end subroutine f_event_fcn_cc
  end interface
  interface
     subroutine f_out_fcn_cc(t,y,dy,n,nevent) bind (c)
       integer :: n,nevent
       double precision :: t
       double precision, dimension(n) :: y, dy
     end subroutine f_out_fcn_cc
  end interface
  interface
     subroutine f_user_trim_get_cc() bind (c)
     end subroutine f_user_trim_get_cc
  end interface

contains

  subroutine integrate_dde_1(&
       n_nvar, nvar, ddes_cc, beta_cc, history_cc,&
       n_tspan, tspan, &
       ! output
       sol_npts, sol_flag, sol_ne, &
       sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr,&
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
    integer(c_int), intent(out) :: sol_npts, sol_flag, sol_ne
    type(c_ptr) :: sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr

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
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     CHANGE_FCN = change_fcn, &
                     OUT_FCN = out_fcn, &
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     CHANGE_FCN = change_fcn, &
                     OUT_FCN = out_fcn)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     CHANGE_FCN = change_fcn, &
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     CHANGE_FCN = change_fcn)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     OUT_FCN = out_fcn, &
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     OUT_FCN = out_fcn)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn,&
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     EVENT_FCN = event_fcn)
             endif
          endif
       endif
    else
       if(have_change) then
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     CHANGE_FCN = change_fcn, &
                     OUT_FCN = out_fcn, &
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     CHANGE_FCN = change_fcn, &
                     OUT_FCN = out_fcn)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     CHANGE_FCN = change_fcn, &
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     CHANGE_FCN = change_fcn)
             endif
          endif
       else
          if(have_out) then
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     OUT_FCN = out_fcn, &
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     OUT_FCN = out_fcn)
             endif
          else
             if(have_trim) then
                sol = dkl_1(nvar,ddes,beta,history,tspan,options,&
                     USER_TRIM_GET = user_trim_get_cc)
             else
                sol = dkl_1(nvar,ddes,beta,history,tspan,options)
             endif
          endif
       endif
    endif

    sol_npts = sol%npts
    sol_flag = sol%flag
    sol_ne = sol%ne
    sol_t_ptr = c_loc(sol%t)
    sol_y_ptr = c_loc(sol%y)
    sol_te_ptr = c_loc(sol%te)
    sol_ye_ptr = c_loc(sol%ye)

    ! call print_stats(sol)
    ! call release_arrays(sol, options)

  contains
    subroutine ddes(t, y, z, dy) bind (c)
      double precision :: t
      double precision, dimension(:) :: y, dy
      double precision, dimension(:,:) :: z
      intent(in):: t,y,z
      intent(out) :: dy
      call ddes_cc(t, size(y), size(z, 2), y, z, dy)
    end subroutine ddes

    subroutine history(t,y) bind (c)
      double precision :: t
      double precision, dimension(:) :: y
      intent(in):: t
      intent(out) :: y
      call history_cc(t, size(y), y)
    end subroutine history

    subroutine beta(t, y, bval) bind (c)
      double precision :: t
      double precision, dimension(:) :: y
      double precision, dimension(:) :: bval
      intent(in):: t,y
      intent(out) :: bval
      call beta_cc(t, size(y), size(bval), y, bval)
    end subroutine beta

    subroutine change_fcn(nevent,tevent,yevent,dyevent,hinit, &
         direction,isterminal,quit) bind (c)
      integer :: nevent
      integer, dimension(:) :: direction
      double precision :: tevent,hinit
      double precision, dimension(:) :: yevent,dyevent
      logical :: quit
      logical, dimension(:) :: isterminal
      intent(in) :: nevent,tevent
      intent(inout) :: yevent,dyevent,hinit,direction,isterminal,quit
      call change_fcn_cc(nevent, tevent, yevent, dyevent, hinit, &
           size(direction), direction,&
           size(isterminal), isterminal, quit)
    end subroutine change_fcn

    subroutine event_fcn(t,y,dydt,z,g) bind (c)
      double precision :: t
      double precision, dimension(:) :: y,dydt
      double precision, dimension(:,:) :: z
      double precision, dimension(:) :: g
      intent(in):: t,y,dydt,z
      intent(out) :: g
      call event_fcn_cc(t, size(y), size(z, 2), y, dydt, z, g)
    end subroutine event_fcn

    subroutine out_fcn(t,y,dy,n,nevent) bind (c)
      integer :: n,nevent
      double precision :: t
      double precision, dimension(:) :: y, dy
    end subroutine out_fcn

  end subroutine integrate_dde_1

end module dde_solver_bind

