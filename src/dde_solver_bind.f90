module dde_solver_bind
  use iso_c_binding

  implicit none

  integer, parameter:: dp=kind(0.d0)

contains
  subroutine integrate_dde_1(&
       re, re_vector, n_re_vector, &
       ae, ae_vector, n_ae_vector, &
       hinit,hmax, n_isterminal, isterminal, &
       n_jumps, n_thit, n_direction, &
       direction,jumps,neutral,track_discontinuities,tracking_level, &
       interpolation,thit_exactly,max_events,max_steps, &
       moving_average,max_delay,trim_frequency,&
       n_nvar, nvar, ddes,beta,history,&
       n_tspan, tspan, event_fcn, &
       change_fcn,out_fcn,user_trim_get, &
       sol_npts, sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr) bind(c)
  use dde_solver_m
  implicit none

    double precision, optional :: ae, hinit, hmax, re, max_delay
    integer, optional :: max_events, max_steps, tracking_level, &
         moving_average, trim_frequency
    logical, optional :: interpolation, neutral, track_discontinuities
    integer, intent(in) :: n_re_vector, n_ae_vector, n_jumps, n_thit
    double precision, optional :: ae_vector(n_ae_vector), re_vector(n_re_vector)
    double precision, optional :: jumps(n_jumps), thit_exactly(n_thit)
    integer, intent(in) :: n_direction
    integer, optional :: direction(n_direction)
    integer, intent(in) :: n_isterminal
    logical, optional :: isterminal(n_isterminal)

    integer, intent(in) :: n_tspan, n_nvar
    double precision :: tspan(n_tspan)
    integer :: nvar(n_nvar)

    optional :: out_fcn, change_fcn, event_fcn, user_trim_get

    interface
       subroutine ddes(t,y,z,dy) bind (c)
         double precision :: t
         double precision, dimension(:) :: y,dy
         double precision, dimension(:,:) :: z
         intent(in):: t,y,z
         intent(out) :: dy
       end subroutine ddes
    end interface
   interface
       subroutine history(t,y) bind (c)
        double precision :: t
        double precision, dimension(:) :: y
        intent(in):: t
        intent(out) :: y
      end subroutine history
   end interface
   interface
      subroutine beta(t,y,bval) bind (c)
        double precision :: t
        double precision, dimension(:) :: y
        double precision, dimension(:) :: bval
        intent(in):: t,y
        intent(out) :: bval
      end subroutine beta
   end interface
    interface
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
       end subroutine change_fcn
    end interface
    interface
       subroutine event_fcn(t,y,dydt,z,g) bind (c)
         double precision :: t
         double precision, dimension(:) :: y,dydt
         double precision, dimension(:,:) :: z
         double precision, dimension(:) :: g
         intent(in):: t,y,dydt,z
         intent(out) :: g
       end subroutine event_fcn
    end interface
    interface
       subroutine out_fcn(t,y,dy,n,nevent) bind (c)
         integer :: n,nevent
         double precision :: t
         double precision, dimension(:) :: y,dy
       end subroutine out_fcn
    end interface
   interface
      subroutine user_trim_get() bind (c)
      end subroutine user_trim_get
   end interface

   integer, intent(out) :: sol_npts
   type(c_ptr), value :: sol_t_ptr, sol_y_ptr, sol_te_ptr, sol_ye_ptr

   ! local
   type (dde_opts) :: options
   type (dde_sol) :: sol
   
   options = dde_set(re,re_vector,ae,ae_vector,hinit,hmax,isterminal, &
        direction,jumps,neutral,track_discontinuities,tracking_level, &
        interpolation,thit_exactly,max_events,max_steps, &
        moving_average,max_delay,trim_frequency)

   sol = dkl_1(nvar,ddes,beta,history,tspan,options,event_fcn, &
        change_fcn,out_fcn,user_trim_get)

   sol_t_ptr = c_loc(sol%t)
   sol_y_ptr = c_loc(sol%y)
   sol_te_ptr = c_loc(sol%te)
   sol_ye_ptr = c_loc(sol%ye)

  end subroutine integrate_dde_1
  
end module dde_solver_bind

