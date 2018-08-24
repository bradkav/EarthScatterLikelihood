module like

use utils, only: count_lines
use DDrate
use expt !E_min, E_max, t_start(), t_end(), A_det, m_det, resolution


implicit none
      
!Comment out before compilation!
!Presumably these are defined elsewhere in MultiNest...
integer :: nest_nPar
parameter (nest_nPar = 2)
      
!List of events and number of events
double precision, allocatable :: events_E(:), events_t(:)
integer :: N_obs

logical :: loaded_data = .False.
      
contains      
      
!Load events from file
subroutine loadEvents()
    integer :: i
    
    if (loaded_data .eqv. .False.) then
    
        !Load data during first evaluation
        N_obs = count_lines("events.txt")
    
        allocate(events_E(N_obs))
        allocate(events_t(N_obs))
    
        OPEN(89,FILE="events.txt")
        do i = 1, N_obs
            READ(89,*) events_t(i), events_E(i)
        end do
        CLOSE (89)
    
        loaded_data = .True.
    end if    
end subroutine loadEvents
      
!=======================================================================

!Full likelihood including energy and times of events
subroutine slikelihood(Cube,slhood)
         
    implicit none
      
    double precision Cube(nest_nPar),slhood
    
    integer ::  i
    double precision :: N_exp
    
    double precision :: m_x, sigma_SI, rho
    double precision :: eventlike
    
    m_x = 0.5d0
    sigma_SI = Cube(1)
    rho = Cube(2)
    
    !Load in events from file - but only the first time
    call loadEvents()
    
    !Initialise modulation (if it hasn't been done already)
    call initialise_modulation()

    
    N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
    
    if (N_exp < 1d-30) then
        slhood = -1d30
        return
    end if
    
    ! Calculate poisson likelihood
    slhood = - N_exp + N_obs*log(N_exp)
    
    ! Calculate event-by-by event likelihood
    do i = 1, N_obs
        eventlike = (rho/rho0)*m_det*dRdE_res(events_E(i), events_t(i), A_det, m_x, sigma_SI)/N_exp
        if (eventlike < 1d-30) then
            slhood = -1d30
            return
        end if
        slhood = slhood + log(eventlike)
    end do
    
    
end subroutine slikelihood
      
!=======================================================================

!Likelihood including only event energies
subroutine slikelihood_Eonly(Cube,slhood)
         
    implicit none
      
    double precision Cube(nest_nPar),slhood
    
    integer ::  i
    double precision :: N_exp
    
    double precision :: m_x, sigma_SI, rho
    double precision :: eventlike
    
    m_x = 0.5d0
    sigma_SI = Cube(1)
    rho = Cube(2)
    
    !Load in events from file - but only the first time
    call loadEvents()
    
    !Initialise modulation (if it hasn't been done already)
    call initialise_modulation()
    
    N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
    
    if (N_exp < 1d-30) then
        slhood = -1d30
        return
    end if
    
    ! Calculate poisson likelihood
    slhood = - N_exp + N_obs*log(N_exp)
    
    ! Calculate event-by-by event likelihood
    do i = 1, N_obs
        eventlike = (rho/rho0)*m_det*dRdE_res_tint(events_E(i), t_start(), t_end(), A_det, m_x, sigma_SI)/N_exp
        if (eventlike < 1d-30) then
            slhood = -1d30
            return
        end if
        slhood = slhood + log(eventlike)
    end do
    
end subroutine slikelihood_Eonly


!=======================================================================

!Likelihood including only event count
subroutine slikelihood_counts(Cube,slhood)
         
    implicit none
      
    double precision Cube(nest_nPar),slhood
    
    double precision :: N_exp
    
    double precision :: m_x, sigma_SI, rho
    
    m_x = 0.5d0
    sigma_SI = Cube(1)
    rho = Cube(2)
    
    !Load in events from file - but only the first time
    call loadEvents()
    
    !Initialise modulation (if it hasn't been done already)
    call initialise_modulation()
    
    N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
    
    if (N_exp < 1d-30) then
        slhood = -1d30
        return
    end if
    
    ! Calculate poisson likelihood
    slhood = - N_exp + N_obs*log(N_exp)
    
end subroutine slikelihood_counts


end module like
