module like

use utils, only: count_lines
use DDrate
use expt !E_min, E_max, t_start(), t_end(), A_det, m_det, resolution


implicit none
      
!Comment out before compilation!
!Presumably these are defined elsewhere in MultiNest...
integer :: nest_nPar
parameter (nest_nPar = 3)
      
!List of events and number of events
double precision, allocatable :: events_E(:), events_t(:)
integer :: N_obs

!List of bin edges and event counts (for Asimov case)
double precision, allocatable :: E_edges(:), t_edges(:), N_binned(:,:)
integer :: N_lines, N_Ebins, N_tbins, i_E, i_t

logical :: loaded_data = .False.
      
contains      
      
!Load events from file
subroutine loadEvents(binned)
    integer :: i
    logical, optional :: binned
    
    !Load data during first evaluation
    if (loaded_data .eqv. .False.) then
    
        if (present(binned)) then
            if (binned) then
                N_lines = count_lines("events_Asimov.txt")
                N_Ebins = nint(sqrt(float(N_lines)))
                N_tbins = 1*N_Ebins
                allocate(E_edges(N_Ebins + 1))
                allocate(t_edges(N_tbins + 1))
                allocate(N_binned(N_tbins, N_Ebins))
                
                OPEN(89,FILE="events_Asimov.txt")
                do i = 1, N_lines
                    i_t = int((i-1)/(N_Ebins)) + 1
                    i_E = mod(i-1,(N_tbins)) + 1
                    READ(89,*) t_edges(i_t), E_edges(i_E), N_binned(i_t, i_E)
                end do
                CLOSE (89)
                !Add the final entry onto the list of bin edges
                t_edges(N_tbins + 1) = 2*t_edges(N_tbins) - t_edges(N_tbins - 1)
                E_edges(N_Ebins + 1) = 10**(2*log10(E_edges(N_Ebins)) - log10(E_edges(N_Ebins - 1)))
                !write(*,*) t_edges
                !write(*,*) E_edges
                
                !N_obs = sum(sum(N_binned,DIM=2))
            end if
        else

            N_obs = count_lines("events.txt")
    
            allocate(events_E(N_obs))
            allocate(events_t(N_obs))
    
            OPEN(89,FILE="events.txt")
            do i = 1, N_obs
                READ(89,*) events_t(i), events_E(i)
            end do
            CLOSE (89)
            
        end if
        loaded_data = .True.
    end if    
end subroutine loadEvents

      
!=======================================================================

!Wrapper for the main likelihood function
subroutine slikelihood(Cube,slhood)
    
    double precision Cube(nest_nPar),slhood
    
    !Here, we call the likelihood function including energy and timing info
    !binned=.True. specifies that we want to use the binned data in 'events_Asimov.txt'
    call loglike(Cube,slhood,binned=.True.)
    
end subroutine slikelihood


!=======================================================================

!Full likelihood including energy and times of events
!'Binned' is an optional argument
subroutine loglike(Cube,slhood,binned)
         
    implicit none
      
    double precision Cube(nest_nPar),slhood
    
    integer ::  i
    double precision :: N_exp
    
    double precision :: m_x, sigma_SI, rho
    double precision :: eventlike
    
    logical, optional :: binned
    
    m_x = Cube(1)
    sigma_SI = Cube(2)
    rho = Cube(3)
    
    !Load in events from file - but only the first time
    call loadEvents(binned)
    
    !Initialise modulation (if it hasn't been done already)
    call initialise_modulation()

    !Calculate binned likelihood
    if (present(binned)) then
        if (binned) then
        
            slhood = 0
            do i_t = 1, N_tbins
                do i_E = 1, N_Ebins
                    !Calculate expected number of events in bin
                    N_exp = t_exp*(rho/rho0)*m_det*Nevents_short(E_edges(i_E), E_edges(i_E + 1), &
                        t_edges(i_t), t_edges(i_t + 1), A_det, m_x, sigma_SI)
            
                    !Add background rate
                    N_exp = N_exp + BG_rate*(t_edges(i_t+1) - t_edges(i_t))* &
                      (E_edges(i_E+1) - E_edges(i_E))*t_exp*m_det
             
                    ! Calculate poisson likelihood
                    slhood = slhood - N_exp + N_binned(i_t, i_E)*log(N_exp)
                end do
            end do
            return
        end if
    end if
    
    !Otherwise, do a full event-by-event likelihood
    N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
    N_exp = N_exp + BG_rate*(E_max - E_th)*t_exp*m_det
    
    if (N_exp < 1d-30) then
        slhood = -1d30
        return
    end if
    
    ! Calculate poisson likelihood
    slhood = - N_exp + N_obs*log(N_exp)
    
    ! Calculate event-by-event likelihood
    do i = 1, N_obs
        
        eventlike = m_det*((rho/rho0)*dRdE_res(events_E(i), events_t(i), A_det, m_x, sigma_SI) + BG_rate)/N_exp
        if (eventlike < 1d-30) then
            slhood = -1d30
            return
        end if
        slhood = slhood + log(eventlike)
    end do
    
    
end subroutine loglike

!=======================================================================

!Likelihood including only event energies
!'Binned' is an optional argument
subroutine loglike_Eonly(Cube,slhood,binned)
         
    implicit none
      
    double precision Cube(nest_nPar),slhood
    
    logical, optional :: binned
    
    integer ::  i
    double precision :: N_exp
    
    double precision :: m_x, sigma_SI, rho
    double precision :: eventlike
    
    m_x = Cube(1)
    sigma_SI = Cube(2)
    rho = Cube(3)
    
    !Load in events from file - but only the first time
    call loadEvents(binned)
    
    !Initialise modulation (if it hasn't been done already)
    call initialise_modulation()
    
    !Calculate binned likelihood
    if (present(binned)) then
        if (binned) then
            
            slhood = 0
            do i_E = 1, N_Ebins
                !Calculate expected number of events in bin
                N_exp = t_exp*(rho/rho0)*m_det*Nevents_short(E_edges(i_E), E_edges(i_E + 1), &
                    t_edges(1), t_edges(N_tbins+1), A_det, m_x, sigma_SI)
        
                !Add background rate
                N_exp = N_exp + BG_rate*(E_edges(i_E+1) - E_edges(i_E))*t_exp*m_det
         
                ! Calculate poisson likelihood
                slhood = slhood - N_exp + sum(N_binned(:, i_E))*log(N_exp)
            end do
            return
        end if
    end if
    
    !Otherwise, do event-by-event likelihood
    N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
    N_exp = N_exp + BG_rate*(E_max - E_th)*t_exp*m_det
    
    if (N_exp < 1d-30) then
        slhood = -1d30
        return
    end if
    
    ! Calculate poisson likelihood
    slhood = - N_exp + N_obs*log(N_exp)
    
    ! Calculate event-by-by event likelihood
    do i = 1, N_obs
        eventlike = m_det*((rho/rho0)*dRdE_res_tint(events_E(i), t_start(), t_end(), A_det, m_x, sigma_SI) + BG_rate)/N_exp
        if (eventlike < 1d-30) then
            slhood = -1d30
            return
        end if
        slhood = slhood + log(eventlike)
    end do
    
end subroutine loglike_Eonly
      

!=======================================================================

!Likelihood including only event count
subroutine loglike_counts(Cube,slhood,binned)
         
    implicit none
      
    double precision Cube(nest_nPar),slhood
    logical, optional :: binned
    
    double precision :: N_exp, N_obs_float
    double precision :: m_x, sigma_SI, rho
    
    m_x = Cube(1)
    sigma_SI = Cube(2)
    rho = Cube(3)
    
    !Load in events from file - but only the first time
    call loadEvents(binned)
    
    !Initialise modulation (if it hasn't been done already)
    call initialise_modulation()
    
    
    N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
    N_exp = N_exp + BG_rate*(E_max - E_th)*t_exp*m_det
        
    if (N_exp < 1d-30) then
        slhood = -1d30
        return
    end if
    
    if (present(binned)) then
        if (binned) then
            N_obs_float = sum(sum(N_binned,DIM=2))
            ! Calculate poisson likelihood
            slhood = - N_exp + N_obs_float*log(N_exp)
            return
        end if
    end if
    
    slhood = - N_exp + N_obs*log(N_exp)
    
end subroutine loglike_counts


end module like
