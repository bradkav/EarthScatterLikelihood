program GenerateEvents_Asimov
    !Run from the command line
    ! ./GenerateEvents_Asimov M_X SIGMA_SI RHO
    ! where 
    !   - M_X is the WIMP mass in GeV
    !   - SIGMA_SI is the WIMP-nucleon xsec in cm^2
    !   - RHO is the local density in GeV/cm^3
    
use ranlib_poisson
use LabFuncs
use DDrate
use modulation
use utils
use expt !E_min, E_max, t_start(), t_end(), A_det, m_det, resolution

implicit none

integer :: N_Ebins
integer :: N_tbins
!Note that right now, the code requires
!that the number of bins be the same in E and t
!I'll fix that later...
parameter (N_Ebins = 12)
parameter (N_tbins = 12)

double precision :: m_x, sigma_SI, rho

double precision :: N_exp(N_tbins,N_Ebins)

double precision :: day

integer :: i_t, i_E, ind

!Initialise bins
double precision :: t_edges(N_tbins + 1)
double precision :: E_edges(N_Ebins + 1)

!Load in cross section and density from the command line
character(len=100) :: input_str

call getarg(1, input_str)
read(input_str, *) m_x

call getarg(2, input_str)
read(input_str,*) sigma_SI

call getarg(3, input_str)
read(input_str,*) rho

write(*,*) "Calculating Asimov data using:"
write(*,*) "    m_x [GeV]:", m_x
write(*,*) "    sigma_SI [cm^2]:", sigma_SI
write(*,*) "    rho [GeV/cm^3]:", rho


day = (23+56.0/60.0 + 4.1/3600)/24.0
t_edges = linspace(t_start(), t_start() + day, N_tbins+1)
E_edges = 10**linspace(log10(E_th), log10(E_max), N_Ebins+1)


!Read in the data tables for the daily modulations
call initialise_modulation

do i_t = 1, N_tbins
    do i_E = 1, N_Ebins
        N_exp(i_t, i_E) = (rho/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
          t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_SI)
          
          !Add background rate
          N_exp(i_t, i_E) = N_exp(i_t, i_E) + BG_rate*(t_edges(i_t+1) - t_edges(i_t))* &
              (E_edges(i_E+1) - E_edges(i_E))*t_exp*m_det
              
        if (N_exp(i_t, i_E) < 0) then
            N_exp(i_t, i_E) = 0
        end if
    end do
end do


!Save to file
write(*,*) "    Saving to file <events_Asimov.txt>..."

open (unit = 4, file = "events_Asimov.txt")
do i_t = 1, N_tbins
    do i_E = 1, N_Ebins
        !ind = i_t + N_tbins*(i_E - 1)
        write (4, *) t_edges(i_t), E_edges(i_E), N_exp(i_t, i_E)
    end do
    
end do
close(4)

end program GenerateEvents_Asimov
