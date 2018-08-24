program GenerateEvents
    !Run from the command line
    ! ./GenerateEvents SIGMA_SI RHO
    ! where 
    !   - SIGMA_SI is the WIMP-nucleon xsec in cm^2
    !   - RHO is the local density in GeV/cm^3
    
use ranlib_poisson
use LabFuncs
use DDrate
use modulation
use utils
use expt !E_min, E_max, t_start(), t_end(), A_det, m_det, resolution

implicit none
 
 
double precision :: N_exp
integer :: i, N_obs 

double precision :: m_x, sigma_SI, rho

double precision, allocatable :: events_E(:), events_t(:)
double precision, allocatable :: reslist(:)

!RNG stuff
integer :: values(1:8), k
integer, dimension(:), allocatable :: seed
double precision, allocatable :: x(:)
double precision :: r 

!Stuff for doing the inverse transform sampling
integer :: NE = 200
integer :: Nt
integer :: i_t, i_E
integer :: i_cut

!Lists for generating the inverse CDFs in (t, E)
double precision, allocatable :: tlist(:), Elist(:)
double precision, allocatable :: rate(:), cum_rate(:)

!Load in cross section and density from the command line
character(len=100) :: input_str
call getarg(1, input_str)
read(input_str,*) sigma_SI

call getarg(2, input_str)
read(input_str,*) rho

write(*,*) "Generating events using:"
write(*,*) "    sigma_SI [cm^2]:", sigma_SI
write(*,*) "    rho [GeV/cm^3]:", rho

!Fix DM mass
m_x = 5d-1


! Initialise seed for rng 
call initialize()
call date_and_time(values=values)
call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)
call random_number(r)
call set_initial_seed ( floor(1d6*r), floor(1d7*r) )

!Read in the data tables for the daily modulations
call initialise_modulation

! Calculate total number of events
N_exp = (rho/rho0)*m_det*Nevents(E_min, E_max, t_start(), t_end(), A_det, m_x, sigma_SI)
write(*,*) "    N_exp:", N_exp
N_obs = IGNPOI(real(N_exp))
write(*,*) "    N_obs:", N_obs

allocate(events_E(N_obs))
allocate(events_t(N_obs))
allocate(x(N_obs))

!----------------
write(*,*) "    Sampling event times..."

!Calculate the rate using roughly 100 time steps per day
Nt = 99*nint(t_end() - t_start()) + 1

allocate(tlist(Nt))
allocate(rate(Nt))
allocate(cum_rate(Nt))

! Sampling event rate as a function of time
! integrated over energy
tlist = linspace(t_start(), t_end(), Nt)
do i_t = 1, Nt
    rate(i_t) = (rho/rho0)*m_det*Nevents_fixedt(E_min, E_max, tlist(i_t), A_det, m_x, sigma_SI)
end do
    
! Inverse transform sample to get event times
cum_rate =  cumtrapz(tlist, rate)
cum_rate = cum_rate/cum_rate(Nt)
call random_number(x)

events_t = linterp(cum_rate, tlist, x)

deallocate(tlist)
deallocate(rate)
deallocate(cum_rate)


!----------------
write(*,*) "    Sampling event energies..."

allocate(Elist(NE))
allocate(reslist(NE))
allocate(rate(NE))
allocate(cum_rate(NE))

Elist = 10**(linspace(log10(E_min), log10(E_max), NE))

do i_E = 1, NE
    reslist(i_E) = resolution(Elist(i_E))
end do 

do i = 1, N_obs
    do i_E = 1, NE
        rate(i_E) = (rho/rho0)*m_det*reslist(i_E)*dRdE(Elist(i_E), events_t(i), A_det, m_x, sigma_SI)
    end do
    
    !Find the point where the rate drops to zero,
    !and cut there, otherwise the interpolation
    !can get screwed up...
    i_cut = 1
    do while (rate(i_cut) > 1d-30)
        i_cut = i_cut + 1
        if (i_cut > NE) then
            exit
        end if
    end do
    !In fact, cut one step earlier, just to be sure...
    i_cut = i_cut - 1
    
    cum_rate(1:i_cut) =  cumtrapz(Elist(1:i_cut), rate(1:i_cut))
    cum_rate(1:i_cut) = cum_rate(1:i_cut)/cum_rate(i_cut)
    
    call random_number(r)
    events_E(i) = linterp_scalar(cum_rate(1:i_cut), Elist(1:i_cut), r)
    
end do

!Save to file
write(*,*) "    Saving to file <events.txt>..."

open (unit = 2, file = "events.txt")
do i = 1, N_obs
    write (2, *) events_t(i), events_E(i)
end do
close(2)

end program GenerateEvents
