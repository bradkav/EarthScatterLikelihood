program testRates

use LabFuncs
use DDrate
use modulation
use utils
use expt

implicit none

integer :: i_t, i_sig

!integer :: N = 100

double precision :: day, dt

integer :: N_tbins = 50
integer :: N_sig = 35

double precision, allocatable :: t_list(:), t_edges(:), R_list(:,:), sigma_list(:)

double precision :: t0, m_x, sigma_SI, A

allocate(t_list(N_tbins))
allocate(t_edges(N_tbins+1))
allocate(R_list(N_tbins, N_sig))
allocate(sigma_list(N_sig))

call initialise_modulation()

m_x = 200d-3

sigma_list = 10**linspace(-37d0, -30d0, N_sig)

day = (23+56.0/60.0 + 4.1/3600)/24.0
t_edges = linspace(t_start(), t_start() + day, N_tbins+1)
dt = t_edges(2) - t_edges(1)


do i_t = 1, N_tbins
    t_list(i_t) = t_edges(i_t) + dt/2.0 - t_edges(1)
    
    do i_sig = 1, N_sig
        R_list(i_t, i_sig) = Nevents_short(E_th, E_max, &
              t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_list(i_sig))
    end do
end do


open (unit = 2, file = "../results/timeseries.txt")
do i_t = 1, N_tbins
    write(2,fmt="(E10.5)", advance="no") t_list(i_t)
    do i_sig = 1, N_sig
        write(2,fmt="(1x,E10.5)", advance = "no") R_list(i_t, i_sig)
    end do
    write(2,*) " "
end do
close(2)


end program testRates