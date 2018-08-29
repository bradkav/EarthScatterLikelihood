program testSpectra

use LabFuncs
use DDrate
use modulation
use utils
use like

implicit none

integer :: i

integer :: N = 100

double precision, allocatable :: Elist(:), dRdE_exact(:), dRdE_smooth(:), dRdE_tint(:)

double precision :: t0, m_x, sigma_SI, A

call initialise_modulation()

allocate(Elist(N))
allocate(dRdE_exact(N))
allocate(dRdE_smooth(N))
allocate(dRdE_tint(N))

write(*,*) "    Calculating some test spectra..."

Elist = 10**linspace(-3d0, log10(2d0), N)

t0 = 1d5
A = 73d0
m_x = 0.5d0
sigma_SI = 1d-34

write(*,*) Nevents_fixedt(1d-4, 2d0, t0, A, m_x, sigma_SI)

do i = 1, N
    dRdE_exact(i) = dRdE(Elist(i), t0, A, m_x, sigma_SI)
    dRdE_smooth(i) = dRdE_res(Elist(i), t0, A, m_x, sigma_SI)
    dRdE_tint(i) = dRdE_res_tint(Elist(i), t0,t0+30d0, A, m_x, sigma_SI)
end do

!Save to file
write(*,*) "    Saving to file <spectra.txt>..."

open (unit = 2, file = "spectra.txt")
do i = 1, N
    write(2,*) Elist(i), dRdE_exact(i), dRdE_smooth(i), dRdE_tint(i)
end do
close(2)


end program testSpectra
