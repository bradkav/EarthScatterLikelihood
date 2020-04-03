program test

use LabFuncs
use DDrate
use modulation
use utils
use like

implicit none

integer N_loop
parameter (N_loop = 200)

double precision siglist(N_loop)

integer :: i_loop
double precision :: slhood(N_loop), slhood_Eonly(N_loop), slhood_counts(N_loop)
 
double precision :: Cube(3)


write(*,*) "    Testing likelihood calculator..."

!Mass, Cross section and local DM density
Cube = (/ 0.35d0, 0d0, 0.3d0 /)


siglist = 10**linspace(-34d0,-30d0,N_loop)




!Calculate likelihoods as a function of cross section
do i_loop = 1, N_loop
    !write(*,*) i_loop
    Cube(2) = siglist(i_loop)
    call loglike(Cube,slhood(i_loop),binned=.True.)
    call loglike_Eonly(Cube,slhood_Eonly(i_loop),binned=.True.)
    call loglike_counts(Cube,slhood_counts(i_loop),binned=.True.)
end do

!do i_loop = 1, N_loop
!    write(*,*) siglist(i_loop), interp_rho_scalar(0.3d0, siglist(i_loop), 25d0)   
!end do

!Save to file
write(*,*) "    Saving to file <likes.txt>..."

open (unit = 2, file = "likes.txt")
do i_loop = 1, N_loop
    write(2,*) siglist(i_loop), slhood(i_loop), slhood_Eonly(i_loop), slhood_counts(i_loop)
end do
close(2)


end program test
