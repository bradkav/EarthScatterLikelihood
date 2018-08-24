program test

use LabFuncs
use DDrate
use modulation
use utils
use like

implicit none

double precision :: siglist(100)

integer :: i_loop, N_loop
double precision :: slhood(100), slhood_Eonly(100), slhood_counts(100)
 
double precision :: Cube(2)


write(*,*) "    Testing likelihood calculator..."

!Cross section and local DM density
Cube = (/ 0d0, 0.3d0 /)

N_loop = 100
siglist = 10**linspace(-37d0,log10(300d-36),N_loop)

!Calculate likelihoods as a function of cross section
do i_loop = 1, N_loop
    Cube(1) = siglist(i_loop)
    call slikelihood(Cube,slhood(i_loop))
    call slikelihood_Eonly(Cube,slhood_Eonly(i_loop))
    call slikelihood_counts(Cube,slhood_counts(i_loop))
end do

!Save to file
write(*,*) "    Saving to file <likes.txt>..."

open (unit = 2, file = "likes.txt")
do i_loop = 1, N_loop
    write(2,*) siglist(i_loop), slhood(i_loop), slhood_Eonly(i_loop), slhood_counts(i_loop)
end do
close(2)


end program test
