program GridLike

use LabFuncs
use DDrate
use modulation
use utils
use like

implicit none

!double precision, allocatable :: siglist(:), rholist(:)

integer :: i_loop, N_loop, i_loop2, N_loop2 ,ind
!double precision, allocatable :: slhood(:), slhood_Eonly(:), slhood_counts(:)
 

 
double precision :: Cube(3)

parameter (N_loop = 50)
parameter (N_loop2 = 50)

double precision :: siglist(N_loop), rholist(N_loop2)

double precision :: slhood(N_loop*N_loop2), slhood_Eonly(N_loop*N_loop2), slhood_counts(N_loop*N_loop2)

write(*,*) "    Testing likelihood calculator..."

!Mass, Cross section and local DM density
Cube = (/ 0.35d0, 0d0, 0.3d0 /)

siglist = 10**linspace(-33d0,-30d0, N_loop)
!rholist(1) = 0.3
rholist = 0.3*10**linspace(-1d0, 1d0, N_loop2)


!Calculate likelihoods as a function of cross section
do i_loop = 1, N_loop
    write(*,*) i_loop, " of ", N_loop
    Cube(2) = siglist(i_loop)
    do i_loop2 = 1, N_loop2
        Cube(3) = rholist(i_loop2)!/siglist(i_loop)
        call loglike(Cube,slhood(i_loop + N_loop*(i_loop2-1)), binned=.True.)
        call loglike_Eonly(Cube,slhood_Eonly(i_loop + N_loop*(i_loop2-1)), binned=.True.)
        !call slikelihood_counts(Cube,slhood_counts(i_loop + N_loop*(i_loop2-1)))
    end do
end do


!Save to file
write(*,*) "    Saving to file <likes_grid.txt>..."

open (unit = 2, file = "likes_grid.txt")
do i_loop = 1, N_loop
    do i_loop2 = 1, N_loop2
        ind = i_loop + N_loop*(i_loop2 - 1)
        write(2,*) siglist(i_loop), rholist(i_loop2), slhood(ind), slhood_Eonly(ind), slhood_counts(ind)
    end do
    
end do
close(2)


end program GridLike
