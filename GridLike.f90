program GridLike

use LabFuncs
use DDrate
use modulation
use utils
use like

implicit none

double precision, allocatable :: siglist(:), rholist(:)

integer :: i_loop, N_loop, i_loop2, N_loop2 ,ind
double precision, allocatable :: slhood(:), slhood_Eonly(:), slhood_counts(:)
 
double precision :: Cube(2)


write(*,*) "    Testing likelihood calculator..."

!Cross section and local DM density
Cube = (/ 0d0, 0.3d0 /)

N_loop = 30
allocate(siglist(N_loop))
siglist = 10**linspace(-37d0,-34d0, N_loop)

N_loop2 = 30
allocate(rholist(N_loop2))
rholist = 10**linspace(log10(1d-36), log10(1d-35), N_loop2)
!rholist = 10**linspace(log10(0.03d0), log10(3.0d0), N_loop2)

allocate(slhood(N_loop*N_loop2))
allocate(slhood_Eonly(N_loop*N_loop2))
allocate(slhood_counts(N_loop*N_loop2))

!Calculate likelihoods as a function of cross section
do i_loop = 1, N_loop
    write(*,*) i_loop, " of ", N_loop
    Cube(1) = siglist(i_loop)
    do i_loop2 = 1, N_loop2
        Cube(2) = rholist(i_loop2)/siglist(i_loop)
        !write(*,*) Cube(2)
        call slikelihood(Cube,slhood(i_loop + N_loop*(i_loop2-1)))
        !slhood_Eonly(i_loop) = 0
        !slhood_counts(i_loop) = 0
        call slikelihood_Eonly(Cube,slhood_Eonly(i_loop + N_loop*(i_loop2-1)))
        call slikelihood_counts(Cube,slhood_counts(i_loop + N_loop*(i_loop2-1)))
    end do
end do


!Save to file
write(*,*) "    Saving to file <likes.txt>..."

open (unit = 2, file = "likes.txt")
do i_loop = 1, N_loop
    do i_loop2 = 1, N_loop2
        ind = i_loop + N_loop*(i_loop2 - 1)
        write(2,*) siglist(i_loop), rholist(i_loop2), slhood(ind), slhood_Eonly(ind), slhood_counts(ind)
    end do
    
end do
close(2)


end program GridLike
