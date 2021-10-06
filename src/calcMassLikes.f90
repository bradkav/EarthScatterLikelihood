program calcMassLikes

use stat

use expt, only: lat_det

implicit none

integer :: i, j, i_E, i_t, N_bench, FIX_MASS

double precision :: x, tol, sigma_b_in, sigma_test_in, rho, par, N_exp, sig_min, sig_max, pval

double precision :: dzero
double precision :: day

character(len=100) :: input_str, i_str, data_str, out_str, m_str, sig_str, ID_str, prof_str


call getarg(1, input_str)
read(input_str,*) m_x

call getarg(2, input_str)
read(input_str,*) sigma_b_in

call getarg(3, input_str)
read(input_str,*) sigma_test_in

call getarg(4, input_str)
read(input_str,*) lat_det

call getarg(5, input_str)
read(input_str,*) out_str



!Define time and energy bins
day = (23+56.0/60.0 + 4.1/3600)/24.0
t_edges = linspace(t_start(), t_start() + day, N_tbins+1)
!E_edges = 10**linspace(log10(E_th), log10(E_max), N_Ebins+1)
E_edges = linspace(E_th, E_max, N_Ebins+1)


!Benchmark local dark matter density
rho_b = 0.4d0

!Mass range
!nmx = 20
!mx_min = 0.1001d0 !with mx_min=0.1 it does not work
!mx_max = 0.5d0


!Number of points in the p-value matrix
nrho = 2000
nmx = 1000

allocate(rho_i(nrho))
allocate(mx_i(nmx))

allocate(likes_grid(nmx,nrho))

allocate(N_events_mx(nmx))

!Output files

write(*,*) "    Benchmark"
write(*,*) "    ---------"
write(*,*) "    m_x [GeV]:", m_x
write(*,*) "    sigma [cm^2]:", sigma_b_in
write(*,*) " "

write(m_str,"(I3)") INT(m_x*1000d0)
write(sig_str, "(F6.2)") LOG10(sigma_b_in)
!do i = 1, N_bench

!write(*,*) outpath

!write(*,*) m_x, sigma_b



write(ID_str,*) "example_mx"//trim(adjustl(m_str))//"_sig"//trim(adjustl(sig_str))//"_"//trim(adjustl(out_str))//".txt"

write(*,*) "    Saving to files <results/>"
write(*,*) "    File ID: "//trim(adjustl(ID_str))



call system('mkdir -p results/')

open (unit = 14, file = "results/"//trim(adjustl(ID_str)) )



!open( unit = 10, file = "p_test.txt")
!open( unit = 11, file = "like_test.txt")

write(*,*) " "
!Read in the data tables for the daily modulations
call initialise_modulation

!m_x = 0.2d0

!x=10.d0**(-36.d0 + 4d0*dble(i-1)/20.d0) 
!x=10.d0**(sig_min + (sig_max-sig_min)*dble(i-1)/dble(N_bench-1)) 


call massLikes(sigma_b_in, sigma_test_in)

do i = 1, nmx
    do j = 1, nrho
        write(14,*) mx_i(i), rho_i(j), likes_grid(i,j)
    end do
end do

!write(7,*) m_x,x,nu_tot   

close(14)

deallocate(likes_grid)
deallocate(rho_i)
deallocate(mx_i)

deallocate(N_events_mx)

end program calcMassLikes

