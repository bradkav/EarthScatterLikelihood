program calcContour

use stat

implicit none

integer :: i, j, nmx, i_E, i_t, N_bench

double precision :: mx_min, mx_max, x, tol, sigma_b, rho, par, N_exp, sig_min, sig_max, pval

double precision :: dzero

double precision :: day

character(len=100) :: input_str, i_str, data_str, outpath, m_str, sig_str, ID_str


call getarg(1, input_str)
read(input_str,*) m_x

call getarg(2, input_str)
read(input_str,*) sigma_b

!Data specifies which method to use!
! 1: Energy + Time
! 2: Time only
! 3: Energy only
call getarg(3, input_str)
read(input_str,*) data

call getarg(4, input_str)
read(input_str,*) outpath



!Define time and energy bins
day = (23+56.0/60.0 + 4.1/3600)/24.0
t_edges = linspace(t_start(), t_start() + day, N_tbins+1)
E_edges = 10**linspace(log10(E_th), log10(E_max), N_Ebins+1)



!Benchmark local dark matter density
rho_b = 0.4d0

!Mass range
!nmx = 20
!mx_min = 0.1001d0 !with mx_min=0.1 it does not work
!mx_max = 0.5d0


!Number of points in the p-value matrix
nrho = 5000
nsig = 500
allocate(p_val(nrho,nsig))
allocate(rho_i(nrho))
allocate(sigma_j(nsig))

!Output files

write(*,*) "    Benchmark"
write(*,*) "    ---------"
write(*,*) "    m_x [GeV]:", m_x
write(*,*) "    sigma [cm^2]:", sigma_b
write(*,*) " "

write(m_str,"(I3)") INT(m_x*1000d0)
write(sig_str, "(F6.2)") LOG10(sigma_b)
!do i = 1, N_bench

!write(*,*) outpath

!write(*,*) m_x, sigma_b

write(ID_str,*) "m"//trim(adjustl(m_str))//"_lsig"//trim(adjustl(sig_str))//".txt"

write(*,*) "    Saving to files <output/"//trim(adjustl(outpath))//"/stat_#.txt>..."
write(*,*) "    File ID: "//trim(adjustl(ID_str))



call system('mkdir -p output/' //trim(adjustl(outpath)))

open (unit = 10, file = "output/"//trim(adjustl(outpath))//"/stat_p_"//trim(adjustl(ID_str)))
open (unit = 9, file = "output/"//trim(adjustl(outpath))//"/stat_rho_"//trim(adjustl(ID_str)))
open (unit = 8, file = "output/"//trim(adjustl(outpath))//"/stat_sigma_"//trim(adjustl(ID_str)))
open (unit = 7, file = "output/"//trim(adjustl(outpath))//"/stat_info_"//trim(adjustl(ID_str)))

write(*,*) " "
!Read in the data tables for the daily modulations
call initialise_modulation

!m_x = 0.2d0

!x=10.d0**(-36.d0 + 4d0*dble(i-1)/20.d0) 
!x=10.d0**(sig_min + (sig_max-sig_min)*dble(i-1)/dble(N_bench-1)) 

write(*,*) "    Calculating p-value grid"
call p_value_profiled(1.0*sigma_b)

do j = 1, nrho
  write(10,'(500E20.12)') p_val(j,1:nsig)
end do

do j = 1, nrho
  write(9,*) rho_i(j)
end do

do j = 1, nsig
  write(8,*) sigma_j(j)
end do

write(7,*) m_x,x,nu_tot   

close(7)
close(8)
close(9)
close(10)


deallocate(p_val)
deallocate(rho_i)
deallocate(sigma_j)


end program calcContour

