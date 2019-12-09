program testStat

use stat

implicit none

integer :: i, j, nmx, i_E, i_t

double precision :: mx_min, mx_max, x, tol, sigma, rho, par, N_exp, sig_min, sig_max, pval

double precision :: dzero

double precision :: day

character(len=100) :: input_str, i_str, data_str

logical :: G_only

call getarg(1, input_str)
read(input_str,*) data

call getarg(2, input_str)
read(input_str,*) G_only

!Define time and energy bins
day = (23+56.0/60.0 + 4.1/3600)/24.0
t_edges = linspace(t_start(), t_start() + day, N_tbins+1)
E_edges = 10**linspace(log10(E_th), log10(E_max), N_Ebins+1)

!Read in the data tables for the daily modulations
call initialise_modulation

if (G_only) go to 50

!Benchmark local dark matter density
rho_b = 0.4d0

!Mass range
nmx = 20
mx_min = 0.1001d0 !with mx_min=0.1 it does not work
mx_max = 0.5d0

!Number of points in the p-value matrix
nrho = 500
nsig = 500
allocate(p_val(nrho,nsig))
allocate(rho_i(nrho))
allocate(sigma_j(nsig))

!Output files
write(*,*) "    Saving to files <stat_#.txt>..."

write(data_str,*) data
do i = 1, 30

   write(i_str,*) i

   open (unit = 10, file = "output/stat_p_"//trim(adjustl(data_str))//"_"//trim(adjustl(i_str))//".txt")
   open (unit = 9, file = "output/stat_rho_"//trim(adjustl(data_str))//"_"//trim(adjustl(i_str))//".txt")
   open (unit = 8, file = "output/stat_sigma_"//trim(adjustl(data_str))//"_"//trim(adjustl(i_str))//".txt")
   open (unit = 7, file = "output/stat_info_"//trim(adjustl(data_str))//"_"//trim(adjustl(i_str))//".txt")

   m_x = 0.2d0

   x=10.d0**(-36.d0 + 4d0*dble(i-1)/20.d0) 
   call p_value(x)
   
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

   write(*,*) x

end do
deallocate(p_val)
deallocate(rho_i)
deallocate(sigma_j)

stop
50 continue

m_x     =  0.2d0
sig_min = -35.d0
sig_max = -31.d0

open(unit=11, file="G_1.txt")
open(unit=22, file="G_2.txt")
open(unit=33, file="G_3.txt")
do j = 1, 1000

   rho    = 0.4d0
   sigma  = 10.d0**(sig_min + (sig_max-sig_min)*dble(j-1)/dble(1000-1))

   par    = 0.d0
   do i_t = 1, N_tbins

      N_exp  = 0.d0
      do i_E = 1, N_Ebins

         par = (rho/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
              t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma)

         N_exp = N_exp + par

      end do

      if (i_t.eq.1) then
         write(11,*) sigma, N_exp/(rho*sigma)
      end if

      if (i_t.eq.6) then
         write(22,*) sigma, N_exp/(rho*sigma)
      end if

      if (i_t.eq.12) then
         write(33,*) sigma, N_exp/(rho*sigma)
      end if

   end do

end do
close(11)
close(22)
close(33)

stop

end program testStat

