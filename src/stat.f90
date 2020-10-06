module stat

  use DDrate
  use expt !E_min, E_max, t_start(), t_end(), A_det, m_det, resolution

  implicit none

  double precision :: df, pnonc

  double precision :: m_x, rho_b, sig_b, nu_tot

  double precision, dimension(:,:), allocatable :: p_val, m_bf, likes_grid, N_events_tot
  double precision, dimension(:), allocatable :: rho_i, sigma_j, p_val_rho

  integer :: data,i_Eg,i_tg

  integer :: nrho,nsig

  integer :: N_Ebins
  integer :: N_tbins                                                                                                                                                                                                  
  parameter (N_Ebins = 12)
  parameter (N_tbins = 12)
  double precision :: t_edges(N_tbins + 1)
  double precision :: E_edges(N_Ebins + 1)

contains

  !NOTES
  ! In order to extend to include uncertainties on DM mass, we just need to profile
  ! I.e., we need to calculate, for each (sigma, rho), the maximum log-L as a function 
  ! of mass - this can be done as a grid scan (take, e.g. 5-10 values of m_x for each)
  ! point. Not *too* much slower. 


  function cdf_nc_chi_minhalf(x)
    double precision :: cdf_nc_chi_minhalf,x,p,q
    !Solve cdf_nc_chi_minhalf(x) = 0 to search for the median of the non central chi square distribution

    !cumchn gives p, the cumulative distribution function at x of the non central chi square distribution
    !df = number of degrees of freedom
    !q  = 1 - p is the p-value
    !pnonc is the noncentrality parameter
    call cumchn(x,df,pnonc,p,q)

    cdf_nc_chi_minhalf=p-0.5d0

  end function cdf_nc_chi_minhalf




  subroutine p_value(x)
    double precision :: x,y,p,q
    double precision :: dzero
    double precision :: sigma_b, sigma, rho, sig_min, sig_max, rho_min, rho_max
    double precision :: sigma_old, N_binned(N_tbins,N_Ebins), N_exp_old(N_tbins,N_Ebins), N_exp, par, par2, lam_A, lam_1D
    double precision :: N_BG(N_Ebins)
    double precision, allocatable :: sig_list(:)
    double precision, allocatable :: N_grid(:,:,:)
    double precision :: Lambda
    real :: start,finish
    integer :: i,j,l,m,i_t,i_E
    !For a dark matter particle mass m_x, it gives a matrix of p-values for rejecting a hypothetized 
    !(rho,sigma) pair when the actual cross section and local density are x and rho_b = 0.4 GeV/cm^3,
    !respectively

    allocate(sig_list(nsig))
    !allocate(mx_list(nmx))
    allocate(N_grid(N_tbins, N_Ebins, nsig))

    sigma_b = x
    sig_b   = sigma_b

    if ( (sigma_b.ge.1.d-38) .and. (sigma_b.le.1.d-30) ) then
       sig_min = 0.1d0 * sigma_b
       sig_max = 10d0 * sigma_b
       rho_min = 0.01d0
       rho_max = 1.d0
    else
       write(*,*) trim('sigma='),sigma_b,trim(' cm^2 is currently out of bounds')
       stop
    end if

    sigma_old = 0.d0
    N_exp_old = 0.d0
    N_binned  = 0.d0

    do i = 1, nsig
       sig_list(i) = 10**(log10(sig_min) + (log10(sig_max)-log10(sig_min))*dble(i-1)/dble(nsig-1))
    end do

    !Calculate background in each bin
    do i_E = 1, N_Ebins
       N_BG(i_E) = BG_R0*BG_E0*(exp(-E_edges(i_E)/BG_E0) - exp(-E_edges(i_E+1)/BG_E0))
       N_BG(i_E) = N_BG(i_E) + BG_flat*(E_edges(i_E + 1) - E_edges(i_E))
    end do
   
    call cpu_time(start)

    !Asimov data
    if (data.eq.1) then !Energy + time

       nu_tot = 0.d0
       par    = 0.d0
       do i_t = 1, N_tbins
          do i_E = 1, N_Ebins

             N_binned(i_t, i_E) = (rho_b/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                  t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_b)

             N_binned(i_t, i_E) = N_binned(i_t, i_E) + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

             if (N_binned(i_t, i_E) < 0) then
                N_binned(i_t, i_E) = 0.d0
             end if

             nu_tot = nu_tot + N_binned(i_t, i_E)

          end do
       end do

    else if (data.eq.2) then !Time only

       nu_tot = 0.d0
       par    = 0.d0
       do i_t = 1, N_tbins

          do i_E = 1, N_Ebins
             par = (rho_b/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                  t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_b)

             par = par + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

             N_binned(i_t, 1) = N_binned(i_t, 1) + par        
          end do

          if (N_binned(i_t, 1) < 0) then
             N_binned(i_t, 1) = 0.d0
          end if

          nu_tot = nu_tot + N_binned(i_t, 1)

       end do

    else if (data.eq.3) then !Energy only

       nu_tot = 0.d0
       par    = 0.d0
       do i_E = 1, N_Ebins

          do i_t = 1, N_tbins
             par = (rho_b/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                  t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_b)

             par = par + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

             N_binned(1, i_E) = N_binned(1, i_E) + par
          end do

          if (N_binned(1, i_E) < 0) then
             N_binned(1, i_E) = 0.d0
          end if

          nu_tot = nu_tot + N_binned(1, i_E)

       end do

    end if

    !Grid of expected event numbers:
    do j = 1, nsig
       !sigma      = sig_min + (sig_max-sig_min)*dble(j-1)/dble(nsig-1)     
       sigma = sig_list(j)

          do i_t = 1, N_tbins
             do i_E = 1, N_Ebins

                !if (sigma.eq.sigma_old) then
                !   N_exp = (rho/rho0)*N_exp_old(i_t, i_E)
                !else
                N_grid(i_t, i_E, j) = m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                     t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma)

                !N_exp_old(i_t, i_E) = (rho0/rho)*N_exp
                !end if

                if (N_grid(i_t, i_E, j) < 0) then
                   N_grid(i_t, i_E, j) = 1d-30
                end if

             end do
          end do

    end do



    !p-value matrix
    do i = 1, nrho

       rho      = rho_min + (rho_max-rho_min)*dble(i-1)/dble(nrho-1)
       rho_i(i) = rho

       lam_1D = -1d30
       do j = 1, nsig

          sigma = sig_list(j)
          sigma_j(j) = sigma

          !Noncentrality parameter generalising Eq. (20) in Cowan et al., Eur. Phys. J. C (2011) 71: 1554
          if (data.eq.1) then !Energy + time

             par    = 0.d0
             N_exp  = 0
             do i_t = 1, N_tbins
                do i_E = 1, N_Ebins

                      N_exp = (rho/rho0)*N_grid(i_t, i_E, j)
                      !write(*,*) N_exp
                      !N_exp_old(i_t, i_E) = (rho0/rho)*N_exp
                      !end if

                      N_exp = N_exp + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

                      if (N_exp < 0) then
                         N_exp = 0d0
                      end if

                      par = par + N_binned(i_t, i_E) - N_exp + N_binned(i_t, i_E)*log(N_exp/N_binned(i_t, i_E))

                end do
             end do

          else if (data.eq.2) then !Time only 

             par    = 0.d0
             par2   = 0.d0
             do i_t = 1, N_tbins

                N_exp  = 0.d0
                do i_E = 1, N_Ebins

                      !if (sigma.eq.sigma_old) then
                      !   par2 = (rho/rho0)*N_exp_old(i_t, i_E)
                      !else
                      par2 = (rho/rho0)*N_grid(i_t, i_E, j)
                      !N_exp_old(i_t, i_E) = (rho0/rho)*par2
                      !end if

                      par2 = par2 + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

                      N_exp = N_exp + par2

                end do

                if (N_exp < 0) then
                      N_exp = 0d0
                end if

                par = par + N_binned(i_t, 1) - N_exp + N_binned(i_t, 1)*log(N_exp/N_binned(i_t, 1))


             end do

          else if (data.eq.3) then !Energy only 

                par    = 0.d0
                par2   = 0.d0
                do i_E = 1, N_Ebins

                   N_exp  = 0.d0
                   do i_t = 1, N_tbins

                      !if (sigma.eq.sigma_old) then
                      !   par2 = (rho/rho0)*N_exp_old(i_t, i_E)
                      !else
                      par2 = (rho/rho0)*N_grid(i_t, i_E, j)
                      !N_exp_old(i_t, i_E) = (rho0/rho)*par2
                      !end if

                      par2 = par2 + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

                      N_exp = N_exp + par2

                   end do

                   if (N_exp < 0) then
                      N_exp = 0d0
                   end if

                   par = par + N_binned(1, i_E) - N_exp + N_binned(1, i_E)*log(N_exp/N_binned(1, i_E)) 

                end do

          end if

          lam_A  = par

          likes_grid(i,j) = lam_A
          
          Lambda = -2.d0 * lam_A

          pnonc = Lambda

          df = 2.d0

          if ( pnonc .gt. 100.d0 ) then

             q = 0.d0

          else

             y  = dzero( cdf_nc_chi_minhalf, 0.d0, 110.d0, 1d-10 )

             !cumchi gives the cumulative distribution function p at y for a chi square distribuition with df degrees of freedom.
             !q = 1 - p is the corresponding p-value
             call cumchi ( y, df, p, q )

          end if

          !Save output for plots
          p_val(i,j)    = q

          !Maximising over sigma 
          if (lam_A.gt.lam_1D) then
             lam_1D = lam_A
          end if

       end do

       Lambda = -2.d0 * lam_1D

       pnonc = Lambda
       df = 1.d0

       if ( pnonc .gt. 100.d0 ) then
          q = 0.d0
       else
          y  = dzero( cdf_nc_chi_minhalf, 0.d0, 110.d0, 1d-10 )
          call cumchi ( y, df, p, q )
       end if

       p_val_rho(i) = q

    end do

    deallocate(sig_list)
    deallocate(N_grid)
    call cpu_time(finish)
    !write(*,*) finish-start

  end subroutine p_value

  !-----------------------------------


  subroutine p_value_profiled(x)

    !This should be nmx = 250 for the full calculations
    integer, parameter :: nmx = 1000
    double precision :: x,y,p,q
    double precision :: dzero
    double precision :: sigma_b, sigma, rho, sig_min, sig_max, rho_min, rho_max, mx_min, mx_max, mx_test, N_tot
    double precision :: sigma_old, N_binned(N_tbins,N_Ebins), N_exp_old(N_tbins,N_Ebins), N_exp, par, par2, lam_A, lam_1D
    double precision :: Lambda, maxLambda
    double precision :: N_BG(N_Ebins)
    double precision, allocatable :: N_grid(:,:,:,:)
    double precision, allocatable :: sig_list(:), mx_list(:)
    double precision, allocatable :: mx_grid(:,:)
    real :: start,finish
    integer :: i,j,l,m,i_t,i_E
    !For a dark matter particle mass m_x, it gives a matrix of p-values for rejecting a hypothetized 
    !(rho,sigma) pair when the actual cross section and local density are x and rho_b = 0.4 GeV/cm^3,
    !respectively

    allocate(sig_list(nsig))
    allocate(mx_list(nmx))

    allocate(N_grid(N_tbins, N_Ebins, nmx, nsig))

    !write(*,*) x

    !write(*,*) m_x, x

    sigma_b = x
    sig_b   = sigma_b


    mx_min = 0.0581d0 !with mx_min=0.1 it does not work
    !write(*,*) "Reset mass limits"
    mx_max = 0.5d0

    if ( (sigma_b.ge.1.d-38) .and. (sigma_b.le.1.d-30) ) then
       !write(*,*) "Check that sigma limits don't matter!"
       sig_min = 0.1d0 * sigma_b
       sig_max = 10d0 * sigma_b
       !rho_min = (0.1d0/0.4d0)*sig_min
       !rho_max = (1.d0/0.4d0)*sig_max
       rho_min = 0.01d0
       rho_max = 1.0d0
    else
       write(*,*) trim('sigma='),sigma_b,trim(' cm^2 is currently out of bounds')
       stop
    end if

    do i = 1, nsig
       sig_list(i) = 10**(log10(sig_min) + (log10(sig_max)-log10(sig_min))*dble(i-1)/dble(nsig-1))
    end do

    !DEFINE M_X GRID HERE!!!
    
    do i = 1, nmx
       mx_list(i) = 10**(log10(mx_min) + (log10(mx_max)-log10(mx_min))*dble(i-1)/dble(nmx-1))
    end do

    
    !write(*,*) sig_list

    !Calculate background in each bin
    do i_E = 1, N_Ebins
       N_BG(i_E) = BG_R0*BG_E0*(exp(-E_edges(i_E)/BG_E0) - exp(-E_edges(i_E+1)/BG_E0))
       N_BG(i_E) = N_BG(i_E) + BG_flat*(E_edges(i_E + 1) - E_edges(i_E))
    end do
    
    
    sigma_old = 0.d0
    N_exp_old = 0.d0
    N_binned  = 0.d0

    call cpu_time(start)

    !First generate a grid of N_binned

    !Asimov data
    if (data.eq.1) then !Energy + time

       nu_tot = 0.d0
       par    = 0.d0
       do i_t = 1, N_tbins
          do i_E = 1, N_Ebins
             !write(*,*) rho_b/rho0, sigma_b
             N_binned(i_t, i_E) = (rho_b/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                  t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_b)

             if (N_binned(i_t, i_E) < 0) then
                N_binned(i_t, i_E) = 1.0d-30
             end if
             
             N_binned(i_t, i_E) = N_binned(i_t, i_E) + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det


             nu_tot = nu_tot + N_binned(i_t, i_E)

          end do
       end do

    else if (data.eq.2) then !Time only

       nu_tot = 0.d0
       par    = 0.d0
       do i_t = 1, N_tbins

          do i_E = 1, N_Ebins
             par = (rho_b/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                  t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_b)

             par = par + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

             N_binned(i_t, 1) = N_binned(i_t, 1) + par        
          end do

          if (N_binned(i_t, 1) < 0) then
             N_binned(i_t, 1) = 0.d0
          end if

          nu_tot = nu_tot + N_binned(i_t, 1)

       end do

    else if (data.eq.3) then !Energy only

       nu_tot = 0.d0
       par    = 0.d0
       do i_E = 1, N_Ebins

          do i_t = 1, N_tbins
             par = (rho_b/rho0)*m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                  t_edges(i_t), t_edges(i_t+1), A_det, m_x, sigma_b)

             par = par + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

             N_binned(1, i_E) = N_binned(1, i_E) + par
          end do

          if (N_binned(1, i_E) < 0) then
             N_binned(1, i_E) = 0.d0
          end if

          nu_tot = nu_tot + N_binned(1, i_E)

       end do

    end if

    write(*,*) "    Total number of signal events:", nu_tot

    !BJK: Use a "REFINE" variable to GOTO here, and use a different set of
    !mass values...

    
    !Grid of expected event numbers:
    do j = 1, nsig
       !sigma      = sig_min + (sig_max-sig_min)*dble(j-1)/dble(nsig-1)     
       sigma = sig_list(j)
       !write(*,*) j
       do l = 1, nmx
          mx_test = mx_list(l)
          !mx_test = mx_min + (mx_max-mx_min)*dble(l-1)/dble(nmx-1) 

          do i_t = 1, N_tbins
             do i_E = 1, N_Ebins

                !if (sigma.eq.sigma_old) then
                !   N_exp = (rho/rho0)*N_exp_old(i_t, i_E)
                !else
                N_grid(i_t, i_E, l, j) = m_det*t_exp*Nevents_short(E_edges(i_E), E_edges(i_E+1), &
                     t_edges(i_t), t_edges(i_t+1), A_det, mx_test, sigma)

                !N_exp_old(i_t, i_E) = (rho0/rho)*N_exp
                !end if

                if (N_grid(i_t, i_E, l, j) < 0) then
                   N_grid(i_t, i_E, l, j) = 1d-30
                end if

             end do
          end do

       end do
    end do



    !p-value matrix
    do i = 1, nrho

       rho      = (rho_min + (rho_max-rho_min)*dble(i-1)/dble(nrho-1))
       rho_i(i) = rho

       lam_1D = -1d30
       do j = 1, nsig

          sigma = sig_list(j)
          sigma_j(j) = sigma

          lam_A = -1d30
          N_tot = 0d0
          do l = 1, nmx

             mx_test = mx_list(l)
             !mx_test = mx_min + (mx_max-mx_min)*dble(l-1)/dble(nmx-1) 

             !Noncentrality parameter generalising Eq. (20) in Cowan et al., Eur. Phys. J. C (2011) 71: 1554
             if (data.eq.1) then !Energy + time

                par    = 0.d0
                N_exp  = 0
                do i_t = 1, N_tbins
                   do i_E = 1, N_Ebins

                      !if (sigma.eq.sigma_old) then
                      !   N_exp = (rho/rho0)*N_exp_old(i_t, i_E)
                      !else

                      N_exp = (rho/rho0)*N_grid(i_t, i_E, l, j)
                      !write(*,*) N_exp
                      !N_exp_old(i_t, i_E) = (rho0/rho)*N_exp
                      !end if

                      if (N_exp < 0) then
                         N_exp = 1.0d-30
                      end if
                      
                      N_exp = N_exp + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

                      N_tot = N_tot + N_exp
                      
                      par = par + N_binned(i_t, i_E) - N_exp + N_binned(i_t, i_E)*log(N_exp/N_binned(i_t, i_E))

                   end do
                end do

             else if (data.eq.2) then !Time only 

                par    = 0.d0
                par2   = 0.d0
                do i_t = 1, N_tbins

                   N_exp  = 0.d0
                   do i_E = 1, N_Ebins

                      !if (sigma.eq.sigma_old) then
                      !   par2 = (rho/rho0)*N_exp_old(i_t, i_E)
                      !else
                      par2 = (rho/rho0)*N_grid(i_t, i_E, l, j)
                      !N_exp_old(i_t, i_E) = (rho0/rho)*par2
                      !end if

                      par2 = par2 + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

                      N_exp = N_exp + par2

                   end do

                   if (N_exp < 0) then
                      N_exp = 0d0
                   end if
                   
                   N_tot = N_tot + N_exp
                   
                   par = par + N_binned(i_t, 1) - N_exp + N_binned(i_t, 1)*log(N_exp/N_binned(i_t, 1))

                end do

             else if (data.eq.3) then !Energy only 

                par    = 0.d0
                par2   = 0.d0
                do i_E = 1, N_Ebins

                   N_exp  = 0.d0
                   do i_t = 1, N_tbins

                      !if (sigma.eq.sigma_old) then
                      !   par2 = (rho/rho0)*N_exp_old(i_t, i_E)
                      !else
                      par2 = (rho/rho0)*N_grid(i_t, i_E, l, j)
                      !N_exp_old(i_t, i_E) = (rho0/rho)*par2
                      !end if

                      par2 = par2 + N_BG(i_E)*(t_edges(i_t+1) - t_edges(i_t))*t_exp*m_det

                      N_exp = N_exp + par2

                   end do

                   if (N_exp < 0) then
                      N_exp = 0d0
                   end if

                   N_tot = N_tot + N_exp
                   
                   par = par + N_binned(1, i_E) - N_exp + N_binned(1, i_E)*log(N_exp/N_binned(1, i_E)) 

                end do

             end if

             !sigma_old = sigma

             if (par.gt.lam_A) then
                m_bf(i,j) = mx_test
                N_events_tot(i,j) = N_tot
                lam_A = par
             end if
             !write(*,*) mx_test, sigma, rho, par, lam_A

             !lam_A  = par !This is indeed delta log L

          end do !End of the loop over mx
          !write(*,*) " "

          likes_grid(i, j) = lam_A
          
          Lambda = -2.d0 * lam_A

          !write(*,*) " "
          pnonc = Lambda

          df = 2.d0

          if ( pnonc .gt. 100.d0 ) then

             q = 0.d0

          else

             y  = dzero( cdf_nc_chi_minhalf, 0.d0, 110.d0, 1d-10 )

             !cumchi gives the cumulative distribution function p at y for a chi square distribuition with df degrees of freedom.
             !q = 1 - p is the corresponding p-value
             call cumchi ( y, df, p, q )

          end if

          !Save output for plots
          p_val(i,j)    = q

          !Maximising over sigma
          if (lam_A.gt.lam_1D) then
             lam_1D = lam_A
          end if

       end do

       Lambda = -2.d0 * lam_1D

       pnonc = Lambda
       df = 1.d0

       if ( pnonc .gt. 100.d0 ) then
          q = 0.d0
       else
          y  = dzero( cdf_nc_chi_minhalf, 0.d0, 110.d0, 1d-10 )
          call cumchi ( y, df, p, q )
       end if

       p_val_rho(i) = q
 
    end do

    deallocate(mx_list)
    deallocate(mx_grid)
    deallocate(sig_list)
    deallocate(N_grid)

    call cpu_time(finish)
    !write(*,*) finish-start

  end subroutine p_value_profiled

end module stat
