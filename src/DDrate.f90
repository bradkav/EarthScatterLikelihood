module DDrate

use utils      
use modulation
use expt, only: resolution, resolution_integrated, lon_det, lat_det, sigma_E, E_min, E_max

implicit none

double precision :: rho0, amu
parameter (rho0 = 4d-1)
parameter (amu = 931.5d3) !Define conversion factor from amu-->keV


integer :: N_smooth
parameter (N_smooth = 20)
      
contains      
      
!=======================================================================

! Minimum velocity 
function vmin(E, A, m_x)
    double precision :: E,A,m_x,vmin
    !double precision :: hour,JulianDay
    
    double precision :: m_A, mu
    m_A = A*0.9315
    mu =  (m_A*m_x)/(m_A+m_x)
    vmin = 3d5*SQRT((E/1d6)*(m_A)/(2*mu*mu))
end function vmin
  
  
!Reduced mass - input A as nucleon number and m_x in GeV
function reduced_m(A, m_x)
    double precision :: A, m_x, m_A, reduced_m
    
    m_A = 0.9315*A
    reduced_m =  (m_A * m_x)/(m_A + m_x)
end function reduced_m  


    
! A helper function for calculating the prefactors to dRdE
function rate_prefactor(m_x)
    double precision :: m_x, rate_prefactor, mu
    
    mu = reduced_m(1d0, m_x)
    rate_prefactor = 4.34d41*rho0/(2.0*m_x*mu*mu)
end function rate_prefactor
 
    

! Convert an NREFT coupling to a cross section
! Coupling in GeV^-2 -> cross section in cm^2
function coupling_to_xsec(c, m_x)
    double precision :: c, m_x, coupling_to_xsec
  
    !0.197 GeV  = 1e13/cm
    ! -> GeV^-1 = 1.97e-14 cm
    coupling_to_xsec = (1.97e-14)**2*c**2*(reduced_m(1d0, m_x)**2)/pi
end function coupling_to_xsec
  
! Standard Helm Form Factor for SI scattering
function calcSIFormFactor(E, A)
    double precision :: calcSIFormFactor, E, A
    
    double precision :: q1, q2, s, a1, c, R1, x, J1, F

    !Convert recoil energy to momentum transfer q in keV
    q1 = SQRT(2*A*amu*E)
    !Convert q into fm^-1
    q2 = q1*(1d-12/1.97d-7)
    
    !Calculate nuclear parameters
    s = 0.9
    a1 = 0.52
    c = 1.23*(A**(1.0/3.0)) - 0.60
    R1 = SQRT(c*c + 7*pi*pi*a1*a1/3.0 - 5*s*s)

 
    x = q2*R1
    J1 = SIN(x)/x**2 - COS(x)/x
    F = 3*J1/x
 
    calcSIFormFactor = (F**2)*(EXP(-(q2*s)**2))
    
end function calcSIFormFactor
  

! Calculation of the free velocity integral
! now takes place in modulation.f90
!function calcEta_free(v)
!    double precision :: aplus, aminus, aesc
!    double precision :: v
!    double precision :: N, eta, calcEta_free
!    
!    aplus = min(v + vlag, vesc)/(sqrt(2d0)*sigmav)
!    aminus = min(v - vlag, vesc)/(sqrt(2d0)*sigmav)
!    aesc = vesc/(sqrt(2d0)*sigmav)
!
!    N = 1.0/(erfun(aesc) - sqrt(2d0/pi)*(vesc/sigmav)*exp(-0.5*(vesc/sigmav)**2))
!    
!    eta = (0.5/vlag)*(erfun(aplus) - erfun(aminus))
!    eta = eta - (1.0/(sqrt(pi)*vlag))*(aplus - aminus)*exp(-0.5*(vesc/sigmav)**2)
!   
!    if (eta < 0) then
!        eta = 0
!    end if
!   
!    calcEta_free = eta*N
!   
!end function calcEta_free 

 
! Spin-Independent recoil rate
! including energy and time
! depending on Earth-scattering effects
function dRdE(E, t, A, m_x, sigma_SI)
    double precision :: E, A, m_x, sigma_SI, dRdE, int_factor
    double precision :: t, angle

    angle = calcIsoAngle(t, lon_det, lat_det)
    !angle = 0.5*pi
    !angle = 0
    
    !For light DM, we set the form factor equal to 1
    !For consistency with DAMASCUS
    int_factor = sigma_SI*A**2 !*calcSIFormFactor(E, A)


    
    !
    dRdE = (interp_rho_scalar(m_x, sigma_SI, angle)/rho0)*rate_prefactor(m_x) &
        *int_factor*interp_eta_full(m_x, sigma_SI, angle, vmin(E, A, m_x))
end function dRdE

!Standard Spin-Independent recoil rate
! as above, but includes resolution correction
function dRdE_res(E, t, A, m_x, sigma_SI)
    double precision :: E, A, m_x, sigma_SI, dRdE_res
    double precision :: t
    
    double precision :: Elist(N_smooth), reslist(N_smooth), dRdElist(N_smooth)
    
    integer :: i
    
    Elist = 10**linspace(log10(max(E - sigma_E*5.0, 1d-4)), log10(E + sigma_E*5.0), N_smooth)
    
    
    do i = 1, N_smooth
        reslist(i) = resolution(Elist(i),E)
        
        if (Elist(i) < 0) then
            dRdElist(i) = 0
        else
            dRdElist(i) = dRdE(Elist(i), t, A, m_x, sigma_SI)
        end if
    end do
    
    dRdE_res = trapz(Elist, reslist*dRdElist)
end function dRdE_res


!Standard Spin-Independent recoil rate
! 'Free' (i.e. no Earth-scattering effects)
function dRdE_free(E, A, m_x, sigma_SI)
    double precision :: E, A, m_x, sigma_SI, dRdE_free, int_factor

    int_factor = sigma_SI*calcSIFormFactor(E, A)*(A**2)
    
    dRdE_free = rate_prefactor(m_x)*int_factor*calcEta_free(vmin(E, A, m_x))
end function dRdE_free

! Recoil rate
! Including resolution and integrating over time [t0, t1]
function dRdE_res_tint(E, t0, t1, A, m_x, sigma_SI)
    double precision :: A, m_x, sigma_SI, dRdE_res_tint
    double precision :: E, t0, t1
    integer :: Nt = 200

    double precision, allocatable :: tlist(:)
    double precision, allocatable :: integvals(:,:), tmp(:)
    
    double precision :: Elist(N_smooth), reslist(N_smooth)
    double precision :: dt
    
    integer :: i_t, i_E
    
    allocate(tlist(Nt))
    allocate(integvals(N_smooth,Nt))
    allocate(tmp(Nt))
    
    
    dt = (t1 - 1d0 - t0)/3d0
    
    !Note: you have to do more than 4 days of exposure for this to work...
    !Do 50 samples on each of 4 separate days
    do i_t = 0,3
        tlist((1+i_t*50):((i_t+1)*50)) = linspace(t0 + i_t*dt, t0 + i_t*dt + 1d0, 50)
    end do
    
    Elist = 10**linspace(log10(max(E - sigma_E*5.0, 1d-4)), log10(E + sigma_E*5.0), N_smooth) 
    
    !write(*,*) Elist
    do i_t = 1, Nt
        do i_E = 1, N_smooth
            integvals(i_E, i_t) = resolution(E, Elist(i_E))*dRdE(Elist(i_E), tlist(i_t), A, m_x, sigma_SI)
        end do
    end do
    
    !dRdE_res_tint = trapz(tlist, integvals)
    !Integrate over energies
    do i_t = 1, Nt
        tmp(i_t) = trapz(Elist, integvals(:, i_t))
    end do
    
    !Integrate over time
    dRdE_res_tint = trapz(tlist, tmp)
        
    deallocate(tlist)
    deallocate(integvals)
    deallocate(tmp)
        
end function dRdE_res_tint


! Total number of events in the energy range [E0, E1]
! and over time interval [t0, t1]
function Nevents(E0, E1, t0, t1, A, m_x, sigma_SI)
    double precision :: A, m_x, sigma_SI, Nevents
    double precision :: E0, E1, t0, t1
    integer :: Nt = 200
    integer :: NE = 100
    double precision, allocatable :: tlist(:), Elist(:), tmp(:)
    double precision, allocatable :: integvals(:,:)
    double precision :: dt
    
    integer :: i_t, i_E
    
    double precision, allocatable :: reslist(:)
    
    allocate(tlist(Nt))
    allocate(Elist(NE))
    allocate(tmp(Nt))
    allocate(integvals(NE, Nt))
    allocate(reslist(NE))
    
    dt = (t1 - 1d0 - t0)/3d0
    
    if ((t1 - t0) < 4d0) then
        write (*,*) "    WARNING! DDrate.f90: Nevents: Calculation may not work for exposures"
        write (*,*) "    shorter than a few days. Try Nevents_short instead..."
    end if
    
    !Note: you have to do more than 4 days of exposure for this to work...
    !Do 50 samples on each of 4 separate days
    do i_t = 0,3
        tlist((1+i_t*50):((i_t+1)*50)) = linspace(t0 + i_t*dt, t0 + i_t*dt + 1d0, 50)
    end do

    Elist = 10**(linspace(log10(E0), log10(E1), NE))
    
    do i_E = 1, NE
        reslist(i_E) = resolution_integrated(Elist(i_E))
    end do 
    
    ! Calculate the rate (including resolution) on a grid of (E, t) values
    do i_E = 1, NE
        do i_t = 1, Nt
            integvals(i_E, i_t) = reslist(i_E)*dRdE(Elist(i_E), tlist(i_t), A, m_x, sigma_SI)
        end do
    end do
    
    !Integrate over energies
    do i_t = 1, Nt
        tmp(i_t) = trapz(Elist, integvals(:, i_t))
    end do
    
    !Integrate over time
    Nevents = trapz(tlist, tmp)
        
    deallocate(tlist)
    deallocate(Elist)
    deallocate(tmp)
    deallocate(integvals)
    deallocate(reslist)
        
end function Nevents


! Total number of events in the energy range [E0, E1]
! and over a (short ~ days) time interval [t0, t1]
function Nevents_short(E0, E1, t0, t1, A, m_x, sigma_SI)
    double precision :: A, m_x, sigma_SI, Nevents_short
    double precision :: E0, E1, t0, t1
    integer :: Nt
    
    parameter (Nt = 12)
   
    
    double precision :: tlist(Nt), Elist(N_smooth), tmp(Nt)
    double precision :: integvals(N_smooth, Nt), reslist(N_smooth)
    double precision :: dt
    
    
    
    integer :: i_t, i_E
    tlist = linspace(t0, t1, Nt)

    Elist = 10**(linspace(log10(E_min), log10(E_max), N_smooth))
    
    do i_E = 1, N_smooth
        reslist(i_E) = resolution_integrated(Elist(i_E), E0, E1)
    end do 
    
    ! Calculate the rate (including resolution) on a grid of (E, t) values
    do i_E = 1, N_smooth
        do i_t = 1, Nt
            integvals(i_E, i_t) = reslist(i_E)*dRdE(Elist(i_E), tlist(i_t), A, m_x, sigma_SI)
        end do
    end do
    
    !Integrate over energies
    do i_t = 1, Nt
        tmp(i_t) = trapz(Elist, integvals(:, i_t))
    end do
    
    !Integrate over time
    Nevents_short = trapz(tlist, tmp)
        
end function Nevents_short

! Number of events in recoil range [E0, E1] at fixed time t
function Nevents_fixedt(E0, E1, t, A, m_x, sigma_SI)
    double precision :: A, m_x, sigma_SI, Nevents_fixedt
    double precision :: E0, E1, t
    integer :: NE = 100
    double precision, allocatable ::  Elist(:)
    double precision, allocatable :: integvals(:)
    
    integer :: i_E

    
    allocate(Elist(NE))
    allocate(integvals(NE))
    
    Elist = 10**(linspace(log10(E0), log10(E1), NE))
    
    do i_E = 1, NE
        integvals(i_E) = resolution_integrated(Elist(i_E))*dRdE(Elist(i_E), t, A, m_x, sigma_SI)
    end do
    
    Nevents_fixedt = trapz(Elist, integvals)
        
    deallocate(Elist)
    deallocate(integvals)
        
end function Nevents_fixedt



end module DDrate
