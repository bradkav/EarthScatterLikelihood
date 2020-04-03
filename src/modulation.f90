module modulation

use utils
use LabFuncs

implicit none
     
character(len=*) :: data_dir
parameter (data_dir = './DaMaSCUS_results/')
      
double precision :: vlag, sigmav, vesc
parameter (vlag = 230d0)
parameter (sigmav = 156d0)
parameter (vesc = 544d0)
      
!Sizes of grids in v, angle, etc.
integer :: Nv, Nang, Nsigma, Nmass
parameter (Nv = 200)
!parameter (Nang = 180)
parameter (Nang = 36)
parameter (Nsigma = 10)
parameter (Nmass = 4)

double precision :: vel_grid(Nv)
double precision :: eta_grid(Nmass,Nsigma,Nang, Nv)
double precision :: angle_list(Nang)
double precision :: rho_list(Nmass,Nsigma,Nang)
     
logical :: modulation_initialised = .False.

! These are the values of sigma used in Timon's simulation (1706.02249)
!parameter (sigma_vals = (/1d-45, 0.521d-36, 4.26d-36, 41.2d-36, 300.0d-36 /) )

double precision :: sigma_vals(Nsigma)
double precision :: mass_vals(Nmass)
parameter (sigma_vals = (/1d-40, 1d-38, 1d-37, 1d-36, 1d-35, 1d-34, 1d-33, 1d-32, 1d-31, 1d-30 /) )
parameter (mass_vals = (/0.100, 0.171, 0.292, 0.500 /) )
     
!character(len=4) :: simID(Nsigma)
!parameter (simID = (/"Free", "SS  ", "MS1 ", "MS10", "MS50"/) )
     
contains      
      
!=======================================================================

!MAKE SURE TO DOUBLE CHECK WHAT HAPPENS WHEN I SET THE CROSS SECTION BELOW 1d-33
!EVENT RATES LOOK HUUUUGE


!Load in the density and velocity tables
!and generate tables of eta (the velocity integral)
subroutine initialise_modulation
       
    integer :: i_m, i_sig, i_ang, i_v
    character(len=3) :: mstr
    character(len=3) :: sigstr
    character(len=17) :: simID
    
    ! Only initialise the first time...
    if (modulation_initialised .eqv. .False.) then

        write(*,*) "    Initialising velocity integrals..."
        !Initialise grid of velocities
        vel_grid = linspace(0.1d0, 775d0, Nv)

    
        !Initialise densities
        do i_m = 1, Nmass
            !First do the 'free' case (sigma << 1e-36)
            rho_list(i_m,1,:) = (/ (0.0, i_ang = 1, Nang) /) + 0.3
            !Then read the other results from file
            do i_sig = 2, Nsigma
                write(sigstr, '(i0)') int(LOG10(sigma_vals(i_sig)))
                write(mstr, '(i0)') nint(1000*mass_vals(i_m))
                simID = 'mDM_' // mstr // '_sig_1e' // sigstr
                !write(*,*) simID  
                call read_density(simID, angle_list, rho_list(i_m,i_sig,:))
            end do
        end do

    
	
        !Initialise velocity integrals
        do i_m = 1, Nmass
            do i_ang = 1, Nang
                ! Initialise 'free' velocity integral from function
                do i_v = 1, Nv
                    eta_grid(i_m,1,i_ang, i_v) = calcEta_free(vel_grid(i_v)) 
                end do
        
                ! Initialise the rest from file
                do i_sig = 2, Nsigma
                    write(sigstr, '(i0)') int(LOG10(sigma_vals(i_sig)))
                    write(mstr, '(i0)') nint(1000*mass_vals(i_m))
                    simID = 'mDM_' // mstr // '_sig_1e' // sigstr
                    call calc_eta(simID, i_ang, eta_grid(i_m,i_sig,i_ang,:))
                end do
        
            end do
        end do

        
        modulation_initialised = .True.
    end if
    
end subroutine initialise_modulation


!----------
subroutine read_density(root, angles, dens)
    character(len=*) :: root
    character(len=200) :: fname 
    double precision :: angles(Nang), dens(Nang)
    integer :: i

    fname = trim(data_dir) // trim(root) // trim('.rho')

    open (unit=99, file=fname, status='old', action='read')
    
    do i=1,Nang,1
        read(99,*) angles(i), dens(i)
     enddo
     close(99)
    
end subroutine read_density

!-----------
subroutine read_hist(root, angle, vlist, flist)
    character(len=*) :: root
    character(len=200) :: fname
    double precision, Pointer :: vlist(:), flist(:)
    integer :: i, nlines
    integer :: angle
    character(len=12) :: int_as_string

    write(int_as_string,'(i0)') angle-1

    fname = trim(data_dir) // trim(root) // '_histograms/speed.' // trim(int_as_string)
    !write(*,*) fname

    nlines = count_lines(fname)
    
    allocate(vlist(nlines))
    allocate(flist(nlines))

     
    open (91, file = fname ,action='read') 
    do i=1,nlines,1
        read(91,*) vlist(i), flist(i)
        !write(*,*) vlist(i), flist(i)
    end do
    close(91)
     
    vlist = vlist*3d5
    flist = flist/3d5
     
end subroutine read_hist


!----------
subroutine calc_eta(root, angle, etalist)
    character(len=*) :: root
    double precision, Pointer :: vlist(:), flist(:)
    double precision ::  fresamplist(Nv), etalist(Nv), eta_temp(Nv)
    integer :: angle

    call read_hist(root, angle, vlist, flist)

	!Resample the velocity distribution to get a grid
	!of fixed dimension
    fresamplist = linterp(vlist, flist, vel_grid)

    !Calculate the cumulative integral of f(v)/v
    eta_temp = cumtrapz(vel_grid, fresamplist/(vel_grid + 1d-20))
    
    deallocate(vlist)
    deallocate(flist)
    
    !Get eta
    etalist = -1.0*eta_temp + eta_temp(size(eta_temp))

end subroutine calc_eta


!----------
!Interpolate the value of rho, given sigma and the
!isodetection angle

function interp_rho_scalar(mass, sigma, angle) result(rho)
    double precision :: angle, sigma, mass
    double precision :: rho

    integer :: ind_ang, ind_sig, ind_m

    double precision :: xd, yd, zd
    
    double precision :: c00, c01, c10, c11, c0, c1
    
    double precision :: lsigma, lm
    
    lsigma = log10(sigma)
    lm = log10(mass)
    
    if ((mass < mass_vals(1)).or.(mass > mass_vals(Nmass))) then
        rho = 0
    end if

    if (sigma > sigma_vals(Nsigma)) then
        rho = 0
    else if (angle > 180) then
        rho = 0
    else
        ind_ang = ceiling(angle/5)

        ind_sig = 1
        do while (sigma_vals(ind_sig+1) < sigma)
            ind_sig = ind_sig + 1
        end do
        
        ind_m = 1
        do while (mass_vals(ind_m+1) < mass)
            ind_m = ind_m + 1
        end do
        
        xd = (lm - log10(mass_vals(ind_m)))/(log10(mass_vals(ind_m + 1)) - log10(mass_vals(ind_m)))
        yd = (lsigma - log10(sigma_vals(ind_sig)))/(log10(sigma_vals(ind_sig + 1)) - log10(sigma_vals(ind_sig)))
        !xd = (sigma - sigma_vals(ind_sig))/(sigma_vals(ind_sig + 1) - sigma_vals(ind_sig))
        zd = (angle - angle_list(ind_ang))/(angle_list(ind_ang + 1) - angle_list(ind_ang))

    
        
        !c0 = rho_list(ind_m, ind_sig, ind_ang)*(1-xd) + rho_list(ind_m, ind_sig+1, ind_ang)*xd
        !c1 = rho_list(ind_m, ind_sig, ind_ang+1)*(1-xd) + rho_list(ind_m, ind_sig+1, ind_ang+1)*xd
    
        c00 = rho_list(ind_m, ind_sig, ind_ang)*(1-xd) + rho_list(ind_m+1, ind_sig, ind_ang)*xd
        c01 = rho_list(ind_m, ind_sig, ind_ang+1)*(1-xd) + rho_list(ind_m+1, ind_sig, ind_ang+1)*xd
        c10 = rho_list(ind_m, ind_sig+1, ind_ang)*(1-xd) + rho_list(ind_m+1, ind_sig+1, ind_ang)*xd
        c11 = rho_list(ind_m, ind_sig+1, ind_ang+1)*(1-xd) + rho_list(ind_m+1, ind_sig+1, ind_ang+1)*xd
    
        c0 = c00*(1-yd) + c10*yd
        c1 = c01*(1-yd) + c11*yd
        
        rho = c0*(1-zd) + c1*zd
    
        !rho = c0*(1-yd) + c1*yd
    
    end if

end function interp_rho_scalar



!---------------
!Return an interpolation of eta...
! Be careful, I haven't done tons of bounds checking, 
! so if sigma < min(sigmavals) or angle < 0,
! I'm not sure what happens


function interp_eta_full(mass, sigma, angle, v)
    double precision :: mass, sigma, angle, v
    double precision :: interp_eta_full
    
    integer :: ind_m
    double precision :: xd
    double precision :: c0, c1
    double precision :: lm
    
    lm = log10(mass)
    
    if (mass > mass_vals(Nmass)) then
        interp_eta_full = 0 
    else if (mass < mass_vals(1)) then
            interp_eta_full = 0
    else
        
        ind_m = 1
        !write(*,*) angle, ind_ang
        do while (mass_vals(ind_m+1) < mass)
            ind_m = ind_m + 1
        end do
        
        xd = (lm - log10(mass_vals(ind_m)))/(log10(mass_vals(ind_m + 1)) - log10(mass_vals(ind_m)))
        !write(*,*) xd
        !write(*,*) ind_m
        interp_eta_full = interp_eta_scalar(ind_m, sigma, angle, v)*(1-xd) + interp_eta_scalar(ind_m+1, sigma, angle, v)*xd
    end if
    
end function

!Here, I should just interp over different values of the mass...
function interp_eta_scalar(ind_m, sigma, angle, v)
    double precision :: sigma, angle, v
    double precision :: interp_eta_scalar
    
    integer :: ind_v, ind_ang, ind_sig, ind_m
    
    double precision :: xd, yd, zd, vmax, vmin
    
    double precision :: c00, c01, c10, c11, c0, c1
    
    double precision :: lsigma
    
    lsigma = log10(sigma)
    
    vmin = vel_grid(1)
    vmax = vel_grid(Nv)
    
    if (sigma > sigma_vals(Nsigma)) then
        interp_eta_scalar = 0
    else if (v > vmax) then
        interp_eta_scalar = 0
    else if (angle > 180) then
        interp_eta_scalar = 0
    else
    
        ind_v = ceiling((v - vmin)*(size(vel_grid)-1)/(vmax - vmin))
        ind_ang = ceiling(angle/5)
        ind_sig = 1
        !write(*,*) angle, ind_ang
        do while (sigma_vals(ind_sig+1) < sigma)
            ind_sig = ind_sig + 1
        end do
        !write(*,*) sigma, ind_sig, sigma_vals(ind_sig)
    
        xd = (lsigma - log10(sigma_vals(ind_sig)))/(log10(sigma_vals(ind_sig + 1)) - log10(sigma_vals(ind_sig)))
        !write(*,*) xd
        !xd = (sigma - sigma_vals(ind_sig))/(sigma_vals(ind_sig + 1) - sigma_vals(ind_sig))
        yd = (angle - angle_list(ind_ang))/(angle_list(ind_ang + 1) - angle_list(ind_ang))
        zd = (v - vel_grid(ind_v))/(vel_grid(ind_v + 1) - vel_grid(ind_v))
        
        !Set up the trilinear interpolation
        c00 = eta_grid(ind_m, ind_sig, ind_ang, ind_v)*(1-xd) + eta_grid(ind_m, ind_sig+1, ind_ang, ind_v)*xd
        c01 = eta_grid(ind_m, ind_sig, ind_ang, ind_v+1)*(1-xd) + eta_grid(ind_m, ind_sig+1, ind_ang, ind_v+1)*xd
        c10 = eta_grid(ind_m, ind_sig, ind_ang+1, ind_v)*(1-xd) + eta_grid(ind_m, ind_sig+1, ind_ang+1, ind_v)*xd
        c11 = eta_grid(ind_m, ind_sig, ind_ang+1, ind_v+1)*(1-xd) + eta_grid(ind_m, ind_sig+1, ind_ang+1, ind_v+1)*xd
        
        c0 = c00*(1-yd) + c10*yd
        c1 = c01*(1-yd) + c11*yd
        
        interp_eta_scalar = c0*(1-zd) + c1*zd
    end if
    
end function interp_eta_scalar


!--------------------
! Calculate the isodetection angle (degrees) given the time
! and detector position on Earth
function calcIsoAngle(t, lon, lat)
	!time t in JulianDay
	
    double precision :: vs(3), vs_hat(3), rdet_hat(3)
    double precision :: t, lat, lon
    double precision :: calcIsoAngle

    vs = -LabVelocity(t, lon, lat)
    vs_hat = vs/sqrt(sum(vs**2))
	
    rdet_hat = (/0d0, 0d0, 1d0/)
    calcIsoAngle = (pi - acos(sum(vs_hat*rdet_hat)))*180d0/pi
  
end function calcIsoAngle


!----------------------
! Calculate the 'free' value of eta
function calcEta_free(v)
    double precision :: aplus, aminus, aesc
    double precision :: v
    double precision :: N, eta, calcEta_free
    
    aplus = min(v + vlag, vesc)/(sqrt(2d0)*sigmav)
    aminus = min(v - vlag, vesc)/(sqrt(2d0)*sigmav)
    aesc = vesc/(sqrt(2d0)*sigmav)

    N = 1.0/(erfun(aesc) - sqrt(2d0/pi)*(vesc/sigmav)*exp(-0.5*(vesc/sigmav)**2))
    
    eta = (0.5/vlag)*(erfun(aplus) - erfun(aminus))
    eta = eta - (1.0/(sqrt(pi)*vlag))*(aplus - aminus)*exp(-0.5*(vesc/sigmav)**2)
   
    if (eta < 0) then
        eta = 0
    end if
   
    calcEta_free = eta*N
   
end function calcEta_free 


end module modulation
