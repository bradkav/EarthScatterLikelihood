module expt

use utils
use LabFuncs

implicit none
      
!Length of the exposure (in days)
!There are subroutines below for 
!calculating the start and end times
!in JulianDay format
double precision :: t_exp
parameter (t_exp = 30d0)

!Physical energy thresholds and resolution
!Minimum and maximum integration limits for the energy
double precision :: E_th, sigma_E, E_min, E_max, Nclip
            
!Here, E_th is the threshold energy
!Sigma_E is the energy resolution
!all in keV
!parameter (E_th = 100d-3)
parameter (E_th = 10d-3)
parameter (sigma_E = 3d-3)
!parameter (E_max = 1d0)
parameter (E_max = 0.25d0)

parameter (Nclip = 2d0)

!Don't change this...
!E_min is the minimum energy required for 
!integrals over all energies...
parameter (E_min = 1d-5)
            
!Background parameters
!Events per keV per kg per day (flat spectrum)
double precision :: BG_flat, BG_R0, BG_E0
parameter (BG_flat = 8.0d3)
!parameter (BG_rate = 10.0)
!parameter (BG_rate = 1d-5)
parameter (BG_R0 = 2.2d8)
parameter (BG_E0 = 25d-3) !25 eV characteristic energy

!Mass numbers and mass fractions
double precision :: A_det(2)
double precision :: mass_frac(2)
integer :: N_elements
!parameter (A_det = 73d0)

parameter (N_elements = 2)
parameter (A_det = (/ 16.0d0, 27.0d0 /))
parameter (mass_frac = (/ 0.53d0, 0.47d0 /))


!Latitude and longitude of detector
double precision :: lon_det = 13.576d0
double precision :: lat_det = 45.454d0

!Here we use the values for LNGS

!parameter (lon_det = 13.576d0)  
!parameter (lat_det = -45.454d0)
!parameter (lat_det = -45.454d0)

!Mass of detector (in kg)
double precision :: m_det
parameter (m_det = 35d-3)


contains

    
!----------
!Gaussian resolution smoothing function
function resolution(E1,E2)
    double precision :: E1, E2, resolution

    if ((E1 < 0).or.(E2 < 0).or.((E1-E2)**2 > (Nclip*sigma_E)**2)) then
        resolution = 0
        return
    end if

    resolution = exp(-(E1-E2)**2/(2.0*sigma_E**2))/sqrt(2.0*pi*sigma_E**2)

    !Clip below zero...
    !if (resolution < 0) then
    !    resolution = 0
    !end if


end function resolution
   
!-----------
!Gaussian resolution integrated over observed energies
function resolution_integrated(E, E1in, E2in)
    double precision :: E, E1, E2, resolution_integrated
    double precision, optional :: E1in, E2in
   
    if (present(E1in)) then
        E1 = E1in
    else
        E1 = E_th    
    end if
   
    if (present(E2in)) then
        E2 = E2in
    else
        E2 = E_max
    end if

    !Now account for truncation
    if (E < (E2 - Nclip*sigma_E)) then
       E2 = E + Nclip*sigma_E
    end if

    if (E > (E1 + Nclip*sigma_E)) then
       E1 = E - Nclip*sigma_E
    end if

    if (E1 < E2) then
       resolution_integrated = 0.5*(erfun((E2 - E)/(sqrt(2.0)*sigma_E)) + erfun((E - E1)/(sqrt(2.0)*sigma_E)))
    else
       resolution_integrated = 0.0
    end if
      
end function resolution_integrated    
  

!----------
function t_start()
    double precision :: t_start
    t_start = JulianDay(1, 1, 2018, 0d0)
end function t_start

!----------
function t_end()
    double precision :: t_end
    t_end = t_start() + t_exp
end function t_end

end module expt
