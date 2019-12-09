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
double precision :: E_th, sigma_E, E_min, E_max
            
!Here, E_th is the threshold energy
!Sigma_E is the energy resolution
!all in keV
parameter (E_th = 100d-3)
parameter (sigma_E = 25d-3)
parameter (E_max = 1d0)

!Don't change this...
!E_min is the minimum energy required for 
!integrals over all energies...
parameter (E_min = 1d-4)
            
!Background parameters
!Events per keV per kg per day (flat spectrum)
double precision :: BG_rate
parameter (BG_rate = 10.0)
!parameter (BG_rate = 1d-5)

!Mass number
double precision :: A_det
parameter (A_det = 73d0)

!Latitude and longitude of detector
double precision :: lon_det, lat_det

!Here we use the values for LNGS
parameter (lon_det = 13.576d0)  
parameter (lat_det = 45.454d0)
!parameter (lat_det = -45.454d0)

!Mass of detector (in kg)
double precision :: m_det
parameter (m_det = 100d-3)


contains
    
!----------
!Gaussian resolution smoothing function
function resolution(E1,E2)
    double precision :: E1, E2, resolution

    if ((E1 < 0).or.(E2 < 0)) then
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
    
    resolution_integrated = 0.5*(erfun((E2 - E)/(sqrt(2.0)*sigma_E)) + erfun((E - E1)/(sqrt(2.0)*sigma_E)))
    
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