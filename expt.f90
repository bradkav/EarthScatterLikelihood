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
!E_min is the minimum energy required for 
!integrals over all energies...
parameter (E_th = 100d-3)
parameter (sigma_E = 20d-3)
parameter (E_min = 1d-4)
parameter (E_max = 2d0)
            
!Mass number
double precision :: A_det

parameter (A_det = 73d0)

!Latitude and longitude of detector
double precision :: lon_det, lat_det

!Here we use the values for LNGS
parameter (lon_det = 13.576d0)  
!parameter (lat_det = 45.454d0)
parameter (lat_det = -45.454d0)

!Mass of detector (in kg)
double precision :: m_det

parameter (m_det = 100d-3)


contains
    
!----------
!Gaussian resolution smoothing function
function resolution(E1,E2)
    double precision :: E1, E2, resolution
    
    resolution = exp(-(E1-E2)**2/(2.0*sigma_E**2))/sqrt(2.0*pi*sigma_E**2)

    !Clip below zero...
    if (resolution < 0) then
        resolution = 0
    end if

end function resolution
    
!-----------
!Gaussian resolution integrated over observed energies
function resolution_integrated(E)
    double precision :: E, resolution_integrated
    
    resolution_integrated = 0.5*(erfun((E_max - E)/(sqrt(2.0)*sigma_E)) + erfun((E - E_th)/(sqrt(2.0)*sigma_E)))
    
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