module expt

use utils
use LabFuncs

implicit none
      
!Length of the exposure (in days)
!There are subroutines below for 
!calculating the start and end times
!in JulianDay format
double precision :: t_exp

parameter (t_exp = 365d0)

!Physical energy thresholds and resolution
!Minimum and maximum integration limits for the energy
double precision :: E_th, sigma_E, E_min, E_max
            
parameter (E_th = 100d-3)
parameter (sigma_E = 18d-3)
parameter (E_min = 1d-3)
parameter (E_max = 2d0)
            
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
    
!-----------
!Gaussian resolution
function resolution(E)
    double precision :: E, resolution
    
    resolution = 0.5*(1 + erfun((E - E_th)/(sqrt(2.0)*sigma_E)))
    
end function resolution    

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