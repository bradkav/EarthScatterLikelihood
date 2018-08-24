! Adapted from code by C. A. J. O'Hare
! Thanks Ciaran!


module LabFuncs

  implicit none
  double precision :: vv_earthrev,eccentricity,eccentricity_deg,orb_long_ecliptic,v_LSR
  double precision, private :: pi
  parameter (vv_earthrev = 29.8d0)
  parameter (eccentricity = 0.016722d0)
  parameter (eccentricity_deg = 0.9574d0)
  parameter (orb_long_ecliptic = 13.0d0+1.0d0)
  parameter (v_LSR = 220.0d0)
  parameter (pi = 3.1415926535d0)
  double precision,dimension(3),parameter :: lat_ecl_gal =  (/-5.5303d0,59.575d0,29.812d0/)
  double precision,dimension(3),parameter :: long_ecl_gal =  (/266.141d0,-13.3485d0,179.3212d0/)
  double precision,dimension(3),parameter :: v_pec = (/11.1d0,12.2d0,7.3d0/)




contains





!=================================Lab Velocity==========================================!

  function JulianDay(month,day,year,hour) ! Calculates time in JD for a given date
    integer :: month,day,year,year_r,month_r
    double precision :: hour,JulianDay
    year_r = year+4800-floor((14-month)/12.0)
    month_r = month+12*floor((14-month)/12.0)-3
    JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r &
         + floor(year_r/4.0) -floor(year_r/100.0) + &
         floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
  end function JulianDay


  function LabVelocity(JD, long, lat) ! Time in Julian days -> v_lab in km/s
    double precision :: JD ! input
    double precision :: LabVelocity(3) ! output
    double precision :: UT,MJD,T_0,t_GAST,t_lab, long, lat
    double precision :: v_galrot(3),v_solar(3),v_earthrot(3)
    double precision :: e,lambda_0,L,g,lambda_sun,beta(3),lambda_i(3),v_earthrev(3)

    ! Lab time conversion
    UT = 24*(JD+0.5-floor(JD+0.5)) ! Universal time
    MJD = JD - 2400000.5 ! Modified Julian Day
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + long/15
    t_lab = 15*t_lab ! Lab time in degrees


    ! Galactic (LSR) Rotation 
    v_galrot = (/0.0d0,v_LSR,0.0d0/)
    call gal2lab(v_galrot,t_lab, lat) ! transform to lab co-ords
    

    ! Peculiar solar Motion
    v_solar = v_pec
    call gal2lab(v_solar,t_lab, lat) ! transform to lab co-ords


    ! Earth's revolution (first calculate in galactic frame then transform)
    e = eccentricity
    lambda_0 = orb_long_ecliptic
    L = 281.0298 + 36000.77*T_0 + 0.04107*UT
    g = 357.9258 + 35999.05*T_0 + 0.04107*UT
    lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0)&
         + 0.020*sin(2*g*pi/180.0)
    beta = lat_ecl_gal
    lambda_i = long_ecl_gal
    v_earthrev = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))*&
         (cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
    call gal2lab(v_earthrev,t_lab, lat) ! transform to lab co-ords


    ! Earth's rotation
    v_earthrot = 0.465102*cos(lat*pi/180)*(/0.0,-1.0,0.0/) ! already in lab co-ords

    ! Total
    LabVelocity = v_earthrot+v_earthrev+v_solar+v_galrot
    
  end function LabVelocity





!=================================Co-ord transformations===============================!

! Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z) [TIME DEPENDENT]
  subroutine eqt2lab(v,t_lab, lat)
    double precision:: v(3),t_lab,t,vp(3),latr, lat
    t = t_lab*pi/180.0
    latr = lat*pi/180.0
    vp = v
    v(1) = -cos(t)*sin(latr)*vp(1) - sin(t)*sin(latr)*vp(2) + cos(latr)*vp(3)
    v(2) = sin(t)*vp(1) - cos(t)*vp(2)
    v(3) = cos(t)*cos(latr)*vp(1) + cos(latr)*sin(t)*vp(2) + sin(latr)*vp(3)
  end subroutine eqt2lab
 
! Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e) 
  subroutine gal2eqt(v)
    double precision :: v(3),vp(3)
    vp = v
    v(1) = -0.06699*vp(1) + 0.4927*vp(2) - 0.8676*vp(3)
    v(2) = -0.8728*vp(1) -0.4503*vp(2) -0.1884*vp(3)
    v(3) = -0.4835*vp(1) + 0.7446*vp(2) + 0.4602*vp(3)
  end subroutine gal2eqt

! Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z) [TIME DEPENDENT]
  subroutine gal2lab(v,t_lab, lat)
    double precision :: v(3),t_lab, lat
    call gal2eqt(v)
    call eqt2lab(v,t_lab, lat)
  end subroutine gal2lab

end module LabFuncs
