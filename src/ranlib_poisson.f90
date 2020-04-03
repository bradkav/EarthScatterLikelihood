module ranlib_poisson
    
contains


function sexpo ( )

    !*****************************************************************************80
    !
    !! SEXPO samples the standard exponential distribution.
    !
    !  Discussion:
    !
    !   This procedure corresponds to algorithm SA in the reference.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 March 2013
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Brown, James Lovato.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Joachim Ahrens, Ulrich Dieter,
    !    Computer Methods for Sampling From the
    !    Exponential and Normal Distributions,
    !    Communications of the ACM,
    !    Volume 15, Number 10, October 1972, pages 873-882.
    !
    !  Parameters:
    !
    !    Output, real ( kind = 4 ) SEXPO, a random deviate from the standard
    !    exponential distribution.
    !
      implicit none

      real ( kind = 4 ) a
      integer ( kind = 4 ) i
      real ( kind = 4 ) q(8)
      real ( kind = 4 ) r4_uni_01
      real ( kind = 4 ) sexpo
      real ( kind = 4 ) u
      real ( kind = 4 ) umin
      real ( kind = 4 ) ustar

      save q

      data q / &
           0.6931472E+00, &
           0.9333737E+00, &
           0.9888778E+00, &
           0.9984959E+00, &
           0.9998293E+00, &
           0.9999833E+00, &
           0.9999986E+00, &
           0.9999999E+00 /

      a = 0.0E+00
      u = r4_uni_01 ( )

      do

        u = u + u

        if ( 1.0E+00 < u ) then
          exit
        end if

        a = a + q(1)

      end do

      u = u - 1.0E+00

      if ( u <= q(1) ) then
        sexpo = a + u
        return
      end if

      i = 1
      ustar = r4_uni_01 ( )
      umin = ustar

      do

        ustar = r4_uni_01 ( )
        umin = min ( umin, ustar )
        i = i + 1

        if ( u <= q(i) ) then
          exit
        end if

      end do

      sexpo = a + umin * q(1)

      return
end function sexpo
    
function snorm ( )

!*****************************************************************************80
!
!! SNORM samples the standard normal distribution.
!
!  Discussion:
!
!    This procedure corresponds to algorithm FL, with M = 5, in the reference.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2013
!
!  Author:
!
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Extensions of Forsythe's Method for Random
!    Sampling from the Normal Distribution,
!    Mathematics of Computation,
!    Volume 27, Number 124, October 1973, page 927-937.
!
!  Parameters:
!
!    Output, real ( kind = 4 ) SNORM, a random deviate from the distribution.
!
  implicit none

  real ( kind = 4 ) a(32)
  real ( kind = 4 ) aa
  real ( kind = 4 ) d(31)
  real ( kind = 4 ) h(31)
  integer ( kind = 4 ) i
  real ( kind = 4 ) r4_uni_01
  real ( kind = 4 ) s
  real ( kind = 4 ) snorm
  real ( kind = 4 ) t(31)
  real ( kind = 4 ) tt
  real ( kind = 4 ) u
  real ( kind = 4 ) ustar
  real ( kind = 4 ) w
  real ( kind = 4 ) y

  save a
  save d
  save h
  save t

  data a / &
        0.0000000E+00, 0.3917609E-01, 0.7841241E-01, 0.1177699E+00, &
        0.1573107E+00, 0.1970991E+00, 0.2372021E+00, 0.2776904E+00, &
        0.3186394E+00, 0.3601299E+00, 0.4022501E+00, 0.4450965E+00, &
        0.4887764E+00, 0.5334097E+00, 0.5791322E+00, 0.6260990E+00, &
        0.6744898E+00, 0.7245144E+00, 0.7764218E+00, 0.8305109E+00, &
        0.8871466E+00, 0.9467818E+00, 1.009990E+00,  1.077516E+00, &
        1.150349E+00,  1.229859E+00,  1.318011E+00,  1.417797E+00, &
        1.534121E+00,  1.675940E+00,  1.862732E+00,  2.153875E+00 /

  data d / &
        0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, &
        0.0000000E+00, 0.2636843E+00, 0.2425085E+00, 0.2255674E+00, &
        0.2116342E+00, 0.1999243E+00, 0.1899108E+00, 0.1812252E+00, &
        0.1736014E+00, 0.1668419E+00, 0.1607967E+00, 0.1553497E+00, &
        0.1504094E+00, 0.1459026E+00, 0.1417700E+00, 0.1379632E+00, &
        0.1344418E+00, 0.1311722E+00, 0.1281260E+00, 0.1252791E+00, &
        0.1226109E+00, 0.1201036E+00, 0.1177417E+00, 0.1155119E+00, &
        0.1134023E+00, 0.1114027E+00, 0.1095039E+00 /

  data h / &
        0.3920617E-01, 0.3932705E-01, 0.3950999E-01, 0.3975703E-01, &
        0.4007093E-01, 0.4045533E-01, 0.4091481E-01, 0.4145507E-01, &
        0.4208311E-01, 0.4280748E-01, 0.4363863E-01, 0.4458932E-01, &
        0.4567523E-01, 0.4691571E-01, 0.4833487E-01, 0.4996298E-01, &
        0.5183859E-01, 0.5401138E-01, 0.5654656E-01, 0.5953130E-01, &
        0.6308489E-01, 0.6737503E-01, 0.7264544E-01, 0.7926471E-01, &
        0.8781922E-01, 0.9930398E-01, 0.1155599E+00, 0.1404344E+00, &
        0.1836142E+00, 0.2790016E+00, 0.7010474E+00 /

  data t / &
        0.7673828E-03, 0.2306870E-02, 0.3860618E-02, 0.5438454E-02, &
        0.7050699E-02, 0.8708396E-02, 0.1042357E-01, 0.1220953E-01, &
        0.1408125E-01, 0.1605579E-01, 0.1815290E-01, 0.2039573E-01, &
        0.2281177E-01, 0.2543407E-01, 0.2830296E-01, 0.3146822E-01, &
        0.3499233E-01, 0.3895483E-01, 0.4345878E-01, 0.4864035E-01, &
        0.5468334E-01, 0.6184222E-01, 0.7047983E-01, 0.8113195E-01, &
        0.9462444E-01, 0.1123001E+00, 0.1364980E+00, 0.1716886E+00, &
        0.2276241E+00, 0.3304980E+00, 0.5847031E+00 /

  u = r4_uni_01 ( )
  if ( u <= 0.5E+00 ) then
    s = 0.0E+00
  else
    s = 1.0E+00
  end if
  u = 2.0E+00 * u - s
  u = 32.0E+00 * u
  i = int ( u )
  if ( i == 32 ) then
    i = 31
  end if
!
!  Center
!
  if ( i /= 0 ) then

    ustar = u - real ( i )
    aa = a(i)

    do

      if ( t(i) < ustar ) then

        w = ( ustar - t(i) ) * h(i)

        y = aa + w

        if ( s /= 1.0E+00 ) then
          snorm = y
        else
          snorm = -y
        end if

        return

      end if

      u = r4_uni_01 ( )
      w = u * ( a(i+1) - aa )
      tt = ( 0.5E+00 * w + aa ) * w

      do

        if ( tt < ustar ) then
          y = aa + w
          if ( s /= 1.0E+00 ) then
            snorm = y
          else
            snorm = -y
          end if
          return
        end if

        u = r4_uni_01 ( )

        if ( ustar < u ) then
          exit
        end if

        tt = u
        ustar = r4_uni_01 ( )

      end do

      ustar = r4_uni_01 ( )

    end do
!
!  Tail
!
  else

    i = 6
    aa = a(32)

    do

      u = u + u

      if ( 1.0E+00 <= u ) then
        exit
      end if

      aa = aa + d(i)
      i = i + 1

    end do

    u = u - 1.0E+00
    w = u * d(i)
    tt = ( 0.5E+00 * w + aa ) * w

    do

      ustar = r4_uni_01 ( )

      if ( tt < ustar ) then
        y = aa + w
        if ( s /= 1.0E+00 ) then
          snorm = y
        else
          snorm = -y
        end if
        return
      end if

      u = r4_uni_01 ( )

      if ( u <= ustar ) then
        tt = u
      else
        u = r4_uni_01 ( )
        w = u * d(i)
        tt = ( 0.5E+00 * w + aa ) * w
      end if

    end do

  end if

end function snorm

function ignpoi ( mu )

!*****************************************************************************80
!
!! IGNPOI generates a Poisson random deviate.
!
!  Discussion:
!
!    This procedure generates a single random deviate from a Poisson
!    distribution with given mean.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by Barry Brown, James Lovato.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Joachim Ahrens, Ulrich Dieter,
!    Computer Generation of Poisson Deviates
!    From Modified Normal Distributions,
!    ACM Transactions on Mathematical Software,
!    Volume 8, Number 2, June 1982, pages 163-179.
!
!  Parameters:
!
!    Input, real ( kind = 4 ) MU, the mean of the Poisson distribution 
!    from which a random deviate is to be generated.
!
!    Output, integer ( kind = 4 ) IGNPOI, a random deviate from
!    the distribution.
!
  implicit none

  real ( kind = 4 ), parameter :: a0 = -0.5E+00
  real ( kind = 4 ), parameter :: a1 =  0.3333333E+00
  real ( kind = 4 ), parameter :: a2 = -0.2500068E+00
  real ( kind = 4 ), parameter :: a3 =  0.2000118E+00
  real ( kind = 4 ), parameter :: a4 = -0.1661269E+00
  real ( kind = 4 ), parameter :: a5 =  0.1421878E+00
  real ( kind = 4 ), parameter :: a6 = -0.1384794E+00
  real ( kind = 4 ), parameter :: a7 =  0.1250060E+00
  real ( kind = 4 ) b1
  real ( kind = 4 ) b2
  real ( kind = 4 ) c
  real ( kind = 4 ) c0
  real ( kind = 4 ) c1
  real ( kind = 4 ) c2
  real ( kind = 4 ) c3
  real ( kind = 4 ) d
  real ( kind = 4 ) del
  real ( kind = 4 ) difmuk
  real ( kind = 4 ) e
  real ( kind = 4 ) fact(10)
  real ( kind = 4 ) fk
  real ( kind = 4 ) fx
  real ( kind = 4 ) fy
  real ( kind = 4 ) g
  integer ( kind = 4 ) ignpoi
  !integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kflag
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 4 ) mu
  !real ( kind = 4 ) muold
  !real ( kind = 4 ) muprev
  real ( kind = 4 ) omega
  real ( kind = 4 ) p
  real ( kind = 4 ) p0
  real ( kind = 4 ) px
  real ( kind = 4 ) py
  real ( kind = 4 ) q
  real ( kind = 4 ) r4_uni_01
  real ( kind = 4 ) s
  !real ( kind = 4 ) sexpo
  !real ( kind = 4 ) snorm
  real ( kind = 4 ) t
  real ( kind = 4 ) u
  real ( kind = 4 ) v
  real ( kind = 4 ) x
  real ( kind = 4 ) xx

  save fact

  data fact / 1.0E+00, 1.0E+00, 2.0E+00, 6.0E+00, 24.0E+00, &
    120.0E+00, 720.0E+00, 5040.0E+00, 40320.0E+00, 362880.0E+00 /
!
!  Start new table and calculate P0.
!
  if ( mu < 10.0E+00 ) then

    m = max ( 1, int ( mu ) )
    p = exp ( - mu )
    q = p
    p0 = p
!
!  Uniform sample for inversion method.
!
    do

      u = r4_uni_01 ( )
      ignpoi = 0

      if ( u <= p0 ) then
        return
      end if

!
!  Creation of new Poisson probabilities.
!
      do k = 1, 35
        p = p * mu / real ( k )
        q = q + p
        if ( u <= q ) then
          ignpoi = k
          return
        end if
      end do

    end do

  else
    s = sqrt ( mu )
    d = 6.0E+00 * mu * mu
    l = int ( mu - 1.1484E+00 )
!
!  Normal sample.
!
    g = mu + s * snorm ( )

    if ( 0.0E+00 <= g ) then

      ignpoi = int ( g )
!
!  Immediate acceptance if large enough.
!
      if ( l <= ignpoi ) then
        return
      end if
!
!  Squeeze acceptance.
!
      fk = real ( ignpoi )
      difmuk = mu - fk
      u = r4_uni_01 ( )

      if ( difmuk * difmuk * difmuk <= d * u ) then
        return
      end if

    end if
!
!  Preparation for steps P and Q.
!
    omega = 0.3989423E+00 / s
    b1 = 0.04166667E+00 / mu
    b2 = 0.3E+00 * b1 * b1
    c3 = 0.1428571E+00 * b1 * b2
    c2 = b2 - 15.0E+00 * c3
    c1 = b1 - 6.0E+00 * b2 + 45.0E+00 * c3
    c0 = 1.0E+00 - b1 + 3.0E+00 * b2 - 15.0E+00 * c3
    c = 0.1069E+00 / mu
    if ( 0.0E+00 <= g ) then

      kflag = 0

      if ( ignpoi < 10 ) then

        px = -mu
        py = mu ** ignpoi / fact(ignpoi+1)

      else

        del = 0.8333333E-01 / fk
        del = del - 4.8E+00 * del * del * del
        v = difmuk / fk

        if ( 0.25E+00 < abs ( v ) ) then
          px = fk * log ( 1.0E+00 + v ) - difmuk - del
        else
          px = fk * v * v * ((((((( a7 &
            * v + a6 ) &
            * v + a5 ) &
            * v + a4 ) &
            * v + a3 ) &
            * v + a2 ) &
            * v + a1 ) &
            * v + a0 ) - del
        end if

        py = 0.3989423E+00 / sqrt ( fk )

      end if
      x = ( 0.5E+00 - difmuk ) / s
      xx = x * x
      fx = -0.5E+00 * xx
      fy = omega * ((( c3 * xx + c2 ) * xx + c1 ) * xx + c0 )

      if ( kflag <= 0 ) then

        if ( fy - u * fy <= py * exp ( px - fx ) ) then
          return
        end if

      else

        if ( c * abs ( u ) <= py * exp ( px + e ) - fy * exp ( fx + e ) ) then
          return
        end if

      end if

    end if
!
!  Exponential sample.
!
    do
      e = sexpo ( )
      u = 2.0E+00 * r4_uni_01 ( ) - 1.0E+00
      if ( u < 0.0E+00 ) then
        t = 1.8E+00 - abs ( e )
      else
        t = 1.8E+00 + abs ( e )
      end if

      if ( t <= -0.6744E+00 ) then
        cycle
      end if

      ignpoi = int ( mu + s * t )
      fk = real ( ignpoi )
      difmuk = mu - fk

      kflag = 1
!
!  Calculation of PX, PY, FX, FY.
!
      if ( ignpoi < 10 ) then

        px = -mu
        py = mu ** ignpoi / fact(ignpoi+1)
 
      else

        del = 0.8333333E-01 / fk
        del = del - 4.8E+00 * del * del * del
        v = difmuk / fk

        if ( 0.25E+00 < abs ( v ) ) then
          px = fk * log ( 1.0E+00 + v ) - difmuk - del
        else
          px = fk * v * v * ((((((( a7 &
            * v + a6 ) &
            * v + a5 ) &
            * v + a4 ) &
            * v + a3 ) &
            * v + a2 ) &
            * v + a1 ) &
            * v + a0 ) - del
        end if

        py = 0.3989423E+00 / sqrt ( fk )

      end if
      x = ( 0.5E+00 - difmuk ) / s
      xx = x * x
      fx = -0.5E+00 * xx
      fy = omega * ((( c3 * xx + c2 ) * xx + c1 ) * xx + c0 )

      if ( kflag <= 0 ) then

        if ( fy - u * fy <= py * exp ( px - fx ) ) then
          return
        end if

      else

        if ( c * abs ( u ) <= py * exp ( px + e ) - fy * exp ( fx + e ) ) then
          return
        end if

      end if

    end do

  end if

end function ignpoi

end module ranlib_poisson
