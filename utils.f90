module utils

implicit none
      
double precision :: pi
parameter (pi = 3.1415926535d0)
	  
contains      
      
!=======================================================================


!--------
! Generate an array of linearly spaced values...
function linspace(startval, endval, npoints)
   double precision :: startval, endval, delta
   integer :: i,npoints
   double precision :: linspace(npoints)

   delta = (endval - startval)/(1.0*npoints - 1)
 
   linspace = startval + (/((i*delta),i=0,(npoints-1))/) 
end function linspace

 !--------
 ! Count the number of lines in a file.
function count_lines(fname)
    character(len=*) :: fname
    integer :: count_lines
	
    character(100) :: tmp
    
    count_lines = 0 
    open (90, file = fname) 
    do 
        read (90,*, END=10) tmp
        if (.not.(tmp(1:1) == '#')) then 
            count_lines = count_lines + 1 
        end if
    end do 
    10 CLOSE (90)
	
end function count_lines



!------------
! error function in double precision
function erfun(x)
      implicit double precision (a - h, o - z)
      double precision :: erfun, y
      integer :: i, k
      dimension a(0 : 64), b(0 : 64)
      data (a(i), i = 0, 12) / &
      0.00000000005958930743d0, -0.00000000113739022964d0, &
      0.00000001466005199839d0, -0.00000016350354461960d0, &
      0.00000164610044809620d0, -0.00001492559551950604d0, &
     0.00012055331122299265d0, -0.00085483269811296660d0, &
     0.00522397762482322257d0, -0.02686617064507733420d0, &
    0.11283791670954881569d0, -0.37612638903183748117d0, &
     1.12837916709551257377d0 / 
      data (a(i), i = 13, 25) / &
         0.00000000002372510631d0, -0.00000000045493253732d0, &
         0.00000000590362766598d0, -0.00000006642090827576d0, &
         0.00000067595634268133d0, -0.00000621188515924000d0, &
         0.00005103883009709690d0, -0.00037015410692956173d0, &
         0.00233307631218880978d0, -0.01254988477182192210d0, &
         0.05657061146827041994d0, -0.21379664776456006580d0, &
         0.84270079294971486929d0 / 
      data (a(i), i = 26, 38) / &
          0.00000000000949905026d0, -0.00000000018310229805d0, &
          0.00000000239463074000d0, -0.00000002721444369609d0, &
         0.00000028045522331686d0, -0.00000261830022482897d0, &
         0.00002195455056768781d0, -0.00016358986921372656d0, &
         0.00107052153564110318d0, -0.00608284718113590151d0, &
         0.02986978465246258244d0, -0.13055593046562267625d0, &
         0.67493323603965504676d0 / 
      data (a(i), i = 39, 51) / &
         0.00000000000382722073d0, -0.00000000007421598602d0, &
         0.00000000097930574080d0, -0.00000001126008898854d0, &
         0.00000011775134830784d0, -0.00000111992758382650d0, &
         0.00000962023443095201d0, -0.00007404402135070773d0, &
         0.00050689993654144881d0, -0.00307553051439272889d0, &
         0.01668977892553165586d0, -0.08548534594781312114d0, &
         0.56909076642393639985d0 / 
      data (a(i), i = 52, 64) / &
         0.00000000000155296588d0, -0.00000000003032205868d0, &
         0.00000000040424830707d0, -0.00000000471135111493d0, &
         0.00000005011915876293d0, -0.00000048722516178974d0, &
         0.00000430683284629395d0, -0.00003445026145385764d0, &
         0.00024879276133931664d0, -0.00162940941748079288d0, &
         0.00988786373932350462d0, -0.05962426839442303805d0, &
         0.49766113250947636708d0 / 
      data (b(i), i = 0, 12) / &
         -0.00000000029734388465d0, 0.00000000269776334046d0, &
         -0.00000000640788827665d0, -0.00000001667820132100d0, &
         -0.00000021854388148686d0, 0.00000266246030457984d0, &
         0.00001612722157047886d0, -0.00025616361025506629d0, &
         0.00015380842432375365d0, 0.00815533022524927908d0, &
         -0.01402283663896319337d0, -0.19746892495383021487d0, & 
         0.71511720328842845913d0 / 
      data (b(i), i = 13, 25) / &
         -0.00000000001951073787d0, -0.00000000032302692214d0, &
         0.00000000522461866919d0, 0.00000000342940918551d0, &
         -0.00000035772874310272d0, 0.00000019999935792654d0, &
         0.00002687044575042908d0, -0.00011843240273775776d0, &
         -0.00080991728956032271d0, 0.00661062970502241174d0, &
         0.00909530922354827295d0, -0.20160072778491013140d0, &
         0.51169696718727644908d0 / 
      data (b(i), i = 26, 38) / &
         0.00000000003147682272d0, -0.00000000048465972408d0, &
         0.00000000063675740242d0, 0.00000003377623323271d0, &
         -0.00000015451139637086d0, -0.00000203340624738438d0, &
         0.00001947204525295057d0, 0.00002854147231653228d0, &
         -0.00101565063152200272d0, 0.00271187003520095655d0, &
         0.02328095035422810727d0, -0.16725021123116877197d0, &
         0.32490054966649436974d0 / 
      data (b(i), i = 39, 51) / &
         0.00000000002319363370d0, -0.00000000006303206648d0, &
         -0.00000000264888267434d0, 0.00000002050708040581d0, &
         0.00000011371857327578d0, -0.00000211211337219663d0, &
         0.00000368797328322935d0, 0.00009823686253424796d0, &
         -0.00065860243990455368d0, -0.00075285814895230877d0, &
         0.02585434424202960464d0, -0.11637092784486193258d0, &
         0.18267336775296612024d0 / 
      data (b(i), i = 52, 64) / &
         -0.00000000000367789363d0, 0.00000000020876046746d0, &
         -0.00000000193319027226d0, -0.00000000435953392472d0, &
         0.00000018006992266137d0, -0.00000078441223763969d0, &
         -0.00000675407647949153d0, 0.00008428418334440096d0, &
         -0.00017604388937031815d0, -0.00239729611435071610d0, &
         0.02064129023876022970d0, -0.06905562880005864105d0, &
         0.09084526782065478489d0 / 
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t + &
             a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t +  &
             a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + &
             a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + &
             a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + &
             b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + &
             b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + &
             b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + &
             b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      erfun = y
end function erfun
!



!---------------------------------
!-------Integration routines----
!---------------------------------


! Trapezoidal integration
pure function trapz(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    double precision, intent(in)  :: x(:)         !! Variable x
    double precision, intent(in)  :: y(size(x))   !! Function y(x)
    double precision             :: r            !! Integral ∫y(x)·dx
    integer :: n
    
    n = size(x)
    ! Integrate using the trapezoidal rule
    r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
 end function trapz

! Cumulative trapezoidal integration
 pure function cumtrapz(x, y) result(r)
    !! Calculates the integral of an array y with respect to x using the trapezoid
    !! approximation. Note that the mesh spacing of x does not have to be uniform.
    double precision, intent(in)  :: x(:)         !! Variable x
    double precision, intent(in)  :: y(size(x))   !! Function y(x)
    double precision             :: r(size(x))           !! Integral ∫y(x)·dx
    integer :: i, n
	
	
    ! Integrate using the trapezoidal rule
    n = size(x)
    do i = 1, n
        r(i) = sum((y(1+1:i-0) + y(1+0:i-1))*(x(1+1:i-0) - x(1+0:i-1)))/2
    end do
	
end function cumtrapz




 !---------------------------------
 !-------Interpolation routines----
 !---------------------------------

 !--------------
 ! Linear interpolation
 ! Note that this routine does not check for monotonically
 ! ascending values...
function linterp(xdata, ydata, xvals) result(yvals)
    double precision, intent(in) :: xdata(:)
    double precision, intent(in) :: ydata(size(xdata))
    double precision :: xvals(:)
    double precision :: yvals(size(xvals))
	 
    double precision :: min_xdata,max_xdata, weight    
	 
    integer :: i, data_index
	 
    min_xdata = xdata(1)
    max_xdata = xdata(size(xdata))
	
    do i = 1, size(xvals)
        if (xvals(i) < min_xdata) then
            yvals(i) = 0
        else if (xvals(i) > max_xdata) then
            yvals(i) = 0 
        else
            data_index = 1 
            do while (xdata(data_index+1) < xvals(i))
                data_index = data_index + 1
            end do
            weight = (xvals(i) - xdata(data_index))/(xdata(data_index+1)-xdata(data_index))
            yvals(i) = ydata(data_index) + weight*(ydata(data_index + 1) - ydata(data_index))
        end if
    end do
	 	 
 end function linterp
 
 ! Same as linterp, but accepts a single, scalar value
 function linterp_scalar(xdata, ydata, xval) result(yval)
     double precision, intent(in) :: xdata(:)
     double precision, intent(in) :: ydata(size(xdata))
     double precision :: xval
     double precision :: yval
	 
     double precision :: min_xdata,max_xdata, weight    
	 
     integer :: data_index
	 
     min_xdata = xdata(1)
     max_xdata = xdata(size(xdata))
	
     
    if (xval < min_xdata) then
        yval = 0
    else if (xval > max_xdata) then
        yval = 0 
    else
        data_index = 1 
        do while (xdata(data_index+1) < xval)
            data_index = data_index + 1
        end do
        weight = (xval - xdata(data_index))/(xdata(data_index+1)-xdata(data_index))
        yval = ydata(data_index) + weight*(ydata(data_index + 1) - ydata(data_index))
    end if
	 	 
  end function linterp_scalar
 
  !---------------
  ! Grid interpolation - interpolation on a regular grid
  ! Should in principle be faster than linterp,
  ! but I don't really use it, it's probably not worth it
  function gridinterp(xdata, ydata, xvals) result(yvals)
      double precision, intent(in) :: xdata(:)
      double precision, intent(in) :: ydata(size(xdata))
      double precision :: xvals(:)
      double precision :: yvals(size(xvals))
     
      double precision :: min_xdata,max_xdata, weight, dx   
	 
      integer :: i, data_index
	 
      min_xdata = xdata(1)
      max_xdata = xdata(size(xdata))
      dx = (max_xdata - min_xdata)/float(size(xdata)-1)

      do i = 1, size(xvals)
          if (xvals(i) < min_xdata) then
              yvals(i) = 0
          else if (xvals(i) > max_xdata) then
              yvals(i) = 0 
          else
              data_index = ceiling((xvals(i) - min_xdata)/dx)
              weight = (xvals(i) - xdata(data_index))/dx
              yvals(i) = ydata(data_index) + weight*(ydata(data_index + 1) - ydata(data_index))
          end if
      end do
	 	 
   end function gridinterp
 
!---------------
function gridinterp_scalar(xdata, ydata, xval) result(yval)
    double precision, intent(in) :: xdata(:)
    double precision, intent(in) :: ydata(size(xdata))
    double precision :: xval
    double precision :: yval
   
    double precision :: min_xdata,max_xdata, weight, dx   
 
    integer :: data_index
 
    min_xdata = xdata(1)
    max_xdata = xdata(size(xdata))
    dx = (max_xdata - min_xdata)/float(size(xdata)-1)

    
 
    if (xval < min_xdata) then
        yval = 0
    else if (xval > max_xdata) then
        yval = 0 
    else
        data_index = ceiling((xval - min_xdata)/dx)
        weight = (xval - xdata(data_index))/dx
        yval = ydata(data_index) + weight*(ydata(data_index + 1) - ydata(data_index))
    end if
 	 
end function gridinterp_scalar
 



end module utils
