* ======================================================================
* NIST Guide to Available Math Software. Not for redistribution.
* Fullsource for module DZERO from package PORT.
* Retrieved from NETLIB on Wed Feb  4 07:35:36 1998.
* ======================================================================
      DOUBLE PRECISION FUNCTION DZERO(F,A,B,T)
C
C  FINDS THE REAL ROOT OF THE FUNCTION F LYING BETWEEN A AND B
C  TO WITHIN A TOLERANCE OF
C
C         6*D1MACH(3) * DABS(DZERO) + 2 * T
C
C  F(A) AND F(B) MUST HAVE OPPOSITE SIGNS
C
C  THIS IS BRENTS ALGORITHM
C
C  A, STORED IN SA, IS THE PREVIOUS BEST APPROXIMATION (I.E. THE OLD B)
C  B, STORED IN SB, IS THE CURRENT BEST APPROXIMATION
C  C IS THE MOST RECENTLY COMPUTED POINT SATISFYING F(B)*F(C) .LT. 0
C  D CONTAINS THE CORRECTION TO THE APPROXIMATION
C  E CONTAINS THE PREVIOUS VALUE OF D
C  M CONTAINS THE BISECTION QUANTITY (C-B)/2
C
      DOUBLE PRECISION F,A,B,T,TT,SA,SB,C,D,E,FA,FB,FC,TOL,M,P,Q,R,S
      EXTERNAL F
      DOUBLE PRECISION D1MACH
C
      TT = T
      IF (T .LE. 0.0D0) TT = 10.D0*D1MACH(1)
C
      SA = A
      SB = B
      FA = F(SA)
      FB = F(SB)
      IF (FA .NE. 0.0D0) GO TO 5
      DZERO = SA
      RETURN
 5    IF (FB .EQ. 0.0D0) GO TO 140
C/6S
C     IF (DSIGN(FA,FB) .EQ. FA) CALL SETERR(
C    1   47H DZERO - F(A) AND F(B) ARE NOT OF OPPOSITE SIGN, 47, 1, 1)
C/7S
!!!!!!!!      IF (DSIGN(FA,FB) .EQ. FA) CALL SETERR(
!!!!!!!!     1 ' DZERO - F(A) AND F(B) ARE NOT OF OPPOSITE SIGN', 47, 1, 1)
      if (DSIGN(FA,FB) .EQ. FA) then
         write(*,*)  'F(A) AND F(B) ARE NOT OF OPPOSITE SIGN'
         stop
      end if
C/
C
 10   C  = SA
      FC = FA
      E  = SB-SA
      D  = E
C
C  INTERCHANGE B AND C IF DABS F(C) .LT. DABS F(B)
C
 20   IF (DABS(FC).GE.DABS(FB)) GO TO 30
      SA = SB
      SB = C
      C  = SA
      FA = FB
      FB = FC
      FC = FA
C
 30   TOL = 2.0D0*D1MACH(4)*DABS(SB)+TT
      M = 0.5D0*(C-SB)
C
C  SUCCESS INDICATED BY M REDUCES TO UNDER TOLERANCE OR
C  BY F(B) = 0
C
      IF ((DABS(M).LE.TOL).OR.(FB.EQ.0.0D0)) GO TO 140
C
C  A BISECTION IS FORCED IF E, THE NEXT-TO-LAST CORRECTION
C  WAS LESS THAN THE TOLERANCE OR IF THE PREVIOUS B GAVE
C  A SMALLER F(B).  OTHERWISE GO TO 40.
C
      IF ((DABS(E).GE.TOL).AND.(DABS(FA).GE.DABS(FB))) GO TO 40
      E = M
      D = E
      GO TO 100
 40   S = FB/FA
C
C  QUADRATIC INTERPOLATION CAN ONLY BE DONE IF A (IN SA)
C  AND C ARE DIFFERENT POINTS.
C  OTHERWISE DO THE FOLLOWING LINEAR INTERPOLATION
C
      IF (SA.NE.C) GO TO 50
      P = 2.0D0*M*S
      Q = 1.0D0-S
      GO TO 60
C
C  INVERSE QUADRATIC INTERPOLATION
C
 50   Q = FA/FC
      R = FB/FC
      P = S*(2.0D0*M*Q*(Q-R)-(SB-SA)*(R-1.0D0))
      Q = (Q-1.0D0)*(R-1.0D0)*(S-1.0D0)
 60   IF (P.LE.0.0D0) GO TO 70
      Q = -Q
      GO TO 80
 70   P = -P
C
C  UPDATE THE QUANTITIES USING THE NEWLY COMPUTED
C  INTERPOLATE UNLESS IT WOULD EITHER FORCE THE
C  NEW POINT TOO FAR TO ONE SIDE OF THE INTERVAL
C  OR WOULD REPRESENT A CORRECTION GREATER THAN
C  HALF THE PREVIOUS CORRECTION.
C
C  IN THESE LAST TWO CASES - DO THE BISECTION
C  BELOW (FROM STATEMENT 90 TO 100)
C
 80   S = E
      E = D
      IF ((2.0D0*P.GE.3.0D0*M*Q-DABS(TOL*Q)).OR.
     1    (P.GE.DABS(0.5D0*S*Q))) GO TO 90
      D = P/Q
      GO TO 100
 90   E = M
      D = E
C
C  SET A TO THE PREVIOUS B
C
 100  SA = SB
      FA = FB
C
C  IF THE CORRECTION TO BE MADE IS SMALLER THAN
C  THE TOLERANCE, JUST TAKE A  DELTA STEP  (DELTA=TOLERANCE)
C         B = B + DELTA * SIGN(M)
C
      IF (DABS(D).LE.TOL) GO TO 110
      SB = SB+D
      GO TO 130
C
 110  IF (M.LE.0.0D0) GO TO 120
      SB = SB+TOL
      GO TO 130
C
 120  SB = SB-TOL
 130  FB = F(SB)
C
C  IF F(B) AND F(C) HAVE THE SAME SIGN ONLY
C  LINEAR INTERPOLATION (NOT INVERSE QUADRATIC)
C  CAN BE DONE
C
      IF ((FB.GT.0.0D0).AND.(FC.GT.0.0D0)) GO TO 10
      IF ((FB.LE.0.0D0).AND.(FC.LE.0.0D0)) GO TO 10
      GO TO 20
C     
C***SUCCESS***
 140  DZERO = SB
      RETURN
      END
