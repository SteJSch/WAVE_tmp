*CMZ :          16/08/2024  10.07.55  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  11.03.55  by  Michael Scheer
*CMZ :  2.63/03 21/05/2008  13.37.36  by  Michael Scheer
*CMZ : 00.00/02 17/08/2004  09.47.26  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.25.29  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_INTEGRAL(X,Y,N,RESULT
     &                                 ,COEF,WORK1,WORK2,WORK3,WORK4)
*KEEP,gplhint.
*KEND.

C---  CALCULATES INTEGRAL OF Y(X) VIA SPLINES

      IMPLICIT NONE

      INTEGER I,N
      REAL*8 X(N),Y(N),RESULT
      REAL*8 COEF(N),WORK1(N),WORK2(N),WORK3(N),WORK4(N)

C---  SPLINE-COEFFICIENTS

      CALL UTIL_SPLINE_COEF(X,Y,N,0.0d0,0.0d0,COEF,WORK1,WORK2,WORK3,WORK4)

C--- INTEGRATION

      RESULT=0.0D0
      DO I=1,N-1

      RESULT=RESULT
     &          +(X(I+1)-X(I))*0.5D0
     &          *(Y(I)+Y(I+1))
     &          -(X(I+1)-X(I))**3/24.D0
     &          *(COEF(I)+COEF(I+1))

      ENDDO

      RETURN
      END
