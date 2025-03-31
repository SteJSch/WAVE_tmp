*CMZ :  4.00/16 10/09/2022  09.39.11  by  Michael Scheer
*CMZ :  3.07/00 14/03/2019  18.19.35  by  Michael Scheer
*CMZ :  3.01/03 19/03/2014  12.24.14  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.69/00 24/10/2012  15.43.29  by  Michael Scheer
*CMZ :  2.66/19 07/06/2011  14.08.31  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  07.43.29  by  Michael Scheer
*CMZ :  2.63/03 07/05/2008  14.17.54  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.36  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.33.46  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.35.44  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  18.42.47  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.57  by  Michael Scheer
*-- Author : Michael Scheer
C***************************************************************
      SUBROUTINE SPLINZY_omp(N,XIN,Y,XA,YA,Y2A,KLO)
*KEEP,gplhint.
*KEND.

      IMPLICIT NONE
      INTEGER NOLD,N,KLO,KHI,KLOLD,K
      DOUBLE PRECISION XIN,Y,X,XA1OLD,XANOLD,H,A,B

      save klold,nold,xa1old,xanold

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.


      INTEGER max
      DOUBLE PRECISION XA(max(NDOBSVZP,NDOBSVYP))
     &      ,YA(max(NDOBSVZP,NDOBSVYP))
     &      ,Y2A(max(NDOBSVZP,NDOBSVYP))

      DATA KLOLD/0/
      DATA NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      IF (DABS(XIN-XA(1)).LT.1D-15) THEN
          X=XA(1)
      ELSEIF (DABS(XIN-XA(N)).LT.1D-15) THEN
          X=XA(N)
      ELSE
          X=XIN
      ENDIF

      IF(     XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &        .OR.
     &         XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1))) THEN
          WRITE(6,*) '*** ERROR IN SPLINZY: ARGUMENT OUT OF RANGE ***'
          STOP
      ENDIF

      KLO=1

      IF (X.LT.XA(KLO+1)) THEN
      KHI=KLO+1
      GOTO 2
      ENDIF

      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

2     H=XA(KHI)-XA(KLO)
      IF (H.LE.0.) THEN
        WRITE(6,*) '*** ERROR IN SPLINZY: BAD XA INPUT'
        STOP
      ENDIF
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)


      RETURN
      END
