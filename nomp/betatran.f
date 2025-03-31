*CMZ :  3.05/23 26/11/2018  21.53.32  by  Michael Scheer
*CMZ :  3.05/11 15/08/2018  15.13.49  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.23.35  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.02  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.35  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BETATRAN

*KEEP,gplhint.
*KEND.

C--- SUBROUTINE CALCULATES HORIZONTAL LINEAR TRANSFER-MATRIX FROM BETA-FUNCTION
C--- CALCULATION IS PERFORMED ONLY FOR POSITIVE X-VALUES, I.E. ONE HALF OF THE
C--- WLS
      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,betawls.
      include 'betawls.cmn'
*KEEP,depola.
      include 'depola.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER I
      DOUBLE PRECISION DELPHI,ALPHA,BETA,ALPHA0,T11,T12,T22,T21,COSPHI,SINPHI,B0,DX

C--- SUM UP

      DELPHI=0.
      DO I=INULLP,MBETA
        IF(I.EQ.1) THEN
          DX=XBETA(2)-XBETA(1)
        ELSE IF(I.EQ.MBETA) THEN
          DX=XBETA(MBETA)-XBETA(MBETA-1)
        ELSE
          DX=0.5*(XBETA(I+1)-XBETA(I-1))
        ENDIF
        DELPHI=DELPHI+1/YBETA(I)*DX
      ENDDO

      COSPHI=DCOS(DELPHI)
      SINPHI=DSIN(DELPHI)
      B0=YBETA(INULLP)
      BETA=YBETA(MBETA)
      ALPHA0=-0.5*YBETAP(INULLP)
      ALPHA=-0.5*YBETAP(MBETA)

      T11=DSQRT(BETA/B0)*(COSPHI+ALPHA0*SINPHI)
      T12=DSQRT(BETA*B0)*SINPHI
      T21=-1.0d0/DSQRT(BETA*B0)*((ALPHA-ALPHA0)*COSPHI+
     &  (1.0d0+ALPHA*ALPHA0)*SINPHI)
      T22=DSQRT(B0/BETA)*(COSPHI-ALPHA*SINPHI)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*) '*** SR BETATRAN ***'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'HORIZONTALE LINEARE TRANSFER-MATRIX OF SECOND HALF OF THE WLS:'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)SNGL(T11),SNGL(T12)
      WRITE(LUNGFO,*)SNGL(T21),SNGL(T22)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'DETERMINANTE:',SNGL(T11*T22-T12*T21)

      RETURN
      END
