*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.36.44  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.45  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE EMINP(R,B3,T0,TK,EMIN)
*KEEP,gplhint.
*KEND.

C     SOLVES EQUATION TO DETERMINE MINIMUM ENERGY FOR GIVEN POLARIZATION LEVEL
C     LOGBUCH SEITE 189

          IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

          DOUBLE PRECISION R,B3,T0,TK,EMIN,E
          DOUBLE PRECISION A,B,EN,DEFUNP,EFUNP
          INTEGER I

C--- ALGORITHM OF NEWTON X(N+1)=X(N)-F(X(N))/F'(X(N))

      E=1.D0
      DO I=1,100
          A=DEFUNP(R,B3,T0,TK,E)
          B=EFUNP(R,B3,T0,TK,E)-A*E
          EN=-B/A
          IF(DABS((EN-E)/EN).LT.1.D-10) GOTO 100
          E=EN
          END DO

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'*** WARNING SR EMINP ***'
      WRITE(LUNGFO,*)'NEWTON ALGORITHM FAILED TO DETERMINE MINIMUM ENERGY'
      WRITE(LUNGFO,*)'CHECK RESULTS CAREFULLY'
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

100   CONTINUE

      EMIN=EN

      RETURN
      END
