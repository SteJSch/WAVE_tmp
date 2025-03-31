*CMZ :  4.00/11 07/05/2021  07.34.04  by  Michael Scheer
*CMZ :  4.00/07 09/01/2020  13.45.07  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.29.14  by  Michael Scheer
*CMZ :  2.14/02 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.13/09 09/03/2000  12.40.42  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.09.58  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.56.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.51.02  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.55  by  Michael Scheer
*-- Author : Chaoen Wang
C...............................................................................
      Subroutine FFT(Z1,M,ISign)
C...............................................................................
*KEEP,gplhint.
*KEND.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      COMPLEX Z1(NDOBSVZP+NDOBSVYP) !SR USMCON2

CMSH  Parameter (MaxPt0=616)
CMSH  CompLex Z1(MaxPt0)

      CompLex ZU,ZW,ZT

      INTEGER K,NM1,NV2,I,J,N,IP,LE,LE1,L,M,ISIGN

      REAL PI

      DATA PI/3.141592653589793D0/

      print*,"FFT IS OBSOLETE, SEE //WAVE/USEM IN WAVE.CMZ"

      Return
      End
