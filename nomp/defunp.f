*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.15/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.53.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.37.11  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.45  by  Michael Scheer
*-- Author : Michael Scheer
*KEEP,gplhint.
*KEND.
C******************************************************************

      DOUBLE PRECISION FUNCTION DEFUNP(R,B3,T0,TK,E)

          IMPLICIT NONE
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.
          DOUBLE PRECISION R,B3,T0,TK,E

      T0=T0
      TK=TK
          DEFUNP=5.*E**4+2.*(CLIGHT1*1.D-9)**3*R**2/(2.*PI1)*B3*E

          RETURN
      END
