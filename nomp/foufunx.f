*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.15/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.14/02 27/04/2000  17.56.53  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.51.21  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.19  by  Michael Scheer
*-- Author : Michael Scheer
C**************************************************************
      FUNCTION FOUFUNX(X,WSIG)
C*************************************************************
*KEEP,gplhint.
*KEND.
      IMPLICIT NONE
      DOUBLE PRECISION FOUFUNX,X,WSIG
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      FOUFUNX=DEXP(-X*X/(WSIG*WSIG)/2.D0)/WSIG/DSQRT(2.D0*PI1)

      RETURN
      END
