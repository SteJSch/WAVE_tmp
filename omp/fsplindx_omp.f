*CMZ :  4.00/16 10/09/2022  09.39.11  by  Michael Scheer
*CMZ :  3.07/00 06/03/2019  13.33.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  07.43.29  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.47/07 14/04/2003  15.17.05  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  16.33.46  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.36.50  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.55  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE FSPLINDX_omp(DX,Y,N,YP1,YPN,Y2)
*KEEP,gplhint.
*KEND.


      IMPLICIT NONE
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.

      INTEGER N,J
      DOUBLE PRECISION DX,Y(N),Y2(N)


      INTEGER max
      DOUBLE PRECISION  C(max(NDOBSVZP,NDOBSVYP))
      DOUBLE PRECISION AA(max(NDOBSVZP,NDOBSVYP))
      DOUBLE PRECISION BB(max(NDOBSVZP,NDOBSVYP))
      DOUBLE PRECISION CC(max(NDOBSVZP,NDOBSVYP))

      DOUBLE PRECISION YP1,YPN,H6,H3

      Y2(1)=YP1
      Y2(N)=YPN

      IF (N.LT.3) RETURN

      C(1)=YP1
      C(N)=YPN

      H6=1./6.D0
      H3=4.D0*H6

      BB(1)=1.D0
      CC(1)=0.D0


      DO J=2,N-1
          AA(J)=H6
          BB(J)=H3
          CC(J)=H6
          C(J)=(Y(J+1)-2.D0*Y(J)+Y(J-1))/DX**2
      ENDDO !J

      DO J=2,N-1

          BB(J)=BB(J)-AA(J)*CC(J-1)
           C(J)= C(J)-AA(J)* C(J-1)
C030414          AA(J)=AA(J)-AA(J)*BB(J-1)

          CC(J)=CC(J)/BB(J)
           C(J)= C(J)/BB(J)
          BB(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
         Y2(J)=C(J)-CC(J)*Y2(J+1)
      ENDDO

      RETURN
      END
