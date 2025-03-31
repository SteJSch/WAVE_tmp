*CMZ :  3.07/00 07/03/2019  10.08.17  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  07.43.29  by  Michael Scheer
*CMZ :  2.48/04 12/03/2004  15.40.31  by  Michael Scheer
*CMZ :  2.47/07 14/04/2003  15.17.05  by  Michael Scheer
*CMZ :  2.37/02 14/11/2001  12.53.09  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  1.03/06 06/08/98  18.15.25  by  Michael Scheer
*CMZ : 00.01/04 28/11/94  18.44.24  by  Michael Scheer
*CMZ : 00.01/02 18/11/94  16.38.52  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.51.45  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.56  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE FSPLPER_omp(DX,Y,N,Y2)
*KEEP,gplhint.
*KEND.

C--- FOR PERIODIC CASE
C--- ONLY FOR EQUIDISTANT X-VALUES

C ROUTINES SOLVES THE FOLLOWING SYSTEM OF LINEAR EQUATIONS:

C  BB( 1 ) CC( 1 )   0     ....... ....... .......   0     ZZ( 1 ) | C( 1 )
C  AA( 2 ) BB( 2 ) CC( 2 )    0    ....... .......   0        0    | C( 2 )
C  ....... AA( 3 ) BB( 3 ) CC( 3 )    0    .......   0        0    | C( 3 )
C              .
C              .
C              .
C     0    ....... .......    0    AA(N-3) BB(N-3) CC(N-3)    0    | C(N-3)
C     0    ....... ....... .......    0    AA(N-2) BB(N-2) ZZ(N-2) | C(N-2)
C  CC(N-1)    0    ....... ....... .......    0    AA(N-1) ZZ(N-1) | C(N-1)


      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEND.

      INTEGER N,J
      DOUBLE PRECISION DX,Y(N),Y2(N)

      DOUBLE PRECISION  C(2*(NDOBSVZP+NDOBSVYP))
      DOUBLE PRECISION AA(2*(NDOBSVZP+NDOBSVYP))
      DOUBLE PRECISION BB(2*(NDOBSVZP+NDOBSVYP))
      DOUBLE PRECISION CC(2*(NDOBSVZP+NDOBSVYP))
      DOUBLE PRECISION ZZ(2*(NDOBSVZP+NDOBSVYP))


      DOUBLE PRECISION H6,H3

      IF (N.LT.4) THEN
        WRITE(6,*) '*** ERROR IN FSPLPER_omp: N TOO LOW ***'
        STOP
      ENDIF

      H6=DX/6.D0
      H3=4.D0*H6

      BB(1)=1.0
      CC(1)=0.25
      ZZ(1)=0.25
      C(1)=( Y(3)-2.D0*Y(2)+Y(1) )/DX/H3

      DO J=2,N-3
          AA(J)=H6
          BB(J)=H3
          CC(J)=H6
          ZZ(J)=0.D0
          C(J)=(Y(J+1+1)-2.D0*Y(J+1)+Y(J-1+1))/DX
      ENDDO !J

      AA(N-2)=H6
      BB(N-2)=H3
      ZZ(N-2)=H6
      C (N-2)=(Y(N-2+1+1)-2.D0*Y(N-2+1)+Y(N-2-1+1))/DX

      AA(N-1)=H6
      ZZ(N-1)=H3
      CC(N-1)=H6
      C (N-1)=(Y(2)-2.D0*Y(N)+Y(N-1))/DX

      DO J=2,N-3

          ZZ(J)=ZZ(J)-AA(J)*ZZ(J-1)
           C(J)= C(J)-AA(J)* C(J-1)
          BB(J)=BB(J)-AA(J)*CC(J-1)
C030414          AA(J)=AA(J)-AA(J)*BB(J-1)

          CC(J)=CC(J)/BB(J)
          ZZ(J)=ZZ(J)/BB(J)
           C(J)= C(J)/BB(J)
          BB(J)=1.D0

      ENDDO !J

      BB(N-2)=BB(N-2)-AA(N-2)*CC(N-3)
      ZZ(N-2)=ZZ(N-2)-AA(N-2)*ZZ(N-3)
      C (N-2)=C(N-2) -AA(N-2)*C (N-3)

      ZZ(N-2)=ZZ(N-2)/BB(N-2)
      C (N-2)=C (N-2)/BB(N-2)
      BB(N-2)=1.D0

      ZZ(N-1)=ZZ(N-1)-CC(N-1)*ZZ(1)
       C(N-1)= C(N-1)-CC(N-1)* C(1)
      CC(N-1)=       -CC(N-1)*CC(1) ! CC MOVES DOWN THE ROW

      DO J=2,N-4
        ZZ(N-1)=ZZ(N-1)-CC(N-1)*ZZ(J)
        C(N-1)=  C(N-1)-CC(N-1)* C(J)
        CC(N-1)=       -CC(N-1)*CC(J)  ! CC MOVES DOWN THE ROW
      ENDDO !J

      AA(N-1)=AA(N-1)-CC(N-1)*CC(N-3)
      ZZ(N-1)=ZZ(N-1)-CC(N-1)*ZZ(N-3)
      C(N-1)=  C(N-1)-CC(N-1)*C(N-3)

      ZZ(N-1)=ZZ(N-1)-AA(N-1)*ZZ(N-2)
       C(N-1)= C(N-1)-AA(N-1)* C(N-2)
      AA(N-1)=AA(N-1)-AA(N-1)*BB(N-2)

       C(N-1)=C(N-1)/ZZ(N-1)
      ZZ(N-1)=1.D0

      Y2(N-1+1)= C(N-1)
      CC(N-2  )=ZZ(N-2)

      Y2(N-2+1)=C(N-2)-CC(N-2)*Y2(N-1+1)
      DO J=N-3,1,-1
         Y2(J+1)=C(J)-CC(J)*Y2(J+1+1)-ZZ(J)*Y2(N-1+1)
      ENDDO

      Y2(1)=Y2(N)

      RETURN
      END
