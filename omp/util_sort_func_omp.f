*CMZ :  3.07/00 05/03/2019  21.43.24  by  Michael Scheer
*CMZ :  3.05/05 12/07/2018  13.19.00  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.00/00 06/06/97  16.44.06  by  Michael Scheer
*CMZ : 00.01/07 08/03/95  10.05.45  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SORT_FUNC_omp(N,RA,YA)
*KEEP,gplhint.
*KEND.

C--- HEAPSORT ROUTINE; SEE NUMERICAL RECIPES 8.2 S 231
C--- ARRAY YA IS FUNCTION OF RA AND SORTED ACCORDINGLY

      IMPLICIT NONE

      INTEGER N,L,IR,I,J

      DOUBLE PRECISION RA(N),RRA
      DOUBLE PRECISION YA(N),YYA

      IF (N.LT.2) RETURN

      L=N/2+1
      IR=N

10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          YYA=YA(L)
        ELSE
          RRA=RA(IR)
          YYA=YA(IR)
          RA(IR)=RA(1)
          YA(IR)=YA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            YA(1)=YYA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            YA(I)=YA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        YA(I)=YYA
      GO TO 10
      END
