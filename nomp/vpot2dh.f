*CMZ :  2.41/08 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.44.31  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.00/00 30/06/97  11.38.41  by  Michael Scheer
*CMZ : 00.02/04 10/02/97  14.07.52  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.50.14  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97
      SUBROUTINE VPOT2DH
*KEEP,gplhint.
*KEND.
     &(NPOI,IFAIL)

*KEEP,bpoly2dhf90u.
      include 'bpoly2dhf90u.cmn'
*KEND.

C---  TO FIT POTENTIAL WITH TRANSVERSAL POLYNOMIAL AND
C     LONGITUDINAL SIN/COS-LIKE ANSATZ
C     OF A MAGNETIC FIELD B=(BX,BY,BZ)=-GRAD(V)
C

C--- INPUT:

C     NPOI  : NUMBER OF DATA POINTS X,Y,Z,BX,BY,BZ

C--- OUTPUT:

C     IFAIL : FAILURE FLAG

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly2dh.
      include 'bpoly2dh.cmn'
*KEND.

      INTEGER NPOI


      INTEGER IFAIL,NARG,NFUN,IPOI,IPAR,NPARP
      PARAMETER (NPARP=NPARTOTP,NFUN=3,NARG=3)

      DOUBLE PRECISION WS(NARG+NFUN)
      DOUBLE PRECISION PARAM(NPARP),A(NPARP,NPARP),T(NFUN,NPARP)

C NOTE CHANGE COORDINATE SYSTEMS!

      ALLOCATE(FUNDATA(NARG+NFUN,NPOI))
      DO IPOI=1,NPOI
         fundata(1,IPOI)=-z(IPOI)
         fundata(2,IPOI)=y(IPOI)
         fundata(3,IPOI)=x(IPOI)
         fundata(4,IPOI)=-BZ(IPOI)
         fundata(5,IPOI)=BY(IPOI)
         fundata(6,IPOI)=BX(IPOI)
      ENDDO

      CALL UTIL_LINEAR_FIT
     &  (IFAIL,NPARTOT,PARAM,NPOI,NPOI,NARG,NFUN,A,T,FUNDATA,WS)
      IF (IFAIL.NE.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN VPOT2DH: FIT FAILED'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'IFAIL:',IFAIL
            WRITE(LUNGFO,*) '(maybe no transversal field gradient'
          WRITE(LUNGFO,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN VPOT2DH: FIT FAILED'
          WRITE(6,*)
          WRITE(6,*)'IFAIL:',IFAIL
            WRITE(6,*) '(maybe no transversal field gradient'
          WRITE(6,*)
      ENDIF

      DO IPAR=1,NPARTOT
          PAR2DH(IPAR)=PARAM(IPAR)
      ENDDO

      DEALLOCATE(FUNDATA)
      RETURN
      END
