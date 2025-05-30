*CMZ :  3.00/01 20/03/2013  10.20.10  by  Michael Scheer
*CMZ :  2.44/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.42/04 29/10/2002  10.25.44  by  Michael Scheer
*CMZ :  2.41/08 12/08/2002  16.38.26  by  Michael Scheer
*CMZ :  2.16/08 29/10/2000  17.38.47  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ :  1.00/00 30/06/97  11.45.36  by  Michael Scheer
*CMZ : 00.02/04 18/02/97  12.21.58  by  Michael Scheer
*CMZ : 00.02/03 04/02/97  16.50.14  by  Michael Scheer
*-- Author :    Michael Scheer   22/01/97
      SUBROUTINE VPOTPH(NPOI,IFAIL)
*KEEP,gplhint.
*KEND.

*KEEP,bpharmf90u.
      include 'bpharmf90u.cmn'
*KEND.

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
*KEEP,bpharm.
      include 'bpharm.cmn'
*KEND.

      INTEGER NPOI

      INTEGER IFAIL,NARG,NFUN,IPOI,ITRANS,IHARM
      PARAMETER (NFUN=3,NARG=3)

      DOUBLE PRECISION WS(NARG+NFUN)
      DOUBLE PRECISION A(NPARPHP,NPARPHP),T(NFUN,NPARPHP)

      ALLOCATE(FUNDATA(NARG+NFUN,NPOI))

C NOTE CHANGE OF COORDINATE SYSTEMS!

      DO IPOI=1,NPOI
        fundata(1,IPOI)=-z(IPOI)
        fundata(2,IPOI)=y(IPOI)
        fundata(3,IPOI)=x(IPOI)
        fundata(4,IPOI)=-BZ(IPOI)
        fundata(5,IPOI)=BY(IPOI)
        fundata(6,IPOI)=BX(IPOI)
      ENDDO

      NPARPH=3
      DO ITRANS=NTRANS0,NTRANS,NTRANSD
        IF (NTRANS0.GT.0) THEN
          DO IHARM=NHARM0,NHARM,NHARMD
            NPARPH=NPARPH+2
          ENDDO
        ENDIF
      ENDDO

      CALL BPHARM_FIT
     &  (IFAIL,NPARPH,PARPH,NPOI,NPOI,NARG,NFUN,A,T,FUNDATA,WS)

      IF (IFAIL.NE.0) THEN
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'*** WARNING IN VPOTPH: FIT FAILED'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)'IFAIL:',IFAIL
        WRITE(LUNGFO,*) '(maybe no transversal field gradient'
        WRITE(LUNGFO,*) '(or XLENCPH, YLENSPH zero)'
        WRITE(LUNGFO,*)
        WRITE(6,*)
        WRITE(6,*)'*** WARNING IN VPOTPH: FIT FAILED'
        WRITE(6,*)
        WRITE(6,*)'IFAIL:',IFAIL
        WRITE(6,*) '(maybe no transversal field gradient'
        WRITE(6,*) '(or XLENCPH, YLENSPH zero)'
        WRITE(6,*)
      ENDIF

      DEALLOCATE(FUNDATA)

      RETURN
      END
