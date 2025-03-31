*CMZ :  3.05/01 04/05/2018  19.58.00  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.10.30  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.50/00 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.20/01 29/11/2000  18.31.06  by  Michael Scheer
*CMZ :  2.16/08 01/11/2000  18.41.44  by  Michael Scheer
*CMZ :  2.15/00 04/05/2000  17.02.01  by  Michael Scheer
*CMZ :  2.12/00 27/05/99  10.26.26  by  Michael Scheer
*CMZ :  2.10/01 24/02/99  10.20.40  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.17.32  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.26  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.11.48  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE WFILINT
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,reargf90u.
      include 'reargf90u.cmn'
*KEND.

C--- DUMP INTEGRAND

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEND.

      INTEGER IPOI,IC,ICAL


      COMPLEX*16 AX,AY,AZ,EXPOM
      DOUBLE PRECISION OM

      DATA ICAL/0/

      OM=FREQ(1)/HBAREV1

C-- LOOP OVER TIME STEPS (ACTUAL INTEGRATION)

      IF (ICAL.EQ.0) THEN  !CV2
C         OPEN(UNIT=LUNINT,FILE=FILEINT,STATUS='NEW')
        ICAL=1
      ENDIF !ICAL     CV2

      DO IPOI=1,IARGUM

        IF (WSOU(1,1,IPOI).GE.XIANF.AND.WSOU(1,1,IPOI).LE.XIEND) THEN

          EXPOM=CDEXP(DCMPLX(0.D0,REARGUM(4,IPOI)*OM))

          AX=DCMPLX(REARGUM(1,IPOI))*EXPOM
          AY=DCMPLX(REARGUM(2,IPOI))*EXPOM
          AZ=DCMPLX(REARGUM(3,IPOI))*EXPOM

          WRITE(LUNINT,*) WSOU(1,1,IPOI)
     &      ,(SNGL(REARGUM(IC,IPOI)),IC=1,3)
     &      ,SNGL(REARGUM(4,IPOI)*OM),SNGL(REARGUM(5,IPOI))
     &      ,SNGL(DREAL(EXPOM)),SNGL(DIMAG(EXPOM))
     &      ,SNGL(DREAL(AX)),SNGL(DIMAG(AX))
     &      ,SNGL(DREAL(AY)),SNGL(DIMAG(AY))
     &      ,SNGL(DREAL(AZ)),SNGL(DIMAG(AZ))

        ENDIF   !XIANF

      ENDDO   !LOOP OVER TIME STEPS

C      IF (NSADD.NE.0) CLOSE(LUNINT)

CV2   CLOSE(LUNINT)

      RETURN
      END
