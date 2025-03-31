*CMZ :  2.41/10 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.33  by  Michael Scheer
*CMZ : 00.01/09 05/10/95  16.34.35  by  Michael Scheer
*CMZ :  1.00/01 25/09/95  14.26.06  by  Michael Scheer
*CMZ :  1.00/00 21/09/95  17.14.08  by  Michael Scheer
*-- Author :    Michael Scheer   21/09/95
      SUBROUTINE VPOT3DINIT(NDIMP,LORD,MORD,NDORD,INDEX,NCOEF,NSTAKP,FSTAK
     &                        ,LUNCOEFF,FILECOEFF,COMMENT)
*KEEP,gplhint.
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      CHARACTER(60) FILECOEFF,COMMENT

      INTEGER NDIMP,INDEX(4,NDIMP),NSTAKP,NCOEF
     &         ,LUNCOEFF,IORD,LORD,MORD,NDORD

      INTEGER IX,IY,IZ,IS,NS,IR,IXYZ,IFOUND

      DOUBLE PRECISION FSTAK(4,NSTAKP,NDIMP)

      NCOEF=0
      OPEN(UNIT=LUNCOEFF,FILE=FILECOEFF,STATUS='OLD')

          READ(LUNCOEFF,'(A60)')COMMENT
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     SR VPOT3DINIT:'
          WRITE(LUNGFO,*)'     comment on COEF-FILE:'
          WRITE(LUNGFO,*)'     ',COMMENT
          WRITE(LUNGFO,*)
100       READ(LUNCOEFF,*,END=900) IX,IY,IZ,NS

          DO IORD=LORD,MORD,NDORD
         IF(IX+IY+IZ.EQ.IORD) GOTO 50
          ENDDO
          GOTO 100

50        NCOEF=NCOEF+1

          IF (NCOEF.GT.NDIMP) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
     &   '*** ERROR IN VPOT3DINIT: DIMENSION NDIMP EXCEEDED ***'
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
     &   '*** ERROR IN VPOT3DINIT: DIMENSION NDIMP EXCEEDED ***'
             WRITE(6,*)
             STOP
          ENDIF

          IF (NS.GT.NSTAKP) THEN
             WRITE(LUNGFO,*)
             WRITE(LUNGFO,*)
     &   '*** ERROR IN VPOT3DINIT: DIMENSION NSTAKP EXCEEDED ***'
             WRITE(LUNGFO,*)
             WRITE(6,*)
             WRITE(6,*)
     &   '*** ERROR IN VPOT3DINIT: DIMENSION NSTAKP EXCEEDED ***'
             WRITE(6,*)
             STOP
          ENDIF

          INDEX(1,NCOEF)=IX
          INDEX(2,NCOEF)=IY
          INDEX(3,NCOEF)=IZ

          IF (NS.GT.0) THEN
         INDEX(4,NCOEF)=NS
         DO IS=1,NS
             READ(LUNCOEFF,*)(FSTAK(IR,IS,NCOEF),IR=1,4)
         ENDDO
          ELSE
         INDEX(4,NCOEF)=1
         FSTAK(1,1,NCOEF)=IX
         FSTAK(2,1,NCOEF)=IY
         FSTAK(3,1,NCOEF)=IZ
         FSTAK(4,1,NCOEF)=1.D0
          ENDIF   !IS

      GOTO 100

900   CLOSE(LUNCOEFF)


      DO IORD=LORD,MORD,NDORD
      IFOUND=0
      DO IXYZ=1,NCOEF
         IX=INDEX(1,IXYZ)
         IY=INDEX(2,IXYZ)
         IZ=INDEX(3,IXYZ)
         IF (IX+IY+IZ.EQ.IORD) IFOUND=1
      ENDDO
      IF (IFOUND.EQ.0) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)
     &'*** WARNING SR VPOT3DINIT: NO COEFFICIENTS FOUND FOR INQUIRED ORDER ***'
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)'ORDER :',IORD
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)
     &'*** WARNING SR VPOT3DINIT: NO COEFFICIENTS FOUND FOR INQUIRED ORDER ***'
         WRITE(6,*)
         WRITE(6,*)'ORDER :',IORD
         WRITE(6,*)
      ENDIF
      ENDDO

      RETURN
      END
