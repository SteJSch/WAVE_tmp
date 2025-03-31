*CMZ :  2.41/10 18/09/2013  12.33.23  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ : 00.01/09 25/10/95  17.57.38  by  Michael Scheer
*-- Author :    Michael Scheer   29/09/95

      SUBROUTINE BPOLY3DINI
*KEEP,gplhint.
*KEND.

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,bpoly3d.
      include 'bpoly3d.cmn'
*KEND.

      INTEGER IX,IY,IZ,NREAD,IORD
      CHARACTER(64) COMMENT,WCOMMENT

      DOUBLE PRECISION CC
     &  ,X3DMIND,Y3DMIND,Z3DMIND
     &  ,X3DMAXD,Y3DMAXD,Z3DMAXD

      IF (LORD3D.LT.1) LORD3D=1

      NREAD=0
      NCIND=0
      OPEN(UNIT=LUN3DFIT,FILE=FILE3DFIT,STATUS='OLD')
          READ(LUN3DFIT,'(A64)')WCOMMENT
          READ(LUN3DFIT,'(A64)')COMMENT
          READ(LUN3DFIT,*)B3DSCALE
          XYZ3DSC=B3DSCALE
          READ(LUN3DFIT,*)X3DMIND,X3DMAXD
          READ(LUN3DFIT,*)Y3DMIND,Y3DMAXD
          READ(LUN3DFIT,*)Z3DMIND,Z3DMAXD
100       READ(LUN3DFIT,*,END=900) IX,IY,IZ,CC
          NREAD=NREAD+1
          IF (IX.LT.1.OR.IX.GT.NDIMC) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)
     &'*** ERROR IN BPOLY3DINI: INDEX ON FILE EXCEEDS DIMENSION ***'
         WRITE(LUNGFO,*)'INDEX IX IS:',IX
         WRITE(LUNGFO,*)'DIMENSION IS:',NDIMC
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)
     &'*** ERROR IN BPOLY3DINI: INDEX ON FILE EXCEEDS DIMENSION ***'
         WRITE(6,*)'INDEX IX IS:',IX
         WRITE(6,*)'DIMENSION IS:',NDIMC
         WRITE(6,*)
         STOP
          ENDIF
          IF (IY.LT.1.OR.IY.GT.NDIMC) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)
     &'*** ERROR IN BPOLY3DINI: INDEX ON FILE EXCEEDS DIMENSION ***'
         WRITE(LUNGFO,*)'INDEX IY IS:',IY
         WRITE(LUNGFO,*)'DIMENSION IS:',NDIMC
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)
     &'*** ERROR IN BPOLY3DINI: INDEX ON FILE EXCEEDS DIMENSION ***'
         WRITE(6,*)'INDEX IY IS:',IY
         WRITE(6,*)'DIMENSION IS:',NDIMC
         WRITE(6,*)
         STOP
          ENDIF
          IF (IZ.LT.1.OR.IZ.GT.NDIMC) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*)
     &'*** ERROR IN BPOLY3DINI: INDEX ON FILE EXCEEDS DIMENSION ***'
         WRITE(LUNGFO,*)'INDEX IZ IS:',IZ
         WRITE(LUNGFO,*)'DIMENSION IS:',NDIMC
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*)
     &'*** ERROR IN BPOLY3DINI: INDEX ON FILE EXCEEDS DIMENSION ***'
         WRITE(6,*)'INDEX IZ IS:',IZ
         WRITE(6,*)'DIMENSION IS:',NDIMC
         WRITE(6,*)
         STOP
          ENDIF
          C(IX,IY,IZ)=0.D0
          IF (
     &               IX-1+IY-1+IZ-1.GE.LORD3D
     &         .AND.
     &               IX-1+IY-1+IZ-1.LE.MORD3D
     &      ) THEN
             C(IX,IY,IZ)=CC
             DO IORD=LORD3D,MORD3D,NDORD3D
                IF (IX+IY+IZ.EQ.IORD+3) THEN
                   NCIND=NCIND+1
                   IF (NCIND.GT.NDIMCIP) THEN
                     WRITE(LUNGFO,*)
                     WRITE(LUNGFO,*)
     &'*** ERROR IN BPOLY3DINI: DIMENSION NDIMCIP (COMMON BPOLY3DC) EXCEEDED ***'
                     WRITE(LUNGFO,*)
                     WRITE(6,*)
                     WRITE(6,*)
     &'*** ERROR IN BPOLY3DINI: DIMENSION NDIMCIP (COMMON BPOLY3DC) EXCEEDED ***'
                     WRITE(6,*)
                     STOP '--- PROGRAM ABORTED ---'
                   ENDIF
                   ICIND(1,NCIND)=IX
                   ICIND(2,NCIND)=IY
                   ICIND(3,NCIND)=IZ
                   CIND(NCIND)=CC
                ENDIF
             ENDDO
          ENDIF
          GOTO 100
900   CLOSE(LUN3DFIT)

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR BPOLY3DINI:'
      WRITE(LUNGFO,*)'     Number of coefficients read:',NREAD
      WRITE(LUNGFO,*)'     from file:'
      WRITE(LUNGFO,*)'     ',FILE3DFIT
      WRITE(LUNGFO,*)'     WAVE comment on file'
      WRITE(LUNGFO,*)'     ',WCOMMENT
      WRITE(LUNGFO,*)'     comment on file'
      WRITE(LUNGFO,*)'     ',COMMENT
      WRITE(LUNGFO,*)'     number of coefficients used:',NCIND
      WRITE(LUNGFO,*)'     scaling factor:',B3DSCALE
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     x-range of fit:',X3DMIND,X3DMAXD
      WRITE(LUNGFO,*)'     y-range of fit:',Y3DMIND,Y3DMAXD
      WRITE(LUNGFO,*)'     z-range of fit:',Z3DMIND,Z3DMAXD
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)

      RETURN
      END
