*CMZ :  3.08/01 01/04/2019  16.57.05  by  Michael Scheer
*CMZ :  3.07/00 15/03/2019  12.57.36  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.10  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.68/00 25/05/2012  16.23.18  by  Michael Scheer
*CMZ :  2.35/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.34/09 26/09/2001  11.56.11  by  Michael Scheer
*CMZ :  2.34/00 11/05/2001  12.20.31  by  Michael Scheer
*CMZ :  2.16/08 24/10/2000  14.15.55  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.08.20  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  14.27.55  by  Michael Scheer
*CMZ :  2.11/01 18/05/99  11.15.36  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  15.04.42  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.13.09  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.50  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.03  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BLENSTO(kfreq)
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- INTEGRATES THE STOKES INTENSITIES

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DOUBLE PRECISION DSUM,RPHI,SUMP

      INTEGER kfreq,IY,IZ,IOBSV
     &  ,ICAL,IR,IP,NR,MR,MP,ISTOK,II,KDUM,IERR

      DOUBLE PRECISION SUMZ(NDOBSVYP),SUM,DIA
      DOUBLE PRECISION SUMY(NDOBSVZP)
      DOUBLE PRECISION R(NDOBSVZP),PHI(NDOBSVYP)
     &  ,FPHI(NDOBSVYP)
     &  ,SZ(NDOBSVZP),SY(NDOBSVYP)
     &  ,SZY(NDOBSVZP,NDOBSVYP)
     &  ,FZ(NDOBSVZP),FY(NDOBSVYP),DPHI,DR,X,Y

      data ical/0/

      if (ipin.eq.3) then
        DO ISTOK=1,4
          if (ipincirc.eq.0) then
            wstokes(istok,kfreq)=stokes(istok,kfreq)*pinh*pinw
          else
            wstokes(istok,kfreq)=stokes(istok,kfreq)*pinr**2*pi1
          endif
        enddo
        return
      endif

      DO ISTOK=1,4

        IF (IPINCIRC.EQ.0) THEN

          IF (IF1DIM.NE.2) THEN

            DO IY=1,NOBSVY

              SUMZ(IY)=0.0d0

              IF(MOBSVZ.GT.1) THEN

                DO IZ=1,NOBSVZ
                  FZ(IZ)=STOKES(ISTOK,(IY-1)*NOBSVZ+IZ+NOBSV*(kfreq-1))
                ENDDO      !IZ

                CALL FSPLINDX(OBSVDZ,FZ,NOBSVZ,0.D0,0.D0,SZ)

                DO IZ=(NOBSVZ-MOBSVZ)/2+1,(NOBSVZ-MOBSVZ)/2+MOBSVZ-1

                  DSUM=OBSVDZ*0.5D0
     &              *(FZ(IZ)+FZ(IZ+1))
     &              -OBSVDZ**3/24.D0*(SZ(IZ)+SZ(IZ+1))

                  SUMZ(IY)=SUMZ(IY)+DSUM

                ENDDO   !IZ

              ELSE  !MOBSVZ

                IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
                SUMZ(IY)=OBSVDZ*STOKES(ISTOK,IOBSV+NOBSV*(kfreq-1))
              ENDIF !MOBSVZ

            ENDDO !IY

          ELSE  !IF1DIM.EQ.2

            DO IY=1,NOBSVY

              SUMZ(IY)=0.0d0

              IOBSV=(IY-1)*NOBSVZ+(NOBSVZ-MOBSVZ)/2+1
              DIA=ABS((PINR-(PINCEN(2)-OBSV(2,IOBSV)))
     &          *(PINR+(PINCEN(2)-OBSV(2,IOBSV))))
              IF (DIA.GT.0.D0) THEN
                DIA=2.D0*SQRT(DIA)
              ELSE
                DIA=0.0D0
              ENDIF
              SUMZ(IY)=DIA*STOKES(ISTOK,IOBSV+NOBSV*(kfreq-1))

            ENDDO !IY

          ENDIF !IF1DIM.EQ.2

C--- INTEGRATION ALONG VERTICAL AXIS Y

          CALL FSPLINDX(OBSVDY,SUMZ,NOBSVY,0.D0,0.D0,SY)

          IF(MOBSVY.GT.1) THEN

            SUM=0.0d0
            DO IY=(NOBSVY-MOBSVY)/2+1,(NOBSVY-MOBSVY)/2+MOBSVY-1

              DSUM=
     &          OBSVDY*0.5D0
     &          *(SUMZ(IY)+SUMZ(IY+1))
     &          -OBSVDY**3/24.D0
     &          *(SY(IY)+SY(IY+1))

              SUM=SUM+DSUM

            ENDDO

          ELSE IF (IF1DIM.EQ.2) THEN

            SUM=PI1*PINR*SUMZ(NOBSVY/2+1)/2.D0

          ELSE

            SUM=OBSVDY*SUMZ(NOBSVY/2+1)

          ENDIF

        ELSE  !PINCIRC

          IF (IRPHI.NE.0) THEN !INTEGRATION WITH RESPECT TO POLAR COORDINATES

C--- INTEGRATION OVER PHI

            IF (ICAL.EQ.0) THEN

              DR=DMIN1(OBSVDZ,OBSVDY)
              MR=NINT(PINR/DR)+1
              DR=PINR/(MR-1)
              MEDGER=MIN( MEDGEZ, MEDGEY)
              NR=MR+MEDGER
              MP=MOBSVY

              IF (MR.LT.1) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** ERROR IN BLENSTO ***'
                WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
                WRITE(LUNGFO,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
                WRITE(LUNGFO,*)
                WRITE(6,*)
                WRITE(6,*)'*** ERROR IN BLENSTO ***'
                WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
                WRITE(6,*)'INCREASE PARAMETER MPINZ IN NAMELIST PINHOLE'
                WRITE(6,*)
              ENDIF
              IF (MP.LT.4) THEN
                WRITE(LUNGFO,*)
                WRITE(LUNGFO,*)'*** ERROR IN BLENSTO ***'
                WRITE(LUNGFO,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
                WRITE(LUNGFO,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
                WRITE(LUNGFO,*)
                WRITE(6,*)
                WRITE(6,*)'*** ERROR IN BLENSTO ***'
                WRITE(6,*)'NOT ENOUGH GRID POINTS FOR CIRCULAR PINHOLE'
                WRITE(6,*)'INCREASE PARAMETER MPINY IN NAMELIST PINHOLE'
                WRITE(6,*)
              ENDIF

              DPHI=2.D0*PI1/(MP-1)
              DO IP=1,MP
                PHI(IP)=(IP-1)*DPHI
              ENDDO

              DO IR=1,NR
                R(IR)=(IR-1)*DR
              ENDDO

              DO IR=2,NR
                DO IP=1,MP
                  XC(IP+(IR-1)*NOBSVY)=R(IR)*DCOS(PHI(IP))+PINCEN(3)
                  YC(IP+(IR-1)*NOBSVY)=R(IR)*DSIN(PHI(IP))+PINCEN(2)
                ENDDO !IP
              ENDDO !IR

              ICAL=1

            ENDIF !ICAL

C--- INTERPOLATION OF INTENSITY FOR CIRCULAR GRID

            DO IY=1,NOBSVY
              DO IZ=1,NOBSVZ
                II=(IY-1)*NOBSVZ+IZ
                FZ(IZ)=STOKES(ISTOK,II+NOBSV*(kfreq-1))
              ENDDO !IZ

              CALL FSPLINDX(OBSVDZ,FZ,NOBSVZ,0.D0,0.D0,SZ)

              DO IZ=1,NOBSVZ
                SZY(IZ,IY)=SZ(IZ)
              ENDDO !IZ

            ENDDO !IY

            DO IR=2,NR
              DO IP=1,MP

                X=XC(IP+(IR-1)*NOBSVY)
                Y=YC(IP+(IR-1)*NOBSVY)

                DO IY=1,NOBSVY

                  DO IZ=1,NOBSVZ
                    FZ(IZ)=STOKES(ISTOK,(IY-1)*NOBSVZ+IZ+NOBSV*(kfreq-1))
                    SZ(IZ)=SZY(IZ,IY)
                  ENDDO      !IZ

                  CALL SPLINZY(NOBSVZ,X,FY(IY),OBSVZ,FZ,SZ,KDUM)

                ENDDO !IY

                CALL FSPLINDX(OBSVDY,FY,NOBSVY,0.D0,0.D0,SY)
                CALL SPLINZY(NOBSVY,Y,FPHIR(IP+(IR-1)*NOBSVY),OBSVY,FY,SY,KDUM)

              ENDDO !IP
            ENDDO !IR

C--- DO THE INTEGRATION OF FPHIR OVER PHI AND R

            SUM=0.0D0
            SUMY(1)=0.0D0
            DO IR=2,NR  !FIRST RADIUS IS ZERO

              DO IP=1,MP
                FPHI(IP)=FPHIR(IP+(IR-1)*NOBSVY)
              ENDDO   !IP

              CALL FSPLPER(DPHI,FPHI,MP,SY)

              SUMY(IR)=0.0D0
              RPHI=R(IR)*DPHI
              DO IP=1,MP-1

                DSUM=
     &            RPHI*0.5D0*(FPHI(IP)+FPHI(IP+1))
     &            -RPHI**3/24.D0*(SY(IP)+SY(IP+1))

                SUMY(IR)=SUMY(IR)+DSUM

              ENDDO   !IP

            ENDDO !IR

            CALL FSPLINDX(DR,SUMY,NR,0.D0,0.D0,SZ)

            SUM=0.0d0
            DO IR=1,MR-1

              DSUM=DR*0.5D0
     &          *(SUMY(IR)+SUMY(IR+1))
     &          -DR**3/24.D0
     &          *(SZ(IR)+SZ(IR+1))

              SUM=SUM+DSUM

            ENDDO

          ELSE  !IRPHI

            DO IOBSV=1,NOBSV
              FPHIR(IOBSV)=STOKES(ISTOK,IOBSV+NOBSV*(kfreq-1))
            ENDDO !IOBSV
            CALL CIRCPIN(NOBSVZ,NOBSVY,MOBSVZ,MOBSVY,FPHIR,SUM,SUMP,-ISTOK,kfreq,IERR)
          ENDIF !IRPHI

        ENDIF !PINCIRC

        WSTOKES(ISTOK,kfreq)=SUM

      ENDDO !ISTO

      RETURN
      END
