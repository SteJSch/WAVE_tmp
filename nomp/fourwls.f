*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.16/04 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.26.25  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  17.02.52  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  10.55.14  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.38.33  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.41  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE FOURWLS
*KEEP,gplhint.
*KEND.

C Expands the magnetic field By on the device axis into a Fourier series
C The field has to be symmetric i.e. By(x)= By(-x).
C The field is calculated for the negativ x-values and also assigned
C to the positiv x.

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,fourier.
      include 'fourier.cmn'
*KEND.

      INTEGER I,K,I1,IP,IM,NFOURD2,MFOUR

C--- DUE TO CERN MAIL DIMENSION INCREASED TO 2 (AT LEAST BY 1) 22.06.95

      COMPLEX CKOEF(NFOURD/2+1+2)
      REAL*4  YFOUR(NFOURD+2+2),AKOEF(NFOURD/2+1+2)
      EQUIVALENCE (CKOEF,YFOUR)

      DOUBLE PRECISION XFOUR(NFOURD+2+2)
      DOUBLE PRECISION XLFOUR,DXFOUR,BX,BY,BZ,AX,AY,AZ

C--- LOOP UBER NFOURD-PUNKTE FUER FAST-FOURIER-TRANSFORMATION

      IF (XSTART+XSTOP.NE.0.0) THEN
         WRITE(LUNGFO,*)
         WRITE(LUNGFO,*) '*** WARNING SR FOURWLS: XSTART.NE.-XSTOP ***'
         WRITE(LUNGFO,*) 'CHECK RESULTS CAREFULLY'
         WRITE(LUNGFO,*)
         WRITE(6,*)
         WRITE(6,*) '*** WARNING SR FOURWLS: XSTART.NE.-XSTOP ***'
         WRITE(6,*) 'CHECK RESULTS CAREFULLY'
         WRITE(6,*)
      ENDIF !XSTART

      IF(NFOURWLS.GT.NFOURD/2) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &      '*** ERROR IN FOURWLS ***'
          WRITE(LUNGFO,*)
     &      'NFOURWLS.GT.NFOURD/2'
          WRITE(LUNGFO,*)
     &      'CHECK NFOURWLS OR INCREASE PARAMETER NFOURD IN FILE CMPARA.CMN'
          WRITE(6,*)
          WRITE(6,*)
     &      '*** ERROR IN FOURWLS ***'
          WRITE(6,*)
     &      'NFOURWLS.GT.NFOURD/2'
          WRITE(6,*)
     &      'CHECK NFOURWLS OR INCREASE PARAMETER NFOURD IN FILE CMPARA.CMN'
           STOP
      ENDIF

C     XLFOUR=DABS(XSTART-XSTOP)
      XLFOUR=DABS(2.*XSTOP)   !WEGEN SR BMAGSEQ
      DXFOUR=XLFOUR/NFOURD
      NFOURD2=NFOURD/2
      MFOUR=NINT(ALOG(FLOAT(NFOURD))/ALOG(2.E0))

      DO I=1,NFOURD2+1           !SYMMETRISCHE X-WERTE
          XFOUR(I)          =-DXFOUR*(NFOURD2+1-I)
          XFOUR(NFOURD+1-I+1)=-XFOUR(I)
      END DO

      DO I=1,NFOURD2+1           !SYMMETRISCHE Y-WERTE

          I1=I-1
          IP=NFOURD2+1+I1
          IM=NFOURD2+1-I1

C141091      CALL BREC00(BX,BY,BZ,XFOUR(IP),0.,0.)
          CALL MYBFELD(XFOUR(IP),0.D0,0.D0,BX,BY,BZ,AX,AY,AZ)

          YFOUR(IP)=BY
          YFOUR(IM)=BY

      END DO



      CALL RFFT(CKOEF,-MFOUR) !FFT MIT CERN-ROUTINE D703


      DO K=1,NFOURD2+1  !REELLE KOEFFIZIENTEN
          AKOEF(K)=(-1.)**(K-1)*2.*REAL(CKOEF(K))
      ENDDO

C--- OUTPUT AUF FILE SCHREIBEN

      OPEN(UNIT=LUNF,FILE=FILEF,STATUS='NEW',FORM='FORMATTED')

      WRITE(LUNF,1000)ICODE,CODE
1000    FORMAT(I10,'  ',1A60)
      WRITE(LUNF,*)XLFOUR
      WRITE(LUNF,*)NFOURWLS

      DO I=1,NFOURWLS
          WRITE(LUNF,*)I,AKOEF(I)
      ENDDO

      CLOSE(LUNF)

      IF(IWFILF.NE.99) THEN
      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR FOURWLS: Coefficients written to file:'
      WRITE(LUNGFO,*)'     ',FILEF
      WRITE(LUNGFO,*)

      DO I=1,NFOURWLS
          WRITE(LUNGFO,*)I,AKOEF(I)
      ENDDO
      WRITE(LUNGFO,*)
      ENDIF !IWFILF

      RETURN
      END
