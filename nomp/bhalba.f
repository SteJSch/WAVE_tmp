*CMZ :  4.00/15 07/04/2022  08.20.53  by  Michael Scheer
*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.67/04 11/05/2012  11.18.26  by  Michael Scheer
*CMZ :  2.64/01 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.36/00 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/10 21/08/96  12.30.33  by  Michael Scheer
*CMZ : 00.01/04 30/11/94  14.09.41  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.47.18  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.57  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BHALBA(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

*KEEP,gplhint.
*KEND.

C Calculates magnetic field according to Halbach's formulas
C Inside the routine the coodinate system corresponds to
C Halbach's convention i.e. Z is the longitudinal coordinate
C Input and Output are converted according to the coordinate
C system of WAVE

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,halbach.
      include 'halbach.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL

      DOUBLE PRECISION XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT,
     &  XKX,YKY,ZKZ,DSNXKX,DCSXKX,DSHYKY,DCHYKY,DSNZKZ,DCSZKZ
     &  ,BXH,BYH,BZH,AXH,AYH,AZH

      DOUBLE PRECISION PARK,WLEN1,EHARM1

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN

        park=pkhalba

        if (nhhalba.ne.0.and.hhalba.ne.0.0d0) then
          if (hhalba.eq.-9999.0d0) then
            if (ifreq2p.eq.1) then
              hhalba=freqlow
            else
              hhalba=(freqlow+freqhig)/2.0d0
            endif
          else if (hhalba.lt.0.0d0) then
            hhalba=-wtoe1/hhalba
          endif
          WLEN1=wtoe1/abs(hhalba/nhhalba)
          park=2.0d0*(wlen1/(zlhalba*1.0D9/2.0d0/DMYGAMMA**2)-1.0d0)
          if (park.lt.0.0d0) then
            write(6,*)
     &        '*** Error in BHALBA:'
            write(6,*)
     &        'Inconsistent values of NHHALBA, HHALBA, and ZLHALBA'
            write(6,*)' '
            write(lungfo,*)
     &        '*** Error in BHALBA:'
            write(lungfo,*)
     &        'Inconsistent values of NHHALBA, HHALBA, and ZLHALBA'
            write(lungfo,*)' '
            stop
          endif
          park=sqrt(park)
          pkhalba=park
        endif

        IF (park.ne.0.0d0) THEN
          B0halba=pkhalba/(echarge1*zlhalba/(2.*PI1*EMASSKG1*CLIGHT1))
        endif

        IF(B0HALBA.EQ.9999.0D0) THEN
          OPEN(UNIT=LUNH,FILE=FILEH,STATUS='OLD',FORM='FORMATTED')
          READ(LUNH,*)B0HALBA,XLHALBA,YLHALBA,ZLHALBA
          ZLENHAL=ZLHALBA
        ENDIF

        IF (PERHAL.LT.0.) THEN
          ZLHALBA=-ZLENHAL/PERHAL
          XLHALBA= XLHALBA/PERHAL
        ENDIF

        XKHALBA=0.0D0
        YKHALBA=0.0D0
        ZKHALBA=0.0D0

        IF (ZLHALBA.NE.0.0D0) ZKHALBA=2.0D0*PI1/ZLHALBA
        IF (YLHALBA.NE.0.0D0) YKHALBA=2.0D0*PI1/YLHALBA
        IF (XLHALBA.NE.0.0D0) XKHALBA=2.0D0*PI1/XLHALBA

C--- ADJUST K-VALUES

        YKHALBA=DSQRT(ZKHALBA*ZKHALBA+XKHALBA*XKHALBA)
        YLHALBA=2.0D0*PI1/YKHALBA

        IF(PERHAL.GE.0) ZLENHAL=PERHAL*ZLHALBA     !29.10.91

        ZLENHAL2=ZLENHAL/2.0D0

C--- BENDING RADIUS

C21.8.96         RHALBA=DMYGAMMA*EMASSE1/(CLIGHT1*B0HALBA)
        IF (B0HALBA.NE.0.0D0) THEN
          RHALBA=EMOM/(CLIGHT1*B0HALBA)
          COHAL=1.0D0/(ZKHALBA**2*RHALBA)   !CLOSED ORBIT OFFSET
          zamphalba=b0halba*clight1/emom/zkhalba**2
        ELSE
          RHALBA=0.0D0
          COHAL=0.0D0   !CLOSED ORBIT OFFSET
        ENDIF  !(B0HALBA.NE.0.0D0)

        PARK=ECHARGE1*DABS(B0HALBA)*ZLHALBA/(2.*PI1*EMASSKG1*CLIGHT1)
        WLEN1=(1+PARK**2/2.)/2./DMYGAMMA**2*ZLHALBA*1.0D9
        IF (WLEN1.NE.0.0) EHARM1=WTOE1/WLEN1

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      SR BHALBA: parameters of Halbach wiggler'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '      B0 [T] and corresponding radius [m]:     ',
     &    SNGL(B0HALBA),SNGL(RHALBA)
        WRITE(LUNGFO,*)
     &    '      Length and half length of device [m]:    '
     &    ,SNGL(ZLENHAL),SNGL(ZLENHAL2)
        WRITE(LUNGFO,*)
     &    '      Approximated path amplitude for zero gradient [m]:',
     &    SNGL(zamphalba)
        WRITE(LUNGFO,*)'      L0, L0X, L0Y [m]:'
     &    ,SNGL(ZLHALBA),SNGL(XLHALBA),SNGL(YLHALBA)
        WRITE(LUNGFO,*)'      K, KX, KY [1/m]: '
     &    ,SNGL(ZKHALBA),SNGL(XKHALBA),SNGL(YKHALBA)
        WRITE(LUNGFO,*)
     &    '      Number of periods, closed orbit offset [m]:'
     &    ,SNGL(PERHAL), SNGL(COHAL)
        WRITE(LUNGFO,*)
     &    '      Deflection parameter K, 1. harmonical [eV]:'
     &    ,SNGL(PARK),SNGL(EHARM1)
        WRITE(LUNGFO,*)'      Critical energy [eV]:'
     &    ,SNGL(ecdipev1*DABS(B0HALBA)*DMYENERGY**2)
        WRITE(LUNGFO,*)

        ICAL =1

      ENDIF

      IF (KHALBA.LT.0.AND.ABS(XIN-XCENHALBA).GT.ZLENHAL2) THEN
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        AXOUT=0.0D0
        AYOUT=0.0D0
        AZOUT=0.0D0
        RETURN
      ENDIF

      XKX=XKHALBA*(-ZIN)
      YKY=YKHALBA*YIN
      ZKZ=ZKHALBA*(XIN-XCENHALBA)

      DSNXKX=DSIN(XKX)
      DCSXKX=DCOS(XKX)
      DSHYKY=DSINH(YKY)
      DCHYKY=DSQRT(1.0D0+DSHYKY*DSHYKY)
      DSNZKZ=DSIN(ZKZ)
      DCSZKZ=DCOS(ZKZ)

      BXH=-XKHALBA/YKHALBA*B0HALBA*DSNXKX*DSHYKY*DCSZKZ
      BYH=                 B0HALBA*DCSXKX*DCHYKY*DCSZKZ
      BZH=-ZKHALBA/YKHALBA*B0HALBA*DCSXKX*DSHYKY*DSNZKZ

      AXH=B0HALBA/ZKHALBA*                DCSXKX*DCHYKY*DSNZKZ
      AYH=B0HALBA/ZKHALBA*XKHALBA/YKHALBA*DSNXKX*DSHYKY*DSNZKZ
      AZH=0.

      BZOUT=-BXH
      BYOUT=BYH
      BXOUT=BZH

      AZOUT=-AXH
      AYOUT= AYH
      AXOUT= AZH

      RETURN
      END
