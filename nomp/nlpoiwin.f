*CMZ :  3.05/06 17/07/2018  11.15.16  by  Michael Scheer
*CMZ :  2.67/04 25/10/2012  15.10.37  by  Michael Scheer
*CMZ :  2.16/08 12/08/2009  08.49.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.00/00 16/12/98  14.36.03  by  Michael Scheer
*-- Author :    Michael Scheer   15/12/98
      SUBROUTINE NLPOIWIN
*KEEP,gplhint.
*KEND.

*KEEP,trackf90u.
      include 'trackf90u.cmn'
*KEND.


C--- NLPOIWIN ESTIMATES VALUES FOR NLPOI AND WGWINFC
C--- FORMULA FROM X-RAY DATA BOOKLET

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,colli.
      include 'colli.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,track.
      include 'track.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER NLPOIMX,NLPERMN,NLPOIMN
      DOUBLE PRECISION B0DIP,YLOW,YHIGH,WGLOW,WGHIGH,ECDIP,PHIDEFL,DEFLEC,DLENG,DNPER

      DATA NLPOIMX/250000/,NLPERMN/50/

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     SR NLPOIWIN CALLED:'
      WRITE(LUNGFO,*)

      PHIDEFL=(PHIMX-PHIMN)/2.
      B0DIP=DSQRT(DABS(BMAXGL2))
      DEFLEC=PHIDEFL*DMYGAMMA
      DLENG=DEFLEC/93.4/B0DIP
      IF (DLENG.EQ.0.0) DLENG=(XMX-XMN)
      DNPER=(XMX-XMN)/DLENG
      IF (DNPER.LT.1.D0) DNPER=1.D0
      NLPOIMN=NLPERMN*DNPER*DEFLEC

      IF (B0DIP.EQ.0.D0) THEN
          WRITE(LUNGFO,*)'*** ERROR IN NLPOIWIN: B0DIP.EQ.0  ***'
          WRITE(LUNGFO,*)'CHECK MAG. FIELD'
          WRITE(LUNGFO,*)'*** PROGRAM WAVE ABORTED ***'
          WRITE(6,*)'*** ERROR IN NLPOIWIN: B0DIP.EQ.0  ***'
          WRITE(6,*)'CHECK MAG. FIELD'
          WRITE(6,*)'*** PROGRAM WAVE ABORTED ***'
          STOP
      ENDIF

      ECDIP=ecdipev1*DMYENERGY**2*B0DIP

      WRITE(LUNGFO,*)
     &'     mag. field [T], Ec [eV]:',SNGL(B0DIP),SNGL(ECDIP)
      WRITE(LUNGFO,*)

      YLOW=FREQ(1)/ECDIP
      YHIGH=FREQ(NFREQ)/ECDIP

      IF (YLOW.LT.1.) THEN
          WGLOW=0.408/DMYENERGY*YLOW**(-0.354)/1000.
      ELSE
          WGLOW=0.408/DMYENERGY*YLOW**(-0.549)/1000.
      ENDIF

      IF (YHIGH.LT.1.) THEN
          WGHIGH=0.408/DMYENERGY*YHIGH**(-0.354)/1000.
      ELSE
          WGHIGH=0.408/DMYENERGY*YHIGH**(-0.549)/1000.
      ENDIF

      WRITE(LUNGFO,*)'      WGWINFC low/high: '
     & ,SNGL(WGLOW*DMYGAMMA),SNGL(WGHIGH*DMYGAMMA)
      WRITE(LUNGFO,*)

      IF (WGWINFC.EQ.9999.) THEN
          WGWINFC=WGLOW*20.*DMYGAMMA
          WRITE(LUNGFO,*)'      WGWINFC set to: ',WGWINFC
      ENDIF

      IF (WGLOW.GT.PHIDEFL) THEN
          WGLOW=PHIDEFL
          WGWINFC=PHIDEFL*DMYGAMMA*1.1 !1.1 TO HAVE SOME OVERLAP
          WRITE(LUNGFO,*)
     &'      WGWINFC corresponds to angle greater than max. deflection'
          WRITE(LUNGFO,*)'      WGWINFC set to: ',WGWINFC
      ENDIF

      IF (WGHIGH.GT.PHIDEFL) THEN
          WGHIGH=PHIDEFL
      ENDIF

      IF (NLPOI.EQ.-9999) THEN
          NLPOI=1000*WGWINFC*WGLOW/WGHIGH
          IF (NLPOI.LT.NLPOIMN) THEN
         NLPOI=NLPOIMN
          ENDIF
          WRITE(LUNGFO,*)'      NLPOI set to: ',NLPOI
          IF (NLPOI.GT.NLPOIMX) THEN
          NLPOI=NLPOIMX
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'*** WARNING IN NLPOIWIN: NLPOI VERY LARGE'
          WRITE(LUNGFO,*)'      NLPOI limited to: ',NLPOI
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)'*** WARNING IN NLPOIWIN: NLPOI VERY LARGE'
          WRITE(6,*)'      NLPOI limited to: ',NLPOI
          WRITE(6,*)
          ENDIF
      ENDIF


      WRITE(LUNGFO,*)


      RETURN
      END
