*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.20/09 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ : 00.01/07 24/02/95  09.42.55  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.47.19  by  Michael Scheer
*CMZ : 00.00/07 01/06/94  10.41.53  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BELLCRS(XIN,Y,Z,BX,BY,BZ,AX,AY,AZ)

*KEEP,gplhint.
*KEND.

C--- TWO TIMES ELLIPTICAL UNDULATOR WITH MODULATOR INBETWEEN

      DOUBLE PRECISION X,Y,Z,BX,BY,BZ,AX,AY,AZ
      DOUBLE PRECISION XIN,XLENEU,XLENHU,TOTLEN,XLENEU2,XLENHU2,TOTLEN2
      DOUBLE PRECISION B0H,B0V

      INTEGER ICAL

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,ellip.
      include 'ellip.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      DATA ICAL/0/

      IF (ICAL.EQ.0) THEN
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)
     &'     SR BELLCRS: VARIABLES OF NAMELIST ELLIP AND HALBASY TAKEN FOR FIRST'
          WRITE(LUNGFO,*)
     &'                 UNDULATOR AND MODULATOR, NAMELIST USERN (USER(1:2) USED'
          WRITE(LUNGFO,*)
     &'                 TO SCALE SECOND UNDULATOR'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)

          B0H=B0ELLIPH
          B0V=B0ELLIPV
          ICAL=1
      ENDIF !ICAL

      XLENEU=DABS((PERELLIP+2.5+ELLSHFT)*XLELLIP)
      IF (FASYM.NE.2.D0) THEN
              XLENHU=ZLHALBASY*(AHWPOL+FASYM)/2.D0
      ELSE
              XLENHU=ZLHALBASY*((AHWPOL-1.D0)/2.+1.D0)
      ENDIF
      TOTLEN=2.D0*XLENEU+XLENHU
      TOTLEN2=TOTLEN/2.D0
      XLENEU2=XLENEU/2.D0
      XLENHU2=XLENHU/2.D0

      IF (XIN.LT.-TOTLEN2) THEN
          BX=0.0
          BY=0.0
          BZ=0.0
          AX=0.0
          AY=0.0
          AZ=0.0
      ELSEIF(XIN.LT.-XLENHU2) THEN  !FIRST ELLIPTICAL UNDULATOR
          X=XIN+XLENHU2+XLENEU2
            B0ELLIPH=B0H
            B0ELLIPV=B0V
          CALL BELLIP(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
      ELSEIF(XIN.LT.XLENHU2) THEN   !MODULATOR
          CALL BHALBASY(XIN,Y,Z,BX,BY,BZ,AX,AY,AZ)
      ELSEIF(XIN.LT.TOTLEN2) THEN   !SECOND ELLIPTICAL UNDULATOR
          X=XIN-XLENHU2-XLENEU2
            B0ELLIPH=B0H*USER(1)
            B0ELLIPV=B0V*USER(2)
          CALL BELLIP(X,Y,Z,BX,BY,BZ,AX,AY,AZ)
      ELSEIF(XIN.GT.TOTLEN2) THEN
          BX=0.0
          BY=0.0
          BZ=0.0
          AX=0.0
          AY=0.0
          AZ=0.0
      ENDIF


      RETURN
      END
