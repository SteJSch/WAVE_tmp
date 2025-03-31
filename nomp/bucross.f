*CMZ :  3.00/01 26/03/2013  19.45.55  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.70/06 04/01/2013  13.24.05  by  Michael Scheer
*CMZ :  2.67/02 03/05/2012  09.35.17  by  Michael Scheer
*CMZ :  2.20/09 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.34  by  Michael Scheer
*CMZ :  2.13/05 08/02/2000  17.24.35  by  Michael Scheer
*CMZ :  1.03/06 10/06/98  14.43.01  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  15.31.37  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.48.51  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.36  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BUCROSS(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,
     &                               AXOUT,AYOUT,AZOUT)
*KEEP,gplhint.
*KEND.


C--- MAGNETIC FIELD OF A CROSSED UNDULATOR. FIRST UNDULATOR IS HORIZONTAL,
C    MODULATOR IS HORIZONTAL, SECOND UNDULATOR IS VERTICAL
C    MODULATOR CENTER IS DEVICE CENTER

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,halbasy.
      include 'halbasy.cmn'
*KEEP,ucross.
      include 'ucross.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      INTEGER ICAL,I

      DOUBLE PRECISION X,Y,Z,BY,BZ,BXOUT,BYOUT,BZOUT,AY,AZ,AXOUT,AYOUT,AZOUT
      DOUBLE PRECISION XIN,YIN,ZIN,xcen

      DOUBLE PRECISION WLEN1,EHARM1(3),PARK(3)
      double precision b0halbasyo,xkhalbasyo,zlhalbasyo,zkhalbasyo,fasymo,
     &  ahwpolo

      DATA ICAL/0/

      b0halbasyo=b0halbasy
      xkhalbasyo=xkhalbasy
      zlhalbasyo=zlhalbasy
      zkhalbasyo=zkhalbasy
      fasymo=fasym
      ahwpolo=ahwpol
      xcen=xcenhal
      xcenhal=0.0d0

      IF (IUCRSAX.NE.0) THEN
          CALL BUCROSS_AX(XIN,YIN,ZIN,BXOUT,BYOUT,BZOUT,
     &      AXOUT,AYOUT,AZOUT)
          goto 9999
      ENDIF

C--- FIRST CALL

      IF (ICAL.EQ.0) THEN

        DO I=1,3
          ZKUHAL(I)=2.*PI1/ZLUHAL(I)
          PARK(I)=
     &      ECHARGE1*DABS(B0UCROSS(I))*ZLUHAL(I)/(2.*PI1*EMASSKG1*CLIGHT1)
          WLEN1=(1+PARK(I)**2/2.)/2./DMYGAMMA**2*ZLUHAL(I)*1.D9
          IF (WLEN1.NE.0.0) EHARM1(I)=WTOE1/WLEN1
          ucharm1=eharm1
        ENDDO   !I


      IF (UASYM(2).NE.2.D0) THEN
          ULIMI(2)=-ZLUHAL(2)/2.*(DFLOAT(NMUPOL(2))+UASYM(2))/2.
      ELSE
          ULIMI(2)=-ZLUHAL(2)*(DFLOAT(NMUPOL(2)-1)/2.D0+1.D0)/2.
      ENDIF
          ULIMI(3)=-ULIMI(2)

      IF (UASYM(1).NE.2.D0) THEN
          ULIMI(1)=ULIMI(2)-ZLUHAL(1)/2.*(DFLOAT(NMUPOL(1))+UASYM(1))
      ELSE
          ULIMI(1)=ULIMI(2)-ZLUHAL(1)*(DFLOAT(NMUPOL(1)-1)/2.D0+1.D0)
      ENDIF

      IF (UASYM(3).NE.2.D0) THEN
          ULIMI(4)=ULIMI(3)+ZLUHAL(3)/2.*(DFLOAT(NMUPOL(3))+UASYM(3))
      ELSE
          ULIMI(4)=ULIMI(3)+ZLUHAL(3)*(DFLOAT(NMUPOL(3)-1)/2.D0+1.D0)
      ENDIF


          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     SR BUCROSS: Parameter of crossed undulator:'
          WRITE(LUNGFO,*)
          WRITE(LUNGFO,*)'     Peak field B0:    ',(SNGL(B0UCROSS(I)),I=1,3)
          WRITE(LUNGFO,*)'     Number of main poles:',(NMUPOL(I),I=1,3)
          WRITE(LUNGFO,*)'     Period lengths:   ',(SNGL(ZLUHAL(I)),I=1,3)
          WRITE(LUNGFO,*)'     Limits of modules:'
          WRITE(LUNGFO,*)'     ',(SNGL(ULIMI(I)),I=1,4)
          WRITE(LUNGFO,*)'     Deflection parameter K:',(SNGL(PARK(I)),I=1,3)
          WRITE(LUNGFO,*)'     First harmonical [eV]: ',(SNGL(EHARM1(I)),I=1,3)
          WRITE(LUNGFO,*)

          ICAL=1

      ENDIF

      IF (XIN.LT.ULIMI(1).OR.XIN.GT.ULIMI(4))THEN

          BXOUT=0.
          BYOUT=0.
          BZOUT=0.
          AXOUT=0.
          AYOUT=0.
          AZOUT=0.

      ELSE IF(XIN.LT.ULIMI(2)) THEN

C HORIZONTAL

      B0HALBASY=B0UCROSS(1)
      XKHALBASY=0.
      ZLHALBASY=ZLUHAL(1)
      ZKHALBASY=ZKUHAL(1)
      YKHALBASY=ZKHALBASY
      FASYM=UASYM(1)
      AHWPOL=DFLOAT(NMUPOL(1))

      X=XIN-ULIMI(1)-(ULIMI(2)-ULIMI(1))/2.
      CALL BHALBASY(X,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

      ELSE IF(XIN.LT.ULIMI(3)) THEN

C MODULATOR

      B0HALBASY=B0UCROSS(2)
      XKHALBASY=0.
      ZLHALBASY=ZLUHAL(2)
      ZKHALBASY=ZKUHAL(2)
      YKHALBASY=ZKHALBASY
      FASYM=UASYM(2)
      AHWPOL=DFLOAT(NMUPOL(2))

      X=XIN
      CALL BHALBASY(X,YIN,ZIN,BXOUT,BYOUT,BZOUT,AXOUT,AYOUT,AZOUT)

      ELSE IF(XIN.LT.ULIMI(4)) THEN

C VERTICAL


      B0HALBASY=B0UCROSS(3)
      XKHALBASY=0.
      ZLHALBASY=ZLUHAL(3)
      ZKHALBASY=ZKUHAL(3)
      YKHALBASY=ZKHALBASY
      FASYM=UASYM(3)
      AHWPOL=DFLOAT(NMUPOL(3))

      X=XIN-ULIMI(3)-(ULIMI(4)-ULIMI(3))/2.
      Z=YIN
      Y=ZIN

      CALL BHALBASY(X,Y,Z,BXOUT,BY,BZ,AXOUT,AY,AZ)

      BZOUT= BY
      BYOUT= BZ
      AZOUT= AY
      AYOUT= AZ

      ENDIF

9999  xcenhal=xcen
      b0halbasy=b0halbasyo
      xkhalbasy=xkhalbasyo
      zlhalbasy=zlhalbasyo
      zkhalbasy=zkhalbasyo
      fasym=fasymo
      ahwpol=ahwpolo

      RETURN
      END
