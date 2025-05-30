*CMZ :  2.13/05 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.58.52  by  Michael Scheer
*CMZ : 00.01/09 06/10/95  16.26.59  by  Michael Scheer
*CMZ :  1.00/01 26/09/95  16.53.15  by  Michael Scheer
*CMZ :  1.00/00 21/09/95  17.14.11  by  Michael Scheer
*-- Author :    Michael Scheer   21/09/95
      SUBROUTINE VPOT3DFIT(NPOI,X,Y,Z,BX,BY,BZ,NDIMP,N,MORD
     &                       ,INDEX,VC,NSTAKP,FSTAK,A,B,WS,WA,WB
     &                       ,XPOW,YPOW,ZPOW,IFAIL)
*KEEP,gplhint.
*KEND.

C--- FITS V=...+VC(IXYZ)*X**(INDEX(1,IXYZ))*Y**(INDEX(2,IXYZ))*Z**(INDEX(3,IXYZ)
C--- WITH (BX,BY,BZ)=-GRAD(V)

C---  INPUT:
C-    NPOI      : NUMBER OF POINTS X,F(X)
C-    NDIMP       :  DIMENSION PARAMETER OF ARRAYS A,B,WS,WA,WB,
C-    N           :  NUMBER OF COEFFICIENTS VC
C-    X,Y,Z     : ARGUMENTS OF FUNCTION F(X,Y,Z)
C-    BX,BY,BZ    :  MAGNETIC FIELD COMPONENTS
C-    MORD      :   HIGHEST ORDER OF FITTED COEFFS
C-    INDEX     :   INDEX(3,N) IS POINTER TO INDICES I,J,K
C-    NSTAKP       : STACK DIMENSION
C-    FSTAK     :   STACK OF COEFFICIENTS FSTAK(4,NSTAKP,N)

C---  OUTPUT:
C-    VC     : COEFFICIENTS VC(N) TO BE FITTED
C-    IFAIL       :  FAILURE FLAG

C---  WORKING SPACE:
C-    A  :  MATRIX OF EQUATION SYSTEM A(NDIMP,NDIMP)
C-    B  :  INHOMOGENITY OF EQUATION SYSTEM B(NDIMP)
C-    WS :  WORKING SPACE WS(2*NDIMP)
C-    WA :  WORKING SPACE WA(NDIMP,NDIMP)
C-    WB :  WORKING SPACE WB(NDIMP)
C-    XPOW  :  WORKING SPACE XPOW(N)
C-    YPOW  :  WORKING SPACE YPOW(N)
C-    ZPOW  :  WORKING SPACE ZPOW(N)

      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      INTEGER NPOI,N,MORD,INDEX(4,N),IFAIL,NDIMP,NSTAKP
      INTEGER IX,IY,IZ,JX,JY,JZ,IXYZ,JXYZ,IPOI,IFIT,JFIT,NFIT,IS,NS,JS,MS
      INTEGER IDELTA

      DOUBLE PRECISION X(NPOI),Y(NPOI),Z(NPOI),BX(NPOI),BY(NPOI),BZ(NPOI)
     &                  ,VC(N),DELTA,DELTAMX
     &                  ,A(NDIMP,NDIMP),WA(NDIMP,NDIMP)
     &                  ,B(NDIMP),WS(2*NDIMP),WB(NDIMP)
     &                  ,XPOW(MORD+1),YPOW(MORD+1),ZPOW(MORD+1)
     &                  ,FSTAK(4,NSTAKP,N)

      DOUBLE PRECISION AX1,AY1,AZ1,AX2,AY2,AZ2,XPOI,YPOI,ZPOI
     &                  ,XPOW1,YPOW1,ZPOW1
     &                  ,XPOW2,YPOW2,ZPOW2
     &         ,BSUM,X99,Y99,Z99
      IFAIL=0

      IF (N.GT.NPOI) THEN
          WRITE(LUNGFO,*)
     &      '*** WARNING SR VPOT3DFIT: NOT ENOUGH DATA'
          WRITE(6,*)
     &      '*** WARNING SR VPOT3DFIT: NOT ENOUGH DATA'
          IFAIL=999
          RETURN
      ENDIF

      DO IXYZ=1,N
          B(IXYZ)=0.D0
          VC(IXYZ)=0.D0
      ENDDO

      DO IY=1,N
      DO IX=1,N
          A(IX,IY)=0.D0
      ENDDO
      ENDDO

C- DO THE FITTING

      NFIT=N

      DO IPOI=1,NPOI

          XPOI=X(IPOI)
          YPOI=Y(IPOI)
          ZPOI=Z(IPOI)

          XPOW(1)=1.D0
          DO IX=2,MORD+1
          XPOW(IX)=XPOW(IX-1)*XPOI
        ENDDO

          YPOW(1)=1.D0
          DO IY=2,MORD+1
          YPOW(IY)=YPOW(IY-1)*YPOI
        ENDDO

          ZPOW(1)=1.D0
          DO IZ=2,MORD+1
          ZPOW(IZ)=ZPOW(IZ-1)*ZPOI
        ENDDO

        IF (BX(IPOI).EQ.-9999.) THEN
          X99=0.0D0
        ELSE
          X99=1.0D0
        ENDIF

        IF (BY(IPOI).EQ.-9999.) THEN
          Y99=0.0D0
        ELSE
          Y99=1.0D0
        ENDIF

        IF (BZ(IPOI).EQ.-9999.) THEN
          Z99=0.0D0
        ELSE
          Z99=1.0D0
        ENDIF


      IFIT=0
      DO IXYZ=1,NFIT

         IFIT=IFIT+1

           NS=INDEX(4,IXYZ)

           BSUM=0.D0
           AX1=0.D0
           AY1=0.D0
           AZ1=0.D0
           DO IS=1,NS

            IX=NINT(FSTAK(1,IS,IXYZ))+1
            IY=NINT(FSTAK(2,IS,IXYZ))+1
            IZ=NINT(FSTAK(3,IS,IXYZ))+1

                      IF(IX.GT.1) THEN
                         XPOW1=XPOW(IX-1)
            ELSE
                         XPOW1=1.D0
                      ENDIF

                      IF(IY.GT.1) THEN
                         YPOW1=YPOW(IY-1)
            ELSE
                         YPOW1=1.D0
                      ENDIF

                      IF(IZ.GT.1) THEN
                         ZPOW1=ZPOW(IZ-1)
            ELSE
                         ZPOW1=1.D0
                      ENDIF

            BSUM=BSUM+FSTAK(4,IS,IXYZ)*(
     &                  -(IX-1)*BX(IPOI)*XPOW1*YPOW(IY)*ZPOW(IZ)*X99
     &         -(IY-1)*BY(IPOI)*YPOW1*XPOW(IX)*ZPOW(IZ)*Y99
     &              -(IZ-1)*BZ(IPOI)*ZPOW1*XPOW(IX)*YPOW(IY)*Z99
     &                  )

                  AX1=AX1
     &                     +X99*FSTAK(4,IS,IXYZ)*(IX-1)*XPOW1*YPOW(IY)*ZPOW(IZ)
                      AY1=AY1
     &                     +Y99*FSTAK(4,IS,IXYZ)*(IY-1)*YPOW1*XPOW(IX)*ZPOW(IZ)
            AZ1=AZ1
     &                      +Z99*FSTAK(4,IS,IXYZ)*(IZ-1)*ZPOW1*XPOW(IX)*YPOW(IY)

         ENDDO !NS

                      B(IFIT)=B(IFIT)+BSUM

      JFIT=0
      DO JXYZ=1,NFIT

         JFIT=JFIT+1

           MS=INDEX(4,JXYZ)

           AX2=0.D0
           AY2=0.D0
           AZ2=0.D0
           DO JS=1,MS

            JX=NINT(FSTAK(1,JS,JXYZ))+1
            JY=NINT(FSTAK(2,JS,JXYZ))+1
            JZ=NINT(FSTAK(3,JS,JXYZ))+1

                      IF(JX.GT.1) THEN
                         XPOW2=XPOW(JX-1)
            ELSE
                         XPOW2=1.D0
                      ENDIF

                      IF(JY.GT.1) THEN
                         YPOW2=YPOW(JY-1)
            ELSE
                         YPOW2=1.D0
                      ENDIF

                      IF(JZ.GT.1) THEN
                         ZPOW2=ZPOW(JZ-1)
            ELSE
                         ZPOW2=1.D0
                      ENDIF

                  AX2=AX2
     &                     +X99*FSTAK(4,JS,JXYZ)*(JX-1)*XPOW2*YPOW(JY)*ZPOW(JZ)
                      AY2=AY2
     &                     +Y99*FSTAK(4,JS,JXYZ)*(JY-1)*YPOW2*XPOW(JX)*ZPOW(JZ)
            AZ2=AZ2
     &                      +Z99*FSTAK(4,JS,JXYZ)*(JZ-1)*ZPOW2*XPOW(JX)*YPOW(JY)

         ENDDO !NS

         A(IFIT,JFIT)=A(IFIT,JFIT)+AX1*AX2+AY1*AY2+AZ1*AZ2

      ENDDO !JXYZ
      ENDDO !IXYZ

      ENDDO !IPOI

      DO IXYZ=1,NFIT
         WB(IXYZ)=B(IXYZ)
      ENDDO

      DO IY=1,NFIT
      DO IX=1,NFIT
              WA(IX,IY)=A(IX,IY)
      ENDDO
      ENDDO

      DO IXYZ=1,NFIT
         WB(IXYZ)=B(IXYZ)
      ENDDO

      DO IY=1,NFIT
      DO IX=1,NFIT
              WA(IX,IY)=A(IX,IY)
      ENDDO
      ENDDO

        CALL DEQN(NFIT,WA,NDIMP,WS,IFAIL,1,WB) !CERN F010

      DO IFIT=1,NFIT
         VC(IFIT)=WB(IFIT)
      ENDDO

      DELTAMX=-1.D30
      DO IX=1,NFIT
      DELTA=0.D0
      DO IY=1,NFIT
          DELTA=DELTA+A(IX,IY)*WB(IY)
      ENDDO
      DELTA=DABS(DELTA-B(IX))
      IF (DELTA.GT.DELTAMX) THEN
          DELTAMX=DELTA
          IDELTA=IX
      ENDIF
      ENDDO

      IF (DELTAMX.GT.1.D-5) THEN
          WRITE(LUNGFO,*)
     &'*** WARNING SR VPOT3DFIT: PROBABLY PROBLEM WITH LINEAR EQUATION SYSTEM'
          WRITE(6,*)
     &'*** WARNING SR VPOT3DFIT: PROBABLY PROBLEM WITH LINEAR EQUATION SYSTEM'
      WRITE(6,*)'     max. error of DEQN and corresponding coefficient:'
      WRITE(6,*)'     ',SNGL(DELTAMX),SNGL(VC(IDELTA))
     &                   ,' = C',INDEX(1,IDELTA)
     &                   ,INDEX(2,IDELTA),INDEX(3,IDELTA)
      ENDIF

      WRITE(LUNGFO,*)
      WRITE(LUNGFO,*)'     max. error of DEQN and corresponding coefficient:'
      WRITE(LUNGFO,*)'     ',SNGL(DELTAMX),SNGL(VC(IDELTA))
     &                       ,' = C',INDEX(1,IDELTA)
     &                        ,INDEX(2,IDELTA),INDEX(3,IDELTA)
      WRITE(LUNGFO,*)

      RETURN
      END
