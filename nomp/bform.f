*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.52/11 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.33  by  Michael Scheer
*CMZ :  2.13/09 08/03/2000  17.52.23  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.13.04  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE BFORM

*KEEP,gplhint.
*KEND.

C--- BERECHNET HARD-EDGE B-FELDKONFIGURATION


      IMPLICIT NONE

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,bfeld.
      include 'bfeld.cmn'
*KEND.

      DOUBLE PRECISION B0LM

      IF (IBSYM.EQ.0) STOP '*** S/R BFORM: IBSYM.EQ.0 ***'

CC    FB0M=FB0MFUN(FB0N)

      BBY1=B0FORM
      BBY2=-B0FORM/FB0M
      BBY3=B0FORM/FB0M/FB0N
      BBY4=-B0FORM/FB0N
      BBY5=BBY3

      BBY6=0.0
        BBY7=0.0

      B0LM=FB0N/2.*B0LP

      XM1=0.0
      XP1=B0LP/8.
      XM2=XP1
      XP2=XM2+B0LP/8.
      XM3=XP2
      XP3=XM3+B0LM/8.
      XM4=XP3
      XP4=XM4+B0LM/4.
      XM5=XP4
      XP5=XM5+B0LM/8.

C250991  XM6=XP5+1.
C250991  XP6=XM6+1.
      XM6=XP5  !C250991
      XP6=XM6  !C250991
        XM7=XP7
        XP7=XM7

      RETURN
      END
