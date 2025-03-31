*CMZ :  4.00/14 30/12/2021  10.08.13  by  Michael Scheer
*CMZ :  3.02/05 22/03/2015  18.54.38  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  15.54.08  by  Michael Scheer
*-- Author :    Michael Scheer   21/01/2012
      subroutine hoperam(id1,coperator,id2,id3,w1,w2)
*KEEP,gplhint.
*KEND.

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      integer id1,id2,id3
      real w1,w2
      character(*) coperator

      call mh_opera(id1,coperator,id2,id3,dble(w1),dble(w2))

      return
      end
