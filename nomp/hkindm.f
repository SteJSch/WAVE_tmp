*CMZ :  4.00/14 31/12/2021  10.47.16  by  Michael Scheer
*CMZ :  3.02/03 10/11/2014  11.44.31  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/07 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/06 17/06/2014  16.39.04  by  Michael Scheer
*CMZ :  2.68/05 24/09/2012  11.16.28  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  13.23.56  by  Michael Scheer
*-- Author :    Michael Scheer   17/02/2012
      subroutine hkindm(id,kind4,chopt)
*KEEP,gplhint.
*KEND.

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      integer*8 kind8(32)
      integer id,kind4(32)
      character(*) chopt

      kind4=0

      if (chopt.eq.'A'.and.iroottrees.ge.0) then
        call mh_kind(id,kind4)
      endif

      return
      end
