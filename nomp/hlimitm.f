*CMZ :  4.00/14 21/12/2021  10.43.42  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  12.18.27  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 20/06/2014  17.35.22  by  Michael Scheer
*CMZ :  3.01/05 11/06/2014  16.01.34  by  Michael Scheer
*-- Author :    Michael Scheer   11/06/2014
      subroutine hlimitm(limit)
*KEEP,gplhint.
*KEND.
      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      integer limit

      call mh_limit(limit)

      return
      end
