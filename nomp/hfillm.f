*CMZ :  4.00/14 22/12/2021  13.41.56  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  12.18.27  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/05 12/06/2014  11.20.34  by  Michael Scheer
*CMZ :  2.70/07 14/01/2013  16.55.40  by  Michael Scheer
*CMZ :  2.68/05 28/09/2012  09.17.41  by  Michael Scheer
*CMZ :  2.67/01 21/02/2012  10.32.09  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  12.51.13  by  Michael Scheer
*-- Author :    Michael Scheer   21/01/2012
      subroutine hfillm(id,x,y,w)
*KEEP,gplhint.
*KEND.

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      integer id,idold
      real x,y
      double precision w

      data idold/-9999/


      if (iroottrees.ne.0) then
      endif

      idold=id


      call mh_hfill(id,dble(x),dble(y),w)

      return
      end
