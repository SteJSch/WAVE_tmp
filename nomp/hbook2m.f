*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 19/12/2021  10.13.31  by  Michael Scheer
*CMZ :  3.02/04 12/12/2014  15.28.39  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  2.69/00 24/10/2012  16.55.08  by  Michael Scheer
*CMZ :  2.67/04 14/05/2012  13.05.02  by  Michael Scheer
*CMZ :  2.67/01 21/02/2012  10.22.54  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  16.23.54  by  Michael Scheer
*-- Author :    Michael Scheer   19/01/2012
      subroutine hbook2m(id,ctitle,nxbins,xmin,xmax,nybins,ymin,ymax,vmx)
*KEEP,gplhint.
*KEND.

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

c creates 2D-histo for root.

      integer id,nxbins,nybins,k,i
      real xmin,xmax,ymin,ymax,vmx
      character(*) ctitle
      character c1
      character(8) cname,cname2




        call mh_book2(id,ctitle,nxbins,dble(xmin),dble(xmax),
     &    nybins,dble(ymin),dble(ymax))



      return
      end
