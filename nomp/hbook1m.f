*CMZ :  4.00/14 30/12/2021  15.41.22  by  Michael Scheer
*CMZ :  4.00/13 16/12/2021  12.18.27  by  Michael Scheer
*CMZ :  3.02/04 12/12/2014  15.28.39  by  Michael Scheer
*CMZ :  3.02/00 01/09/2014  11.11.30  by  Michael Scheer
*CMZ :  3.01/06 23/06/2014  16.20.53  by  Michael Scheer
*CMZ :  3.01/05 11/06/2014  16.12.58  by  Michael Scheer
*CMZ :  2.69/00 24/10/2012  16.50.54  by  Michael Scheer
*CMZ :  2.67/00 17/02/2012  16.23.54  by  Michael Scheer
*-- Author :    Michael Scheer   19/01/2012
      subroutine hbook1m(id,ctitle,nbins,xmin,xmax,vmx)
*KEEP,gplhint.
*KEND.

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

c creates 1D-histo for root.

      integer id,nbins,k,i,null
      real xmin,xmax,vmx
      character(*) ctitle
      character c1
      character(8) cname,cname2
      equivalence (c1,null)

      data null/0/






      call mh_book1(id,ctitle,nbins,dble(xmin),dble(xmax))
      return
      end

