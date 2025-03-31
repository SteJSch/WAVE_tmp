*CMZ :  2.68/04 03/09/2012  13.22.46  by  Michael Scheer
*CMZ :  2.68/03 22/08/2012  14.13.20  by  Michael Scheer
*CMZ :  2.66/21 22/11/2011  13.51.05  by  Michael Scheer
*CMZ :  2.66/20 22/11/2011  10.17.35  by  Michael Scheer
*-- Author :    Michael Scheer   18/11/2011
      subroutine estep(x,y,z,vx,vy,vz,dt,gamma,dgamma)
*KEEP,gplhint.
*KEND.

c Simple approach to take electrical fields into account
c The energy and gamma of the particle is not changed, but returned in dgamma

*KEEP,primkin.
      include 'primkin.cmn'
*KEEP,efield.
      include 'efield.cmn'
*KEEP,b0scglob.
      include 'b0scglob.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      double precision x,y,z,vx,vy,vz,dt,gamma,dgamma
      double precision px,py,pz,qdt,eg,vn,pn,
     &  px0,py0,pz0,pn0,pn2,pn02,de,dpx,dpy,dpz,dtm

      if (efieldx.eq.0.0d0.and.efieldy.eq.0.0d0.and.efieldz.eq.0.0d0) then
        dgamma=0.0d0
        return
      endif

      qdt=icharge*echarge1*dt
      eg=emasskg1*gamma

      px0=vx*eg !SI-units
      py0=vy*eg
      pz0=vz*eg
      pn02=px0*px0+py0*py0+pz0*pz0
      pn0=sqrt(pn02)

      dpx=efieldx*qdt !SI-units
      dpy=efieldy*qdt
      dpz=efieldz*qdt

      px=px0+dpx !SI-units
      py=py0+dpy
      pz=pz0+dpz

      vn=sqrt(vx*vx+vy*vy+vz*vz)
      pn2=px*px+py*py+pz*pz
      pn=sqrt(pn2)

c total momentum and energy are kept!!
      vx=px/pn*vn
      vy=py/pn*vn
      vz=pz/pn*vn

      dtm=0.5d0*dt/(emasskg1*gamma) ! F=dp/dt, m=m_e*gamma, a=F/m*dp/dt/m_e/gamma

      x=x+dtm*dpx
      y=y+dtm*dpy
      z=z+dtm*dpz

      de=(vx*dpx+dpy*vy+dpz*vz)/echarge1/1.0d9 !GeV
      dgamma=de/dmyenergyP*dmygammaP

      return
      end
