*CMZ :          24/08/2024  12.47.01  by  Michael Scheer
*CMZ :  4.01/05 26/04/2024  07.41.13  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  13.39.24  by  Michael Scheer
*CMZ :  4.01/02 14/05/2023  11.47.49  by  Michael Scheer
*CMZ :  4.01/00 22/02/2023  14.34.04  by  Michael Scheer
*CMZ :  4.00/17 05/12/2022  10.30.41  by  Michael Scheer
*CMZ :  4.00/16 17/09/2022  15.46.32  by  Michael Scheer
*CMZ :  4.00/15 02/06/2022  09.45.10  by  Michael Scheer
*CMZ :  4.00/11 28/06/2021  10.33.06  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_amprep(modewave)

      use omp_lib
      use uradphasemod

      implicit none

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEEP,track.
      include 'track.cmn'
*KEND.
cc+seq,uservar.

      complex*16 :: cde,czero=(0.0d0,0.0d0)
      double precision :: h2,ddist,wlen,dphi,phase0,cjvsto(4,3)

      double complex , dimension (:,:), allocatable :: aradbuff
      double complex , dimension (:,:,:), allocatable :: arad

      double precision, dimension (:), allocatable :: frq
      double precision, dimension (:,:), allocatable :: wsstokes,pow
      double precision, dimension (:,:,:,:,:), allocatable :: stokesprop
      double precision, dimension (:,:,:), allocatable :: fbunch,stokes

      complex*16, dimension (:), allocatable :: expphiran
      real, dimension (:), allocatable :: pherr,pherrc,phiran
      real, dimension(:,:), allocatable :: pranall,eall

      real eran(6),pran(3),rr(2)

      double complex :: apol,amp0(6),damp(6),amp(6),zexp,
     &  apolh,apolr,apoll,apol45,stokesv(4,3),cero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double precision :: t,udgamtot,upow,vf0,vn,vx0,vx2,vxf0,vxi,vy0,vy2,vyf0,
     &  vyi,vz0,vz2,vzf0,vzi,wlen1,x0,x2,xf0,xi,xlell,y0,y2,yf0,yi,ypi,yy,yyp,
     &  z0,z2,zf0,zi,zpi,zz,zzp,fillb(41),stok1,stok2,stok3,stok4,speknor,
     &  sqnbunch,sqnphsp,specnor,sbnor,rpin,r00(3),xph0,
     &  r(3),r0(3),pw,ph,phsum,pkerr,ppin,parke,pc(3),pcbrill(3),om1,
     &  park,pr,hbarev,obs(3),om,fhigh,flow,gamma,eix,eiy,eiz,emassg,
     &  efx,efy,efz,eharm1,ecdipev,ebeam,dtpho,dt,dtelec,dd0,debeam,
     &  drn0(3),drn00(3),ds,dr0(3),dr00(3),drn(3),dpp,dph,dist,dist0,dobs(3),
     &  bunnor,clight,bunchx,beta,beff,spow,
     &  zp0,yp0,rph,anor,fsum,smax,zob,yob,
     &  xkellip,zampell,yampell,parkv,parkh,zpampell,ypampell,emom,dzpin,dypin,zmin,ymin,phgsh

      double precision xprop,yprop(npinyprop_u),zprop(npinzprop_u),dy,dz,pinwprop,pinhprop
      double complex, dimension(:,:,:,:,:), allocatable :: fprop
      !double complex, dimension(:,:,:), allocatable :: fpriv
      double complex :: fpriv(3,npinzprop_u,npinyprop_u)

      double complex, dimension (:), allocatable ::
     &  uampex,uampey,uampez,uampbx,uampby,uampbz
      double precision, dimension (:,:), allocatable :: utraxyz,ustokes

      double complex :: rea(3),expsh

      integer :: kfreq,iobsv,i,np2,nelec,mbunch,meinbunch,ibu,jbun,
     &  kran=6,icbrill,ilo,kobsv,i1,i2,n,
     &  ifail,ndimu,nstepu,ith,noespread,noemit,jbunch,jubunch,jhbunch,
     &  jcharge=-1,lmodeph,nclo,jeneloss=0,iamppin,
     &  iamppincirc=0,ifrob,iobfr,isub,jvelofield=0,nlbu=0,nepho,ielo,
     &  modewave,iepho,ipobfr,ifieldprop,nzprop,nyprop,im,izm,iym,ifix

      integer, dimension (:), allocatable :: lnbunch

      integer :: idebug=0, lbunch=0, ierr=0, ielec=0
      integer ibunch,ihbunch,mthreads,nobsv,nobsvo,iemit,noranone,iz,iy,ipz,ipy,nobsvz,nobsvy
      integer iobm,iobp,iobfrm,iobfrp

c      integer iuser
c      iuser=user(3)

      nelec_u=max(1,nelec_u)
      mthreads_u=max(1,mthreads_u)
      nelec=nelec_u

      nobsvy=npiny_u
      nobsvz=npinz_u
      nobsvo=npinzo_u*npinyo_u
      dzpin=pinw_u/max(1,npinzo_u-1)
      dypin=pinh_u/max(1,npinyo_u-1)
      zmin=pincen_u(3)-pinw_u/2.0d0
      ymin=pincen_u(2)-pinh_u/2.0d0

      if (nelec_u.gt.1) then
        if (nelec_u.lt.mthreads_u.and.nelec_u.gt.1) then
          mthreads_u=nelec_u
        else
          nelec_u=max(mthreads_u,nelec_u/mthreads_u*mthreads_u)
        endif
      endif

      noranone=noranone_u

      if (nelec_u.eq.1.and.noranone.ne.0) then
        ibunch=0
      else
        ibunch=1
      endif

      ihbunch=ihbunch_u
      mthreads=mthreads_u

      if (modepin_u.eq.1) then
        iamppin=3
        nobsv=1
      else
        iamppin=1
        nobsv=npiny_u*npinz_u
      endif

      icbrill=nobsv/2+1

c      jhbunch=max(0,ihbunch)
      jhbunch=ihbunch
      meinbunch=nelec_u

      lbunch=0

c      if (jhbunch.ne.0) then
c        jhbunch=max(1,jhbunch)
c      endif

      nepho=nepho_u
      if (ifieldprop_u.eq.2) then
        allocate(
c     &    fpriv(3,npinzprop_u,npinyprop_u),
     &    fprop(3,npinzprop_u,npinyprop_u,nepho_u,mthreads),
     &    stokesprop(4,npinzprop_u,npinyprop_u,nepho_u,mthreads))
        stokesprop=0.0d0
      endif

      if (modepin_u.ne.0) then
        allocate(fieldbunch(7,npinzo_u,npinyo_u,nepho_u),stat=ierr)
        if (ierr.ne.0) then
          print*,""
          print*,"*** Warning in urad_amprep: Could not allocate buffer for beam Ntuple ***"
          print*,""
          return
        endif
        fieldbunch=czero
      endif

      call urad_field_ini(perlen_u,shift_u,beffv_u,beffh_u,modewave)

      if (perlen_u.ne.0.0d0) then
        emom=emasse1*dsqrt((gamma_u-1.0d0)*(gamma_u+1.0d0))
c*** OBSOLITE, SEE z0= further down
        xkellip=twopi1/perlen_u
        zampell=beffv_u*clight1/emom/xkellip**2
        yampell=beffh_u*clight1/emom/xkellip**2
c        zampell=zmx
c        yampell=ymx
        parkh=echarge1*dabs(beffh_u)*perlen_u/(twopi1*emasskg1*clight1)
        parkv=echarge1*dabs(beffv_u)*perlen_u/(twopi1*emasskg1*clight1)
        zpampell=parkv/gamma_u
c        print*,zpampell
c        ypampell=parkh/gamma_u
c        zpampell=tan(phimx)
c        print*,zpampell
c        stop
      else
        print*,''
        print*,'*** Error in urad_amprep: Zero period-length of undulator ***'
        print*,''
        stop
      endif

      dr00=[1.0d0,0.0d0,0.0d0]
      drn00=dr00/norm2(dr00)
      dr00=drn00*perlen_u
      r00=[0.0d0,0.0d0,0.0d0]

      x0=-perlen_u/2.0d0
      y0=0.0d0
      z0=0.0d0

      beta=dsqrt((1.0d0-1.0d0/gamma_u)*(1.0d0+1.0d0/gamma_u))

      clight=clight1
      hbarev=hbarev1
      ecdipev=ecdipev1
      emassg=emassg1

      if (modewave.eq.0) then
        z0=-zampell*cos(shift_u/2.0d0/perlen_u*twopi1)
        zp0=zpampell*sin(shift_u/2.0d0/perlen_u*twopi1)
        y0=-yampell*cos(shift_u/2.0d0/perlen_u*twopi1)
        yp0=-ypampell*sin(shift_u/2.0d0/perlen_u*twopi1)
      else
        y0=ytrack
        z0=ztrack
        zp0=vztrack/vxtrack
        yp0=vytrack/vxtrack
      endif

      xf0=-x0
      yf0=y0
      zf0=z0

      vn=clight*beta

      vx0=vn/sqrt(1.0d0+(zp0**2+yp0**2))
      vy0=vn*yp0
      vz0=vn*zp0

      vxf0=vx0
      vyf0=vy0
      vzf0=vz0

      vxi=vx0
      vyi=vy0
      vzi=vz0

      r0=r00
      dr0=dr00
      drn0=drn00
      r=r0
      drn=drn0

      vf0=norm2([vxf0,vyf0,vzf0])
      efx=vxf0/vf0
      efy=vyf0/vf0
      efz=vzf0/vf0

      nclo=nint(perlen_u/step_u)+1

      ds=step_u
c      dtim0=ds/beta

      ndimu=nint(nclo*1.1)

      r0=[x0,y0,z0]
      dr0=[xf0-x0,yf0-y0,zf0-z0]
      dr0=[efx,efy,efz]*perlen_u
      r0=r0+dr0/2.0d0

      allocate(frq(nepho_u),
     &  uampex(nepho_u),uampey(nepho_u),uampez(nepho_u),
     &  uampbx(nepho_u),uampby(nepho_u),uampbz(nepho_u),pow(nobsv,mthreads),
     &  utraxyz(14,ndimu),ustokes(4,nepho_u))

      pow=0.0d0
      frq=epho_u

      flow=frq(1)
      fhigh=frq(nepho_u)

      beff=sqrt(beffv_u**2+beffh_u**2)
      park=echarge1*beff*perlen_u/(2.*pi1*emasskg1*clight)
      wlen1=(1+park**2/2.)/2./gamma_u**2*perlen_u*1.0d9

      if (wlen1.ne.0.0) then
        eharm1=wtoe1/wlen1
      else
        eharm1=0.0d0
      endif

      dtpho=perlen_u/clight

      allocate(pherrc(nper_u),pherr(nper_u),arad(6,nepho_u*nobsv,mthreads),
     &  expphiran(max(1,nelec_u)))

      allocate(pranall(2,nelec_u))
      do i=1,nelec_u
        call util_random(3,pran)
        pranall(:,i)=pran(1:2)
        expphiran(i)=exp(dcmplx(0.0d0,twopi1*pran(3)))
      enddo

      if (ibunch.eq.0.or.
     &    emith_u.eq.0.0d0.and.emitv_u.eq.0.0d0.and.espread_u.eq.0.0d0) then
        iemit=0
        ibunch=0
        nelec=1
        nelec_u=1
      else
        iemit=1
      endif

      if (nelec_u.ne.nelec) then
        print*,''
        print*,'--- Warning in urad_amprep: Nelec adjusted to multiple of number of threads:',nelec_u
        print*,''
      endif

      if (iemit.ne.0) then
        allocate(eall(6,nelec_u))
        do i=1,nelec_u
          xi=x0
          if (modepin_u.ne.2) then
            call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &        emith_u,emitv_u,
     &        disph_u,dispph_u,dispv_u,disppv_u,
     &        espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
          else
            ! espread only for folding procedure
            call util_get_electron(xbeta_u,betah_u,alphah_u,betav_u,alphav_u,
     &        0.0d0,0.0d0,
     &        disph_u,dispph_u,dispv_u,disppv_u,
     &        espread_u,bunchlen_u,xi,yi,zi,ypi,zpi,dpp,modebunch_u)
          endif
          eall(1,i)=xi-x0
          eall(2,i)=yi
          eall(3,i)=zi
          eall(4,i)=ypi
          eall(5,i)=zpi
          eall(6,i)=dpp
        enddo
        if (noranone.ne.0) eall(:,1)=0.0
      endif

      !allocate(affe(6,nepho_u*nobsv))

      allocate(wsstokes(4,nepho_u*nobsv),stokes(4,nepho_u*nobsv,mthreads))
      stokes=0.0d0
      arad=(0.0d0,0.0d0)

      np2=nper_u/2

      call util_random_gauss_omp(nper_u,pherr,rr)
      pherrc=pherr

      lmodeph=modeph_u

      if (pherror_u.ne.0.0d0.and.(lmodeph.lt.0.or.lmodeph.gt.2)) then
        write(6,*) ""
        write(6,*) "*** Error in urad_amprep: MODEPH must be 0,1, or 2 ***"
        write(6,*) "*** Program aborted ***"
      endif

      if (lmodeph.eq.0.and.eharm1.ne.0.0d0) then
        om1=eharm1/hbarev
        pherr=sngl(pherrc*pherror_u/360.0d0*twopi1/om1)
      else if (lmodeph.eq.1) then
        pherr=sngl(pherr*pherror_u)
      else if (lmodeph.eq.2) then
        pherr(nper_u)=0.0
        phsum=0.0d0
        do i=1,nper_u-1
          pherr(i)=pherr(i)+pherrc(i)
          pherr(i+1)=pherr(i+1)-pherrc(i)
          phsum=phsum+pherr(i)
        enddo
        phsum=phsum+pherr(nper_u)
      else
        pherr=0.0d0
      endif !(lmodeph.eq.0)

      mbunch=max(1,nelec_u)
      nelec=nelec_u

      if (ibunch.ne.0.and.bunchcharge_u.ne.0.0d0) then
        sqnbunch=mbunch
        sqnphsp=sqrt(bunchcharge_u/echarge1)
     &    *meinbunch
     &    /(bunchcharge_u/echarge1)
        bunnor=1.0d0/mbunch
      else
        sqnbunch=mbunch
        sqnphsp=sqrt(dble(nelec_u))
        bunnor=1.0d0/mbunch
      endif

      beff=sqrt(beffv_u**2+beffh_u**2)
      parke=echarge1*beff*perlen_u/(2.*pi1*emasskg1*clight)
      xlell=perlen_u

      ielec=0

      pow=0.0d0
      noemit=0
      noespread=0
      jbunch=ibunch
      jubunch=0
      ebeam=ebeam_u
      debeam=espread_u
      stokesv=vstokes
      specnor=
     &  banwid_u
     &  /(4.0d0*pi1**2*clight*hbarev)
     &  /(4.0d0*pi1*eps01)
     &  *curr_u
      sbnor=specnor*bunnor
      speknor=specnor
      jeneloss=0
      pw=pinw_u
      ph=pinh_u
      !pr=pinr
      pc=pincen_u
      do iobsv=1,nobsv
        if (abs(obsv_u(2,iobsv)).lt.1.0d-9) obsv_u(2,iobsv)=0.0d0
        if (abs(obsv_u(3,iobsv)).lt.1.0d-9) obsv_u(3,iobsv)=0.0d0
      enddo
      pcbrill=obsv_u(:,icbrill)
      phgsh=phgshift_u

      rea=(0.0d0,0.0d0)
      expsh=(1.0d0,0.0d0)

      if (ifieldprop_u.eq.2) then
        if (npinzprop_u.eq.1) then
          zprop(1)=0.0d0
        else
          dz=pinwprop_u/(npinzprop_u-1)/1000.0d0
          zprop(1)=-pinwprop_u/2.0d0/1000.0d0
          do i=2,npinzprop_u
            zprop(i)=zprop(i-1)+dz
          enddo
        endif

        xprop=pinxprop_u

        if (npinyprop_u.eq.1) then
          yprop(1)=0.0d0
        else
          dy=pinhprop_u/(npinyprop_u-1)/1000.0d0
          yprop(1)=-pinhprop_u/2.0d0/1000.0d0
          do i=2,npinyprop_u
            yprop(i)=yprop(i-1)+dy
          enddo
        endif
      endif !ifieldprop

      ifieldprop=ifieldprop_u
      nyprop=npinyprop_u
      nzprop=npinzprop_u
      cjvsto=dconjg(vstokes)

      pinwprop=pinwprop_u
      pinhprop=pinhprop_u

      ifix=ifixphase_u

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& FIRSTPRIVATE(nepho,nobsvz,nobsvy,nobsv,nelec,frq,nper_u,np2,perlen_u,clight,hbarev,
!$OMP& ifieldprop,xprop,nyprop,nzprop,yprop,zprop,cjvsto,fpriv,pinwprop,pinhprop,
!$OMP& flow,fhigh,czero,cone,rea,expsh,ifix,zob,yob,
!$OMP& x0,y0,z0,xf0,yf0,zf0,vx0,vy0,vz0,vxf0,vyf0,vzf0,gamma_u,sbnor,speknor,
!$OMP& efx,efy,efz,ds,ndimu,curr_u,xlell,parke,amp,amp0,
!$OMP& uampex,uampey,uampez,uampbx,uampby,uampbz,
!$OMP& lmodeph,zp0,yp0,modewave,
!$OMP& jbunch,jubunch,jhbunch,noespread,noemit,ebeam,
!$OMP& stokesv,icbrill,obsv_u,emassg,debeam,dispv_u,disppv_u,
!$OMP& betah_u,alphah_u,betav_u,alphav_u,emith_u,emitv_u,disph_u,dispph_u,
!$OMP& pran,pranall,eall,fillb,r0,dr0,iamppin,iamppincirc,pc,phase0,pr,banwid_u,
!$OMP& pw,ph,idebug,pcbrill,wsstokes,vn,bunchlen_u,modebunch_u,icohere_u)
!$OMP& SHARED(mthreads,stokes,pherr,expphiran,lbunch,lnbunch,modepin_u,fieldbunch,npinzo_u,nobsvo,dzpin,dypin,
!$OMP& fbunch_u,jcharge,jeneloss,jvelofield,iemit,noranone,arad,pow,zmin,ymin,phgsh,fprop,stokesprop)

      jbun=1
      isub=0
      iobsv=0
      ielo=0
      xph0=-perlen_u*dble(nper_u)/2.0d0

!$OMP DO

c      do ilo=1,nelec*nobsv
      do ielec=1,nelec
      do iobsv=1,nobsv

        wsstokes=0.0d0

        !affe=(0.0D0,0.0D0)
        spow=0.0d0

        ith=OMP_GET_THREAD_NUM()+1

c        iobsv=mod(ilo-1,nobsv)+1
c        ibu=(ilo-1)/nobsv+1
        ibu=ielec
        jbun=ibu

        iy=(iobsv-1)/nobsvz+1
        iz=mod(iobsv-1,nobsvz)+1

        !if (iz.gt.nobsvz/2+1) call til_break

c        ielec=ibu

        xi=x0
        yi=y0
        zi=z0

        zpi=vz0/vx0
        ypi=vy0/vx0

        x2=xf0
        y2=yf0
        z2=zf0

        vx2=vxf0
        vy2=vyf0
        vz2=vzf0

        gamma=gamma_u

        dpp=0.0d0

        if (iemit.ne.0) then

          if (noranone.eq.0.or.ielec.ne.1) then

            bunchx=eall(1,ielec)

            xi=xi+bunchx
            yy=eall(2,ielec)
            zz=eall(3,ielec)

            yyp=eall(4,ielec)
            zzp=eall(5,ielec)

            dpp=eall(6,ielec)
            gamma=(1.0d0+dpp)*gamma_u

            ! assume beta(s)=beta0(s)+s**2/beta(0) and alpha0=-s/beta(0)
            ! and a drift transfer-matrix ((1,s),(1,0))

            zi=zz-x0*zzp !inverse transformation
            zpi=zzp

            yi=yy-x0*yyp
            ypi=yyp

            ! simple treatment of closed orbit, assume small angles

            zi=zi+z0
            zpi=zpi+zp0

            yi=yi+y0
            ypi=ypi+yp0

          else

            xi=x0
            yi=y0
            zi=z0
            ypi=yp0
            zpi=zp0

            bunchx=0.0d0

          endif

        else

          xi=x0
          yi=y0
          zi=z0
          ypi=yp0
          zpi=zp0

          bunchx=0.0d0

        endif !iemit

c+self,if=old.
c        zi=zi+dpp*di0
c        zpi=zpi+dpp*dd0
c+self.
        vn=clight*dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))

        vxi=vn/sqrt(1.0d0+ypi**2+zpi**2)
        vyi=vxi*ypi
        vzi=vxi*zpi

        obs=obsv_u(1:3,iobsv)

        if (noranone.eq.0.or.ielec.ne.1.or.iobsv.ne.icbrill) then
          if (iamppin.eq.3) then
            !call util_random(2,pran)
            pran(1:2)=pranall(:,ielec)
            if (iamppincirc.eq.0) then
              obs(2)=pc(2)+(pran(1)-0.5)*ph
              obs(3)=pc(3)+(pran(2)-0.5)*pw
            else
              rpin=(pran(1)-0.5)*pr
              ppin=pran(2)*twopi1
              obs(2)=pc(2)+rpin*cos(ppin)
              obs(3)=pc(3)+rpin*sin(ppin)
            endif
          endif
        endif

        vn=norm2([vxi,vyi,vzi])
        eix=vxi/vn
        eiy=vyi/vn
        eiz=vzi/vn

        h2=((obs(2)-yi)**2+(obs(3)-zi)**2)/(obs(1)-xph0)**2
        if (h2.lt.0.01) then
          rph=abs(obs(1)-xph0)*(1.0d0+(((((-0.0205078125D0*h2+0.02734375D0)*h2
     &      -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2)
        else
          rph=sqrt((obs(1)-xph0)**2+((obs(2)-yi)**2+(obs(3)-zi)**2))
        endif

        phase0=(rph-(obsv_u(1,icbrill)-xph0))/clight

        call urad_e_b_field(
     &    jcharge,curr_u,
     &    gamma,udgamtot,
     &    xi,yi,zi,vxi,vyi,vzi,
     &    xf0,yf0,zf0,efx,efy,efz,
     &    x2,y2,z2,vx2,vy2,vz2,dtelec,ds,
     &      0,nstepu,ndimu,utraxyz,phase0,
     &    obs(1),obs(2),obs(3),flow,fhigh,
     &    nepho,frq,uampex,uampey,uampez,uampbx,uampby,uampbz,
     &    ustokes,upow,
     &    jeneloss,jvelofield,ifail,ith,banwid_u,modewave
     &    )

        r0=[xi,yi,zi]
        dr0=[x2-xi,y2-yi,z2-zi]

        drn=dr0/norm2(dr0)
        r0=r0+dr0/2.0d0

        do kfreq=1,nepho

          iobfr=iobsv+nobsv*(kfreq-1)

          om=frq(kfreq)/hbarev

          if (modewave.eq.0) then
            amp0=[
     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
     &        ]*1.0d3/sqrt(speknor/curr_u*0.10d0) !urad
          else
            amp0=[
     &        uampex(kfreq),uampey(kfreq),uampez(kfreq),
     &        uampbx(kfreq),uampby(kfreq),uampbz(kfreq)
     &        ]*1.0d3/sqrt(speknor) !urad
          endif

c          call util_random(1,pran)
c          amp0=amp0*dcmplx(0.0d0,dble(pran(1)*twopi1))

          amp=(0.0d0,0.0d0)
          t=bunchx/vn

          do i=1,nper_u

            r=r0+(i-np2-1)*dr0
            dobs=obs-r
            dist0=norm2(obs-r0)
            dist=norm2(dobs)

            if (kfreq.eq.1) then
              spow=spow+upow*(dist0/dist)**2
              pow(iobsv,ith)=pow(iobsv,ith)+upow*(dist0/dist)**2
            endif

            if (lmodeph.eq.0) then
!!!!!                dt=xlell/clight*((1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
!!!!!     &            (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)
              h2=
     &          (ypi-yp0-dobs(2)/dobs(1))**2 +
     &          (zpi-zp0-dobs(3)/dobs(1))**2
c26.4.2024     &          ((ypi-yp0-dobs(2))/dobs(1))**2 +
c26.4.2024     &          ((zpi-zp0-dobs(3))/dobs(1))**2

c              dt=xlell/clight*
c     &          (
c     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+h2/2.0d0-h2**2/8.0d0
c     &          )

              dph=om*(t+pherr(i))

              dt=xlell/clight*
     &          (
     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
     &          (((((-0.0205078125D0*h2+0.02734375D0)*h2
     &          -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2
     &          )

              t=t+dt
            else if (lmodeph.eq.1.or.lmodeph.eq.2) then
              dph=om*t
              pkerr=parke*(1.0d0+pherr(i))
!!!!!                dt=xlell/clight*((1.0d0+pkerr**2/2.0d0)/2.0d0/gamma**2+
!!!!!     &            (((ypi-dobs(2)/dobs(1))**2+(zpi-dobs(3)/dobs(1))**2))/2.0d0)
              h2=
     &          (ypi-yp0-dobs(2)/dobs(1))**2 +
     &          (zpi-zp0-dobs(3)/dobs(1))**2
c26.4.2024              h2=((ypi-yp0-dobs(2))**2+(zpi-zp0-dobs(3))**2)/dobs(1)**2
              dt=xlell/clight*
     &          (
c25.4.2024     &          (1.0d0+parke**2/2.0d0)/2.0d0/gamma**2+
     &          (1.0d0+pkerr**2/2.0d0)/2.0d0/gamma**2+
     &          (((((-0.0205078125D0*h2+0.02734375D0)*h2
     &          -0.0390625D0)*h2+0.0625D0)*h2-0.125D0)*h2+0.5D0)*h2
     &          )
              t=t+dt
            endif !lmodeph

            zexp=cdexp(dcmplx(0.0d0,dph))
            damp=amp0*zexp*dist0/dist
            amp=amp+damp

            if (jhbunch.ne.0) then

              if (
     &            ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
     &            mod(ielec,jhbunch).eq.0) .or.
     &            (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

                if (i.eq.1) then
                  fillb(5)=r(1)
                  fillb(6)=r(2)
                  fillb(7)=r(3)
                  fillb(8)=ypi
                  fillb(9)=zpi
                else if (i.eq.nper_u) then

                  if (abs(phgsh).eq.9999.0d0) then
                    rea=amp(1:3)
                    expsh=rea(3)/abs(rea(3))
                    if (phgsh.eq.-9999.0d0) expsh=expsh*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
                    amp=amp/expsh
                  else if (phgsh.ne.0.0d0) then
                    expsh=cdexp(dcmplx(0.0d0,phgsh))*1.0d3
                    amp=amp/expsh
                  endif

                  if (phgsh.eq.9999.0d0) then
                    rea=amp(1:3)
                    expsh=rea(3)/abs(rea(3))
                    amp=amp/expsh
                  else if (phgsh.eq.-9999.0d0) then
                    rea=amp(1:3)
                    expsh=rea(3)/abs(rea(3))*cdexp(dcmplx(0.0d0,-pi1/2.0d0))
                    amp=amp/expsh
                  else if (phgsh.ne.0.0d0) then
                    expsh=cdexp(dcmplx(0.0d0,phgsh))*1.0d3
                    amp=amp/expsh
                  endif

                  fillb(10:12)=r
                  fillb(13)=ypi
                  fillb(14)=zpi
                  fillb(30)=dreal(amp(1))
                  fillb(31)=dimag(amp(1))
                  fillb(32)=dreal(amp(2))
                  fillb(33)=dimag(amp(2))
                  fillb(34)=dreal(amp(3))
                  fillb(35)=dimag(amp(3))
                  fillb(36)=dreal(amp(4))
                  fillb(37)=dimag(amp(4))
                  fillb(38)=dreal(amp(5))
                  fillb(39)=dimag(amp(5))
                  fillb(40)=dreal(amp(6))
                  fillb(41)=dimag(amp(6))
                endif

              endif

            endif

          enddo !nper_u

          if (ifix.eq.2) then
            amp=amp*expphiran(ielec)
          endif

          if (modepin_u.ne.0) then
            iy=int((obs(2)-ymin)/dypin)+1
            iz=int((obs(3)-zmin)/dzpin)+1
            !print*,ilo,ith,obs(3),zmin,dzpin,iz
            fieldbunch(1:6,iz,iy,kfreq)=fieldbunch(1:6,iz,iy,kfreq)+amp(1:6)
            fieldbunch(7,iz,iy,kfreq)=fieldbunch(7,iz,iy,kfreq)+cone
          endif

          apolh=
     &      amp(1)*conjg(stokesv(1,1))
     &      +amp(2)*conjg(stokesv(1,2))
     &      +amp(3)*conjg(stokesv(1,3))

          apolr=
     &      amp(1)*conjg(stokesv(2,1))
     &      +amp(2)*conjg(stokesv(2,2))
     &      +amp(3)*conjg(stokesv(2,3))

          apoll=
     &      amp(1)*conjg(stokesv(3,1))
     &      +amp(2)*conjg(stokesv(3,2))
     &      +amp(3)*conjg(stokesv(3,3))

          apol45=
     &      amp(1)*conjg(stokesv(4,1))
     &      +amp(2)*conjg(stokesv(4,2))
     &      +amp(3)*conjg(stokesv(4,3))

          stok1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
          stok2=dreal(-stok1+2.0d0*apolh*conjg(apolh))
          stok3=dreal(2.0d0*apol45*conjg(apol45)-stok1)
          stok4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

          wsstokes(1,iobfr)=wsstokes(1,iobfr)+stok1*sbnor
          wsstokes(2,iobfr)=wsstokes(2,iobfr)+stok2*sbnor
          wsstokes(3,iobfr)=wsstokes(3,iobfr)+stok3*sbnor
          wsstokes(4,iobfr)=wsstokes(4,iobfr)+stok4*sbnor

          stokes(1:4,iobfr,ith)=stokes(1:4,iobfr,ith)+wsstokes(1:4,iobfr)

          !affe(:,iobfr)=affe(:,iobfr)+amp
          !arad(:,iobfr,ith)=arad(:,iobfr,ith)+affe(:,iobfr)
          arad(:,iobfr,ith)=arad(:,iobfr,ith)+amp

c          if (
c     &        ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
c     &        mod(ielec,jhbunch).eq.0) .or.
c     &        (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then
c            if (ielec.gt.4) then
c              print*,jhbunch,ith,ilo,ielec
c            endif
c          endif

          if (jhbunch.ne.0) then

            if (
     &          ((iamppin.eq.3.or.iobsv.eq.icbrill).and.jhbunch.gt.0.and.
     &          mod(ielec,jhbunch).eq.0) .or.
     &          (jhbunch.lt.0.and.mod(ielec,-jhbunch).eq.0)) then

c              print*,jhbunch,ith,ilo,jbun,isub,ibu
              fillb(1)=jbun
              fillb(2)=isub
              fillb(3)=ibu
              fillb(4)=bunchx
              fillb(15)=gamma*emassg
              fillb(16)=udgamtot*emassg
              fillb(17)=obs(1)
              fillb(18)=obs(2)
              fillb(19)=obs(3)
              fillb(20)=kfreq
              fillb(21)=frq(kfreq)

              fillb(22)=wsstokes(1,iobfr)*nelec

              fillb(23)=wsstokes(1,iobfr)*nelec
              fillb(24)=wsstokes(2,iobfr)*nelec
              fillb(25)=wsstokes(3,iobfr)*nelec
              fillb(26)=wsstokes(4,iobfr)*nelec

              fillb(27)=spow
              fillb(28)=1
              fillb(29)=dtelec

              fillb(30)=dreal(amp(1))
              fillb(31)=dimag(amp(1))
              fillb(32)=dreal(amp(2))
              fillb(33)=dimag(amp(2))
              fillb(34)=dreal(amp(3))
              fillb(35)=dimag(amp(3))
              fillb(36)=dreal(amp(4))
              fillb(37)=dimag(amp(4))
              fillb(38)=dreal(amp(5))
              fillb(39)=dimag(amp(5))
              fillb(40)=dreal(amp(6))
              fillb(41)=dimag(amp(6))
              lbunch=lbunch+1
              fbunch_u(:,lbunch)=fillb(:)
            endif !fill

          endif !jhbunch

          if (ifieldprop.eq.2) then
            if (iobsv.eq.1) then
              fprop(1:3,1:nzprop,1:nyprop,kfreq,ith)=(0.0d0,0.0d0)
            endif
            call urad_phase_prop_point(obs,amp(1:3),nzprop,nyprop,
     &        xprop,yprop,zprop,pinwprop,pinhprop,frq(kfreq),fpriv)
            fprop(:,:,:,kfreq,ith)=fprop(:,:,:,kfreq,ith)+fpriv(:,:,:)
            if (iobsv.eq.nobsv) then
              i=0
              do ipy=1,nyprop
                do ipz=1,nzprop
                  i=i+1+nzprop*nyprop*(kfreq-1)
                  apolh=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(1,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(1,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(1,3)

                  apolr=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(2,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(2,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(2,3)

                  apoll=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(3,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(3,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(3,3)

                  apol45=
     &              fprop(1,ipz,ipy,kfreq,ith)*cjvsto(4,1)+
     &              fprop(2,ipz,ipy,kfreq,ith)*cjvsto(4,2)+
     &              fprop(3,ipz,ipy,kfreq,ith)*cjvsto(4,3)

                  stok1=dreal(
     &              apolr*conjg(apolr)+
     &              apoll*conjg(apoll))

                  stok2=-stok1+
     &              dreal(2.*apolh*conjg(apolh))

                  stok3=
     &              dreal(2.*apol45*conjg(apol45))-
     &              stok1

                  stok4=dreal(
     &              apolr*conjg(apolr)-
     &              apoll*conjg(apoll))

                  stokesprop(1,ipz,ipy,kfreq,ith)=stokesprop(1,ipz,ipy,kfreq,ith)+stok1
                  stokesprop(2,ipz,ipy,kfreq,ith)=stokesprop(2,ipz,ipy,kfreq,ith)+stok2
                  stokesprop(3,ipz,ipy,kfreq,ith)=stokesprop(3,ipz,ipy,kfreq,ith)+stok3
                  stokesprop(4,ipz,ipy,kfreq,ith)=stokesprop(4,ipz,ipy,kfreq,ith)+stok4

                enddo
              enddo
            endif
          endif

        enddo !kfreq

      enddo !iobsv
      enddo !nelec
c      enddo !ilo

!$OMP END DO
!$OMP END PARALLEL

      do ith=1,mthreads
        pow_u(:)=pow_u(:)+pow(:,ith)
        arad_u(:,:)=arad_u(:,:)+arad(:,:,ith)
      enddo
      pow_u=pow_u/sqnbunch

      if (ifieldprop_u.eq.2) then
        smax=0.0d0
        do ith=1,mthreads
          do kfreq=1,nepho_u
            i=0
            do iy=1,nyprop
              do iz=1,nzprop
                i=i+1+nobsvprop_u*(kfreq-1)
                stokesprop_u(1:4,i)=stokesprop_u(1:4,i)+stokesprop(1:4,iz,iy,kfreq,ith)
                if(abs(stokesprop_u(1,i)).gt.smax) then
                  smax=abs(stokesprop_u(1,i))
                  izm=iz
                  iym=iy
                  im=i
                endif
              enddo
            enddo
          enddo
        enddo
      endif

      if (icohere_u.eq.0) then

        arad_u=arad_u/sqnbunch

        do ith=1,mthreads
          stokes_u(:,:)=stokes_u(:,:)+stokes(:,:,ith)
        enddo

      else

        do iobsv=1,nobsv
          do kfreq=1,nepho

            iobfr=iobsv+nobsv*(kfreq-1)

            amp(1:3)=arad_u(1:3,iobfr) !/sqnphsp

            apolh=
     &        amp(1)*conjg(stokesv(1,1))
     &        +amp(2)*conjg(stokesv(1,2))
     &        +amp(3)*conjg(stokesv(1,3))

            apolr=
     &        amp(1)*conjg(stokesv(2,1))
     &        +amp(2)*conjg(stokesv(2,2))
     &        +amp(3)*conjg(stokesv(2,3))

            apoll=
     &        amp(1)*conjg(stokesv(3,1))
     &        +amp(2)*conjg(stokesv(3,2))
     &        +amp(3)*conjg(stokesv(3,3))

            apol45=
     &        amp(1)*conjg(stokesv(4,1))
     &        +amp(2)*conjg(stokesv(4,2))
     &        +amp(3)*conjg(stokesv(4,3))

            stok1=dreal(apolr*conjg(apolr)+apoll*conjg(apoll))
            stok2=dreal(-stok1+2.0d0*apolh*conjg(apolh))
            stok3=dreal(2.0d0*apol45*conjg(apol45)-stok1)
            stok4=dreal(apolr*conjg(apolr)-apoll*conjg(apoll))

            stokes_u(1,iobfr)=stok1*sbnor
            stokes_u(2,iobfr)=stok2*sbnor
            stokes_u(3,iobfr)=stok3*sbnor
            stokes_u(4,iobfr)=stok4*sbnor

          enddo
        enddo

      endif !icohere_u

c      if (ihbunch.ne.0) then
c        n=0
c        do i=1,nlbu
c          do ith=1,mthreads_u
c            if (fbunch(21,i,ith).ne.0.0d0) then
c              n=n+1
c              fbunch_u(:,n)=fbunch(:,i,ith)
c            endif
c          enddo
c        enddo
c        deallocate(fbunch)
c      endif

      !deallocate(affe)
      deallocate(frq,uampex,uampey,uampez,uampbx,uampby,uampbz,utraxyz,
     &  pherrc,pherr,expphiran,arad,pow,pranall,wsstokes,stokes)

      if (iemit.ne.0) deallocate(eall)

      iobfr=nobsv_u*nepho_u/2+1
      amp(1:3)=arad_u(1:3,iobfr)

      anor=sqrt(stokes_u(1,iobfr)/
     &  (amp(1)*dconjg(amp(1))+amp(2)*dconjg(amp(2))+amp(3)*dconjg(amp(3))))
      arad_u=arad_u*anor

      if (modepin_u.ne.0) then
        do iepho=1,nepho_u
          do iy=1,npinyo_u-1
            do iz=1,npinzo_u-1
              fsum=max(1.0d0,dreal(fieldbunch(7,iz,iy,iepho)))
              fieldbunch(1:6,iz,iy,iepho)=fieldbunch(1:6,iz,iy,iepho)/fsum
c            fieldbunch(1:6,iz,iy,iepho)=fieldbunch(1:6,iz,iy,iepho)
            enddo
          enddo
        enddo
      endif

      return
      end
