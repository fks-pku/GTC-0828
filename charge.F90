subroutine charge(species_name)
  use global_parameters,only: ihybrid,nhybrid
  use particle_array
  use interfaces,only: gkChargeParticle,hybridChargeParticle
  use field_array,only: sfluidne
  implicit none

  integer ierror
  character(*),intent(in) :: species_name

  if(species_name=='thermal-ion')then
    !ideal MHD
    if(iload==0)then
      densityi=0.0
      flowi=0.0
    else
      call gkChargeParticle(zion,wpion,wtion0,wtion1,densityi,flowi,jtion0,&
        jtion1,wzion,zonali,zonalci,markeri,markerit,pmarki,qion,aion,&
        iload,ngyroi,mi,switch="w/o density modification")
    endif
  elseif(species_name=='fast-ion')then
    call gkChargeParticle(zfast,wpfast,wtfast0,wtfast1,densityf,flowf,&
      jtfast0,jtfast1,wzfast,zonalf,zonalcf,markerf,markerft,pmarkf,&
      qfast,afast,fload,ngyrof,mf)
  elseif(species_name=='fast-electron')then
    call gkChargeParticle(zfaste,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,&
      jtfaste0,jtfaste1,wzfaste,zonalfe,zonalcfe,markerfe,markerfet,pmarkfe,&
      qfaste,afaste,feload,ngyrofe,mfe)
  elseif(species_name=='thermal-electron')then
    call hybridChargeParticle(zelectron,wpelectron,wtelectron0,wtelectron1,&
      densitye,flowe,jtelectron0,jtelectron1,wzelectron,zonale,zonalce,&
      markere,markeret,pmarke,qelectron,aelectron,eload,ngyroe,me,&
      ihybrid=ihybrid,nhybrid=nhybrid,pressurepara=pressureepara,&
      pressureperp=pressureeperp,sfluidn=sfluidne,dnsave=dnesave)
  else
    write(*,*)'push.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine charge

subroutine locate(species_name)
  use particle_array
  use interfaces,only: locateParticle
  implicit none

  integer ierror
  character(*),intent(in) :: species_name

  if(species_name=='thermal-ion')then
    call locateParticle(zion,wzion,wpion,wtion0,wtion1,jtion0,jtion1,qion,&
      aion,ngyroi,mi)
  elseif(species_name=='thermal-electron')then
    call locateParticle(zelectron,wzelectron,wpelectron,wtelectron0,&
      wtelectron1,jtelectron0,jtelectron1,qelectron,aelectron,ngyroe,me)
  elseif(species_name=='fast-ion')then
    call locateParticle(zfast,wzfast,wpfast,wtfast0,wtfast1,jtfast0,jtfast1,&
      qfast,afast,ngyrof,mf)
  elseif(species_name=='fast-electron')then
    call locateParticle(zfaste,wzfaste,wpfaste,wtfaste0,wtfaste1,jtfaste0,&
      jtfaste1,qfaste,afaste,ngyrofe,mfe)
  else
    write(*,*)'push.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine locate

! find ion location for interpolation in both gathering (chargei.F90) and scattering (pushi.F90)
subroutine locateParticle(zpart,wzpart,wppart,wtpart0,wtpart1,jtpart0,&
    jtpart1,qpart,apart,ngyro,mp)
  use precision
  use global_parameters,only: rho0,mpsilow,rg0,psi0,psi1,pi2,fielddir,mpsi,&
    mype,istep
  use field_array,only: deltar,deltat,deltap,deltaz,psimesh,zeta0,mtheta,&
    igrid,pgyro,pgyro2,tgyro,tgyro2,qmesh,qtinv,solvermethod
  use equilibrium,only: lsp,spdpsi_inv,spdpsi,rgpsi
  implicit none

  !declaration of the dummy arguments
  real(wp) qpart,apart
  real(wp),dimension(:) :: wzpart
  real(wp),dimension(:,:) :: wppart,wtpart0,wtpart1
  real(wp),dimension(:,:) :: zpart
  integer,dimension(:,:) :: jtpart0,jtpart1
  integer ngyro,mp

  integer m,igyro,ip,jt,ij0,ii,im,j00,j01,isp
  real(wp) pdum0,tdum0,zdum,rho,rg,tflr,pdum,tdum,dpx,delr,delt(0:mpsi),delp(mpsi),delz,rhoc,xdum,ydum,sqrtq0_inv,paxis
#ifndef GPU_UM
  !$acc declare create(delt(0:mpsi),delp(mpsi))
#endif

  delr=1.0/deltar
  delt=1.0/deltat
  delz=1.0/deltaz
  delp=1.0/deltap
#ifndef GPU_UM
  !$acc update device(delt(0:mpsi),delp(mpsi))
#endif
  rhoc=sqrt(apart)/(rho0*qpart) !apart & qpart are needed to calculate ion gyroradius
  if(ngyro==1)rhoc=0.0

  sqrtq0_inv=1.0_wp/sqrt(qmesh(0))
  paxis=0.0
  if(psi0==0.0_wp)paxis=12.5*rho0*rho0*sqrtq0_inv!Use x-y locate for r < 5 rho0


#ifdef _OPENACC
  !$acc parallel loop gang vector
#else
  !$omp parallel do private(m,igyro,pdum0,tdum0,zdum,rho,isp,dpx,&
  !$omp& rg,ip,jt,ij0,pdum,ii,tflr,im,tdum,j00,j01)
#endif
  do m=1,mp
     pdum0=zpart(1,m)
     tdum0=zpart(2,m)
     zdum=zpart(3,m)
     rho=zpart(6,m)*rhoc
     wzpart(m)=(zdum-zeta0)*delz  !weight for upper toroidal grid

! Guiding Center location
     isp=max(1,min(lsp-1,ceiling(pdum0*spdpsi_inv)))
     dpx=pdum0-spdpsi*real(isp-1)
! radaial spline of rg avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dpx*dpx
     ip=max(0,min(mpsi,floor((rg-rg0)*delr+0.5)))       !radial grid closest to GC
     jt=max(0,min(mtheta(ip),floor(tdum0*delt(ip)+0.5))) !poloidal grid closest to GC
     ij0=igrid(ip)+jt                                   !GC grid index in magnetic coordinates

     do igyro=1,ngyro
!        pdum=max(psimesh(0),min(psimesh(mpsi),pdum0+rho*pgyro(igyro,ij0))) !particle position
        pdum=pdum0+rho*pgyro(igyro,ij0)+rho*rho*pgyro2(igyro,ij0) !particle position
        !if(pdum<0)print *,"pdum,bounded pdum:",pdum,modulo(pdum-psi0,psi1-psi0)+psi0
        !pdum=modulo(pdum-psi0,psi1-psi0)+psi0 !periodic BC
        if(pdum < psi0)pdum=2.0*psi0-pdum !reflective BC
        if(pdum > psi1)pdum=2.0*psi1-pdum
! particle position in theta
        tflr=tdum0+rho*tgyro(igyro,ij0)+rho*rho*tgyro2(igyro,ij0)
! Near axis gyroaveraging
        if(pdum0<paxis)then
          pdum = sqrt(pdum0)
          xdum = pdum*cos(tdum0)+rho*rho0*sqrtq0_inv*cos(tdum0+pi2*real(igyro-1)/4.0) 
          ydum = pdum*sin(tdum0)+rho*rho0*sqrtq0_inv*sin(tdum0+pi2*real(igyro-1)/4.0) 
          pdum = max(1.0e-8_wp*psi1,xdum*xdum+ydum*ydum)
          tflr = sign(1.0_wp,ydum)*acos(max(-1.0_wp,min(1.0_wp,xdum/&
            sqrt(zpart(1,m)))))
        endif
        isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
        dpx=pdum-spdpsi*real(isp-1)
! radaial spline of rg avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
        if(isp==1)dpx=sqrt(dpx)
        rg=rgpsi(1,isp)+rgpsi(2,isp)*dpx+rgpsi(3,isp)*dpx*dpx
        ii=max(0,min(mpsi-1,floor((rg-rg0)*delr)))    !radial grid on inner flux surface
        !solvermethod for fluxtube simulations:
        !0:fluxtube
        !1:semispectral
#ifdef _FRC
        if(solvermethod==0)then
        ! assuming mpsi=2, for fluxtube, then flux surfaces = [0,1,2]
           if(ii<1)then
             wppart(igyro,m)=1.0 !all weight on outer flux surface between [0,1]
           else
             wppart(igyro,m)=0.0 !all weight on inner flux surface between [1,2]
           endif
        else
           wppart(igyro,m)=(pdum-psimesh(ii))*delp(ii+1) !weight for outer flux surface
        endif
#else
        wppart(igyro,m)=(pdum-psimesh(ii))*delp(ii+1) !weight for outer flux surface
#endif
! inner flux surface
        im=ii
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(tflr-(zdum-pi2)*qtinv(im),pi2)*delt(im)
        else
          tdum=modulo(tflr-zdum*qtinv(im),pi2)*delt(im)
        endif
        j00=max(0,min(mtheta(im)-1,floor(tdum))) 
        jtpart0(igyro,m)=igrid(im)+j00  !lower poloidal grid on inner flux surface
        wtpart0(igyro,m)=tdum-real(j00) !weight for upper poloidal grid on inner flux surface
        if(ii == 0 .and. psi0== 0.0)wtpart0(igyro,m)=1.0
! outer flux surface
        im=ii+1
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(tflr-(zdum-pi2)*qtinv(im),pi2)*delt(im)
        else
          tdum=modulo(tflr-zdum*qtinv(im),pi2)*delt(im)
        endif
        j01=max(0,min(mtheta(im)-1,floor(tdum)))
        jtpart1(igyro,m)=igrid(im)+j01  !lower poloidal grid on outer flux surface
        wtpart1(igyro,m)=tdum-real(j01) !weight for upper poloidal grid on outer flux surface
     enddo
  enddo
  !$acc end parallel
end subroutine locateParticle
