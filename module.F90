module system_env
  use mpi
  use omp_lib
#ifdef _OPENACC
  use openacc
  use cudafor
#endif
end module system_env

module precision
  use system_env
  integer, parameter :: doubleprec=selected_real_kind(12),&
    singleprec=selected_real_kind(6),defaultprec=kind(0.0)
#ifdef DOUBLE_PRECISION
  integer, parameter :: wp=doubleprec,mpi_Rsize=MPI_DOUBLE_PRECISION,&
                        mpi_Csize=MPI_DOUBLE_COMPLEX
#else
  integer, parameter :: wp=singleprec,mpi_Rsize=MPI_REAL,&
                        mpi_Csize=MPI_COMPLEX
#endif
  real,parameter:: machineEpsilon=10.0_wp*epsilon(1.0_wp)
  save
end module precision

module magnetic_island
  use precision   
  integer,parameter :: l=1 !island number  
  integer wi,wo,ires,qs
  integer,dimension(:),allocatable :: isn,ism
  real(wp),dimension(:,:),allocatable :: hmesh_total,hmesh_perturb,hangle,alphaIs !helical flux , alpha_zero, ksi in Eq(15), alpha in Eq(11)
  real(wp),dimension(:,:,:),allocatable :: gradalphaIs !gradient of alphais
#ifndef GPU_UM
  !$acc declare create(gradalphaIs,alphaIs)
#endif
  save
end module magnetic_island

module global_parameters
  use precision
  integer,parameter :: gtcout=11, num_modes=8

! control parameters
  integer track_particles,ndata3d,mgrid,mpsi,mthetamax,mtoroidal,iodiag,iodata1d,&
       nspecies,istep,ndiag,msnap,mstep,mstepall,izonal,nbound,nboundR,irun,irk,idiag,&
       ncyclee,ncyclefe,mtdiag,nhybrid,ihybrid,nparami,magnetic,nonlinear,nfilter,numereq,&
       neop,neot,neoz,eqcurrent,fielddir,nparamf,nparamfe,nparame,mpsilow,mgridlow,&
       mpsihigh,mgridhigh,bcond,fieldmodel,ilaplacian,antenna,tearingmode,myth0,&
       n_modes(num_modes),m_modes(num_modes),island,fem,mgrid_fem,trilength,ismooth,inorm,irestore,&
       irotation,idiagonal
  real(wp) paranl,psi0,psi1,rg0,pi,pi2,pi2_inv,tstep,ulength,utime,rho0,maxwell(100001),&
       r0,b0,etemp0,eden0,hypr1,hypr2,omega_antenna,eta,toroidaln
#ifdef _FRC
  real(16) torbound
#else
  real(wp) torbound
#endif
#ifndef GPU_UM
  !$acc declare create(nparami,nparame,nparamf,nparamfe)
#endif
! MPI toroidal and particle decomposion
  integer mype,numberpe,npartdom,toroidal_comm,partd_comm,nproc_partd,myrank_partd,&
       nproc_toroidal,myrank_toroidal,left_pe,right_pe,&
       toroidal_domain_location,particle_domain_location

#ifdef _OPENMP
  integer nthreads
#endif

!!XY rolling restart
  integer irest,FileExit
  character(len=10) restart_dir1,restart_dir2
  save
end module global_parameters

module equilibrium
  integer lsp,lst
#ifdef _TOROIDAL3D
  integer,parameter :: spdim=27 ! non-axisymmetric equilibrium
#else
  integer,parameter :: spdim=9
#endif
  real psiw,ped,spdpsi,spdtheta,spdrg,spdtor,spdpsi_inv,spdtheta_inv,spdrg_inv,spdtor_inv
  real,dimension(:),allocatable :: stpp,mipp,mapp,mesher
  real,dimension(:,:),allocatable :: qpsi,gpsi,ppsi,rpsi,torpsi,tpsi,npsi,nepp,tepp,&
       tipp,nipp,tfpp,nfpp,tfepp,nfepp,zepp,ropp,erpp,spcos,spsin,rgpsi,psirg,psitor,cpsi,xygrid
  real,dimension(:,:,:),allocatable :: bsp,xsp,zsp,gsp,fsp,rd,nu,dl,ha,hb
#ifndef GPU_UM
  !$acc declare create(qpsi,gpsi,ropp,erpp,rgpsi,cpsi)
  !$acc declare create(bsp,xsp,mesher)
#endif
  save
end module equilibrium

module particle_array
  use precision
! particle diagnostics: # of quantities per species in history.out and data1d.out
  integer,parameter :: mpdiag=10,mpdata1d=3

! electron
  integer me,me1,memax,ngyroe,eload,etrap,ecoll
  real(wp) qelectron,aelectron,betae,tauee,tauei
  real(wp),dimension(mpdiag) :: diagelectron
  integer,dimension(:,:),allocatable :: jtelectron0,jtelectron1
  real(wp),dimension(:),allocatable :: wzelectron,meshte,meshne,kapane,kapate,&
    dtndpsi,pmarke,zonale,zonalce,markere,rdteme,pfluxe,markeret
  real(wp),dimension(:,:),allocatable :: zelectron,zelectron0,zelectron1,&
    wpelectron,wtelectron0,wtelectron1,densitye,flowe,phit,data1de,dnet,&
    pressureepara,pressureeperp
#ifndef GPU_UM
  !$acc declare create(diagelectron,jtelectron0,jtelectron1)
  !$acc declare create(wzelectron,meshte,meshne,kapane,kapate,rdteme)
  !$acc declare create(meshte,meshne,kapane,kapate)
  !$acc declare create(zelectron,zelectron0,zelectron1,wpelectron,wtelectron0)
  !$acc declare create(wtelectron1,densitye,flowe,data1de,phit,dnet,pressureepara,pressureeperp)
#endif
  real(wp),dimension(:,:,:),allocatable :: phisave,dnesave
  real(wp) tfracn
  real(wp),dimension(:),allocatable :: tfrac

! themal ion (main ion)
  integer mi,mimax,ngyroi,iload,icoll
  real(wp) qion,aion,tauii
  real(wp),dimension(mpdiag) :: diagion
  integer,dimension(:,:),allocatable :: jtion0,jtion1
  real(wp),dimension(:),allocatable :: wzion,meshti,meshni,kapani,kapati,&
       jacobianpsi,pmarki,zonali,zonalci,markeri,rdtemi,pfluxi,markerit
  real(wp),dimension(:,:),allocatable :: zion,zion0,wpion,wtion0,wtion1,&
    densityi,flowi,data1di,pressureipara,pressureiperp
#ifndef GPU_UM
  !$acc declare create(diagion,jtion0,jtion1)
  !$acc declare create(wzion,meshti,meshni,kapani,kapati,rdtemi)
  !$acc declare create(zion,zion0,wpion,wtion0,wtion1,densityi,flowi,data1di,pressureipara,pressureiperp)
#endif

! fast ion
  integer mf,mfmax,ngyrof,fload
  real(wp) qfast,afast,sd_v0,sd_vc,sd_l0,sd_widthInv
  real(wp),dimension(mpdiag) :: diagfast
  integer,dimension(:,:),allocatable :: jtfast0,jtfast1
  real(wp),dimension(:),allocatable :: wzfast,meshtf,meshnf,kapanf,kapatf,&
       pmarkf,zonalf,zonalcf,markerf,rdtemf,pfluxf,markerft
  real(wp),dimension(:,:),allocatable :: zfast,zfast0,wpfast,wtfast0,wtfast1,densityf,flowf,data1df
#ifndef GPU_UM
  !$acc declare create(diagfast,jtfast0,jtfast1)
  !$acc declare create(wzfast,meshtf,meshnf,kapanf,kapatf,rdtemf)
  !$acc declare create(zfast,zfast0,wpfast,wtfast0,wtfast1,densityf,flowf,data1df)
#endif

!fast electron
  integer mfe,mfe1,mfemax,ngyrofe,feload,fetrap
  real(wp) qfaste,afaste
  real(wp),dimension(mpdiag) :: diagfaste
  integer,dimension(:,:),allocatable :: jtfaste0,jtfaste1
  real(wp),dimension(:),allocatable :: wzfaste,meshtfe,meshnfe,kapanfe,kapatfe,&
       pmarkfe,zonalfe,zonalcfe,markerfe,rdtemfe,pfluxfe,markerfet
  real(wp),dimension(:,:),allocatable :: zfaste,zfaste0,zfaste1,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,data1dfe
  real(wp) fetfracn
  real(wp),dimension(:),allocatable :: fetfrac
#ifndef GPU_UM
  !$acc declare create(diagfaste,jtfaste0,jtfaste1)
  !$acc declare create(wzfaste,meshtfe,meshnfe,kapanfe,kapatfe,rdtemfe)
  !$acc declare create(zfaste,zfaste0,zfaste1,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,data1dfe)
#endif
  save
end module particle_array

module particle_tracking
  use precision
  real(wp),dimension(:,:),allocatable :: ptrackedi,ptrackede,ptrackedf
  integer,dimension(3) :: ntrackp
  save
end module particle_tracking

module field_array
  use precision
! PIC global fieldline following mesh
  real(wp) deltar,deltaz,zeta1,zeta0
  integer,dimension(:),allocatable :: itran,igrid,mtheta,nshift,igrid_fem
  real(wp),dimension(:),allocatable :: deltap,deltat,qtinv,qmesh,bmesh,jmesh,psimesh,tormesh,meshze
#ifndef GPU_UM
  !$acc declare create(igrid,mtheta)
  !$acc declare create(deltap,deltat,qtinv,psimesh)
#endif
  real(wp),dimension(:),allocatable :: gpsi200,b2m00,wtshift
  real(wp),dimension(:),allocatable :: gupp,gupt,gutt,guzz,gupz,gutz,rhom
  real(wp),dimension(:),allocatable :: gdpp,gdpt,gdtt,gdzz,gdpz,gdtz
  real(wp),dimension(:),allocatable :: spectraln

! fields on mesh: phi, apara, fluidne, fluidue, zonal and gradient
  real(wp),dimension(:),allocatable :: phi00,phip00,apara00,apara00nl,apara00nl0,fluidne00,d4fluidne,d2apara
  real(wp),dimension(:,:),allocatable :: phi,apara,fluidne,fluidue,apara0,fluidne0,&
                                         deltapsi,deltapsi0,sdeltapsi,&
                                         sapara,sfluidne,sdelapara,MHDxi_mag,phi_zero
  real(wp),dimension(:,:,:),allocatable :: gradphi,gradapara,gradue,gradne,gradext,gradpsi,gradphieff,&
                                           gradkne,gradpepara,gradpeperp,gradgaugef,MHDxiperp
#ifndef GPU_UM
  !$acc declare create(phip00,sapara)
  !$acc declare create(gradphi,gradapara,gradpsi,gradphieff)
#endif

! external field
  real(wp),dimension(:),allocatable :: omega
  real(wp),dimension(:,:),allocatable :: phiext,dn_ext

! diagnostics and filtering
  integer iflux,modes,solvermethod
  integer,dimension(:),allocatable :: nmodes,mmodes
  integer nmode

! radial interpolation
  integer,dimension(:,:),allocatable :: jtp1,jtp2
  real(wp),dimension(:,:),allocatable :: wtp1,wtp2

! laplacian coefficients
  integer mindexlap,mindexlap2,mindex_fem
  integer,dimension(:),allocatable :: nindexlap,nindexlap2,nindex_fem
  integer,dimension(:,:),allocatable :: indexplap,indexlap2,indexp_fem
  integer,dimension(:,:),allocatable :: trilist
  real(wp),dimension(:,:),allocatable :: lapmat,lapmat2,lmat,dmat

! theta of max and min b-field for particle baoundary cond.
  real(wp) maxbfield(0:1),minbfield(0:1),thetabmin(0:1),thetabmax(0:1)
  integer,dimension(:,:),allocatable :: thetaupp,thetadwn

! gyro averaging
  real(wp),dimension(:,:),allocatable :: pgyro,tgyro,pgyro2,tgyro2
#ifndef GPU_UM
  !$acc declare create(pgyro,tgyro,pgyro2,tgyro2)
#endif
  save
end module field_array

module petsc_array
  use precision
  integer newcomm,nproc_newcomm,myrank_newcomm
  integer,dimension(:),allocatable :: userp,users,luserp,lusers,luserp2,lusers2
  real(wp),dimension(:),allocatable :: usera,userb,userx,lusera,luserb,luserx,lusera2,luserb2,luserx2
  save
end module petsc_array

module data_type
  integer,parameter :: bp_char=0,bp_short=1,&
                bp_int=2,bp_long=3,bp_longlong=4,bp_float=5,&
                bp_double=6,bp_longdouble=7,bp_pointer=8,&
                bp_string=9,bp_uchar=50,bp_ushort=51,bp_uint=52,&
                bp_ulong=53,bp_ulonglong=54
  save
end module data_type

module spline_normalization
  contains
    subroutine normalizeSpline2d(spline,physicalUnit,fluxUnit)
      implicit none
  
      ! arguments
      real,dimension(:,:,:),intent(inout) :: spline
      real,intent(in) :: physicalUnit,fluxUnit
      ! internal variables
      integer :: i,pow
      real :: a

      spline=spline/physicalUnit

      do i = 1,9
        a = 1.0
        do pow = 1,modulo(i-1,3)
          a = a*fluxUnit
        enddo
        spline(i,:,:)=spline(i,:,:)*a
        spline(i,1,:)=spline(i,1,:)/sqrt(a)
      enddo
    end subroutine normalizeSpline2d

    subroutine normalizeSpline1d(spline,physicalUnit,fluxUnit)
      implicit none
  
      ! arguments
      real,dimension(:,:),intent(inout) :: spline
      real,intent(in) :: physicalUnit,fluxUnit
      ! internal variables
      integer :: i,pow
      real :: a
  
      spline=spline/physicalUnit

      do i = 1,3
        a = 1.0
        do pow = 1,i-1
          a = a*fluxUnit
        enddo
        spline(i,:)=spline(i,:)*a
      enddo
    end subroutine normalizeSpline1d
end module spline_normalization

module interfaces
  interface
    subroutine gkChargeParticle(zpart,wppart,wtpart0,wtpart1,density,&
        flow,jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,&
        qpart,apart,pload,ngyro,mp,ihybrid,nhybrid,&
        pressurepara,pressureperp,sfluidn,dnsave,&
        switch)
      use precision
      implicit none

      integer pload,ngyro,mp
      integer,optional :: ihybrid,nhybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart,marker,markert
      real(wp),dimension(0:) :: zonal,zonalc,pmark
      real(wp),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: density,flow
      real(wp),dimension(:,:) :: zpart
      real(wp),dimension(:,:),optional :: pressurepara,pressureperp,sfluidn
      real(wp),dimension(:,:,:),optional :: dnsave
      character(*),intent(in),optional :: switch
    end subroutine gkChargeParticle

    subroutine hybridChargeParticle(zpart,wppart,wtpart0,wtpart1,density,&
        flow,jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,&
        qpart,apart,pload,ngyro,mp,ihybrid,nhybrid,&
        pressurepara,pressureperp,sfluidn,dnsave,&
        switch)
      use precision
      implicit none

      integer pload,ngyro,mp
      integer,optional :: ihybrid,nhybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart,marker,markert
      real(wp),dimension(0:) :: zonal,zonalc,pmark
      real(wp),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: density,flow
      real(wp),dimension(:,:) :: zpart
      real(wp),dimension(:,:),optional :: pressurepara,pressureperp,sfluidn
      real(wp),dimension(:,:,:),optional :: dnsave
      character(*),intent(in),optional :: switch
    end subroutine hybridChargeParticle

    subroutine locateParticle(zpart,wzpart,wppart,wtpart0,wtpart1,&
        jtpart0,jtpart1,qpart,apart,ngyro,mp)
      use precision
      implicit none

      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart
      real(wp),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(wp),dimension(:,:) :: zpart
      integer,dimension(:,:) :: jtpart0,jtpart1
      integer ngyro,mp
    end subroutine locateParticle

    subroutine gkLoad(species_name,w_initial)
      use precision
      implicit none
      character(*),intent(in) :: species_name
      real(wp) w_initial
    end subroutine

    subroutine gkLoadParticle(meshn,tppp,nppp,zpart,zpart0,wppart,&
        wtpart0,wtpart1,density,flow,data1d,jtpart0,jtpart1,wzpart,&
        marker,markert,pmark,pflux,rdtem,zonal,zonalc,mesht,qpart,apart,&
        pload,nparam,ngyro,mp,mpmax,w_initial,ncyclep,zpart1,&
        trap,trapfracn,trapfrac,hybrid,pressurepara,pressureperp,phit,dnt,&
        phisave,dnsave)

      use precision
      implicit none
      integer pload,nparam,ngyro,mp,mpmax
      integer,optional :: ncyclep,trap,hybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart,w_initial
      real(wp),optional :: trapfracn
      real(wp),dimension(0:) :: meshn,mesht,pmark,pflux,rdtem,zonal,zonalc
      real(wp),dimension(:) :: wzpart,marker,markert
      real(wp),dimension(:),allocatable,optional :: trapfrac
      real(wp),dimension(:,:) :: tppp,nppp
      real(wp),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: density,flow,data1d
      real(wp),dimension(:,:),optional :: zpart1
      real(wp),dimension(0:,:),optional :: pressurepara,&
        pressureperp,phit,dnt
      real(wp),dimension(0:,:,:),optional :: phisave,dnsave
    end subroutine gkLoadParticle

    subroutine push(species_name,icycle,irksub,ihybrid)
      implicit none
      character(*),intent(in) :: species_name
      integer,intent(in),optional :: icycle,irksub,ihybrid
    end subroutine push

    subroutine gkPushParticle(zpart,zpart0,wppart,wtpart0,wtpart1,&
        data1d,jtpart0,jtpart1,diagpart,wzpart,rdtem,kapan,kapat,meshn,mesht,&
        pflux,qpart,apart,pload,ngyro,mpmax,mp,nparam,ncyclep,icycle,irksub,&
        mp1,zpart1,ihybrid,phit,dnt)
      use precision
      use particle_array,only: mpdiag
      implicit none

      !declaration of the dummy arguments
      integer pload,ngyro,mpmax,mp,nparam
      integer,optional :: ncyclep,icycle,irksub,ihybrid,mp1
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart
      real(wp),dimension(0:) :: kapan,kapat,meshn,mesht,pflux,rdtem
      real(wp),dimension(mpdiag) :: diagpart
      real(wp),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
      real(wp),dimension(0:,:) :: data1d
      real(wp),dimension(:,:),optional :: zpart1
      real(wp),dimension(0:,:),optional :: phit,dnt
    end subroutine gkPushParticle

    subroutine hybridPushParticle(zpart,zpart0,wppart,wtpart0,wtpart1,&
        data1d,jtpart0,jtpart1,diagpart,wzpart,rdtem,kapan,kapat,meshn,mesht,&
        pflux,qpart,apart,pload,ngyro,mpmax,mp,nparam,ncyclep,icycle,irksub,&
        mp1,zpart1,ihybrid,phit,dnt)
      use precision
      use particle_array,only: mpdiag
      implicit none

      integer pload,ngyro,mpmax,mp,nparam
      integer,optional :: ncyclep,icycle,irksub,ihybrid,mp1
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(wp) qpart,apart
      real(wp),dimension(:) :: wzpart
      real(wp),dimension(0:) :: kapan,kapat,meshn,mesht,pflux,rdtem
      real(wp),dimension(mpdiag) :: diagpart
      real(wp),dimension(:,:) :: zpart0,wppart,wtpart0,wtpart1
      real(wp),dimension(:,:) :: zpart
      real(wp),dimension(0:,:) :: data1d
      real(wp),dimension(:,:),optional :: zpart1
      real(wp),dimension(0:,:),optional :: phit,dnt
    end subroutine hybridPushParticle

    subroutine laplacian(scalar,scalar_out,none0bound)
      use precision
      use global_parameters,only:mgrid
      implicit none

      real(wp),dimension(mgrid),intent(in) :: scalar
      real(wp),dimension(mgrid),intent(out) :: scalar_out
      integer, optional, intent (in) :: none0bound
    end subroutine laplacian

    subroutine shiftParticle(zpart,zpart0,mpmax,mp,nparam,mtoroidal,&
        mype,left_pe,right_pe,myrank_toroidal,toroidal_comm0,numberpe,pi,&
        zeta0,zeta1)
      use precision
      implicit none
      
      integer mpmax,mp,nparam,mtoroidal,mype,left_pe,right_pe,myrank_toroidal,&
        toroidal_comm0,numberpe
      real(wp),dimension(:,:) :: zpart,zpart0
      real pi,zeta0,zeta1
    end subroutine shiftParticle

    subroutine axisExtrapolate(farray)
      use precision
      use global_parameters,only: mgrid

      real(wp),dimension(mgrid) :: farray
    end subroutine axisExtrapolate


#ifdef _OPENACC
    subroutine shift_cuda(zpart,zpart0,mpmax,mp,nparam,mtoroidal,mype,&
        left_pe,right_pe,myrank_toroidal,toroidal_comm0,numberpe,pi,zeta0,&
        zeta1) bind(C, name='shift_cuda')
      use iso_c_binding
      real(4),device :: zpart(nparam,*),zpart0(nparam,*)
      integer :: mp
      integer,value :: mpmax,nparam,mtoroidal,mype,left_pe,right_pe,&
        myrank_toroidal,toroidal_comm0,numberpe
      real(4),value :: pi,zeta0,zeta1
    end subroutine shift_cuda
#endif

    subroutine setupgpu(verbose) bind(C, name='setupgpu')
      integer,value :: verbose
    end subroutine setupgpu

    subroutine collision_fokkerplanck(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                                      ncyclep)
      use precision
      use global_parameters
      use particle_array
      use field_array
      use equilibrium
      implicit none

      !declaration of the dummy arguments
      integer mp,pcoll
      real(wp) qpart,apart,taupp
      real(wp),dimension(0:) :: mesht,meshn
      real(wp),dimension(:,:) :: tppp
      real(wp),dimension(:,:) :: zpart
      integer,optional :: ncyclep
    end subroutine collision_fokkerplanck

    subroutine collision_pitch(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                               ncyclep)
      use precision
      use global_parameters
      use particle_array
      use field_array
      use equilibrium
      implicit none

      !declaration of the dummy arguments
      integer mp,pcoll
      real(wp) qpart,apart,taupp
      real(wp),dimension(0:) :: mesht,meshn
      real(wp),dimension(:,:) :: tppp
      real(wp),dimension(:,:) :: zpart
      integer,optional :: ncyclep
    end subroutine collision_pitch

  end interface
end module interfaces
