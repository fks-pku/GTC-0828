!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                            !
!                 Gyrokinetic Toroidal Code (GTC)                            !
!                          Version 3, 2009                                   !
!                             GTC Team                                       !
!                                                                            !
! Developers/users:                                                          !
! Z. Lin (zhihongl@uci.edu), I. Holod, Y. Xiao, W.L. Zhang, S. Klasky,       !
! S. Ethier, W.J. Deng, H.S. Zhang, Z.X. Wang, D. Fulton, J. McClenaghan,    !
! A. Kuley, P. Jiang, J. Bao, G. Dong, C.K. Lau, D.J. Liu, D. Spong,         !
! S. Taimourzadeh, X. Cheng, X. S. Wei, Y. Q. Liu, P. Wang, W. Joubert       !
!                                                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program GTC
  use global_parameters
  use particle_array
  use field_array
  use interfaces,only:push
  use magnetic_island,only:gradalphaIs,alphaIs
  use equilibrium
  implicit none

  integer i
  real(doubleprec) timecpu(10)
  character(len=16) cdate(4)

! Initialize MPI, OpenMP, Adios, PETSc, Random number generator etc
  CALL INITIAL(cdate,timecpu)
! input parameters, equilibrium & profile data, coefficients for gyroaveraging & interpolation
  CALL SETUP
! initialize particle position and velocity for all species
  call load
! calculate ion gather-scatter coefficients for chargei and pushi
  call locate("thermal-ion")
  if(fload>0)call locate("fast-ion")
  if(feload>0)call locate("fast-electron")
! write initialization snapshot 
  call snapshot
if(mstepall+istep==0)then
open(123456,file='phinspec.out',status='replace')
endif
  call timer(timecpu(10),timecpu(7))
! main time loop
  do istep=1,mstep
     do irk=1,2
! idiag=0: do time history diagnosis at irk=1 & istep=ndiag*N
        idiag=mod(irk+1,2)+mod(istep,ndiag)

#ifdef _BIGDATA
        if(idiag==0)then
           if(ndata3d==1)call DATAOUT3D !write 3D fluid data
           if(track_particles==1)then
             CALL PARTOUT !write particle data
             call timer(timecpu(10),timecpu(8))
           endif
        endif
#endif

! gradients of phi, phiind, fluidue, apara; and zonal flow
        CALL FIELD_GRADIENT
        call timer(timecpu(10),timecpu(6))
! push ion
#ifndef GPU_UM
        !$acc update device(gradphi,gradapara,gradphieff,sapara)
        !$acc update device(gradalphaIs,alphaIs)
#endif
        if(iload==9)then
          call pushifk_boris
        else
          call push("thermal-ion")
        endif
        call timer(timecpu(10),timecpu(1))
! redistribute ion across PEs
        call shift("thermal-ion")
        call timer(timecpu(10),timecpu(2))
! ion perturbed density
        call locate("thermal-ion")
        call charge("thermal-ion")
        call timer(timecpu(10),timecpu(3))
! push fast ion
        if(fload>0)then
           call push("fast-ion")
           call shift("fast-ion")
           call locate("fast-ion")
           call charge("fast-ion")
           call timer(timecpu(10),timecpu(5))
        endif

! push fast electron
        if(feload>0)then
#ifdef _FRC
          call push("fast-electron")
          call shift("fast-electron")
          call locate("fast-electron")
          call charge("fast-electron")
          if(irk==2.and.ecoll>0)then
             call lorentz('fast-electron')
             !call fokkerplanck('fast-electron')
          endif
#else
          do i=1,ncyclefe*irk
             call push("fast-electron",i,1)
             call shift("fast-electron")
             call push("fast-electron",i,2)
             call shift("fast-electron")
          enddo
          call locate("fast-electron")
          call charge("fast-electron")
#endif
        call timer(timecpu(10),timecpu(5))
        endif

        if(magnetic==0)then
! solve GK Poisson equation for phi using adiabatic electron
#ifdef _FRC
           CALL POISSON_FRC("ES-DKE")
#else
           CALL POISSON_SOLVER("ES-DKE")
#endif
        else
! push fields of apara, fluidne
           CALL PUSHFIELD
! solver fields of phi, phiind, fluidue
           CALL POISSON_SOLVER("EM-hybrid")
           CALL FIELD_SOLVER
        endif
        call timer(timecpu(10),timecpu(6))

        do ihybrid=1,nhybrid
! time derivative of effective potential, dphieff needs to be re-written for EM with KE
           CALL DPHIEFF
           call timer(timecpu(10),timecpu(6))           
#ifndef GPU_UM
           !$acc update device(dnet,gradpsi,phit,phip00)
#endif
           !push electron, sub-cycling
           do i=1,ncyclee*irk
! 1st RK step
              call push("thermal-electron",icycle=i,irksub=1,ihybrid=ihybrid)
              call shift("thermal-electron")
! 2nd RK step
              call push("thermal-electron",icycle=i,irksub=2,ihybrid=ihybrid)
              call shift("thermal-electron")

! electron-ion and electron-electron collision
              if(ecoll>0 .and. mod(i,ecoll)==0)then
                 call lorentz('thermal-electron')
     !            call fokkerplanck('thermal-electron')
              endif
           enddo

! nonadiabatic electron charge density
           call locate("thermal-electron")
           call charge("thermal-electron")
           call timer(timecpu(10),timecpu(4))
! solve GK Poisson equation using non-adiabatic electron
           if (magnetic==0) then
                CALL POISSON_SOLVER("ES-hybrid")
               ! CALL POISSON_SOLVER("ES-hybrid-old")
           endif
           call timer(timecpu(10),timecpu(6))
        enddo

!write diagnostics and 1D data
        if(idiag==0)CALL DIAGNOSIS
        call timer(timecpu(10),timecpu(8))
     enddo

! i-i collisions
     if(icoll>0 .and. mod(istep,icoll)==0)then
       call fokkerplanck('thermal-ion')
       call timer(timecpu(10),timecpu(1))
     endif
! profile snapshots, write particle information to restart file
     if(msnap/=0)then
       if(mod(istep,ndiag*2)==0 .or. istep*(1-irun)==ndiag)CALL MYSNAP
        if(mod(istep,mstep/msnap)==0 .or. istep*(1-irun)==ndiag)CALL SNAPSHOT
      !if(mod(istep,mstep/msnap)==0 .or. istep+irun==1)CALL PHISHOT
       if(mod(istep,mstep/msnap)==0)CALL RESTART_IO("write")
       call timer(timecpu(10),timecpu(8))
     endif
  enddo

! Program finalize
 CALL FINAL(cdate,timecpu)

end program GTC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
