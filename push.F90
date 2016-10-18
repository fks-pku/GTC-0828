subroutine push(species_name,icycle,irksub,ihybrid)
  use global_parameters,only: nparami,nparame,nparamf,nparamfe,ncyclee,ncyclefe
  use particle_array
  use interfaces,only: gkPushParticle,hybridPushParticle
  use global_parameters,only: mype
  implicit none

  integer ierror
  character(*),intent(in) :: species_name
  integer,intent(in),optional :: icycle,irksub,ihybrid

  if(species_name=="thermal-ion")then
    call gkPushParticle(zion,zion0,wpion,wtion0,wtion1,data1di,&
      jtion0,jtion1,diagion,wzion,rdtemi,kapani,kapati,meshni,meshti,pfluxi,&
      qion,aion,iload,ngyroi,mimax,mi,nparami)
  elseif(species_name=="thermal-electron")then
    call hybridPushParticle(zelectron,zelectron0,wpelectron,wtelectron0,&
      wtelectron1,data1de,jtelectron0,jtelectron1,diagelectron,wzelectron,&
      rdteme,kapane,kapate,meshne,meshte,pfluxe,qelectron,aelectron,eload,&
      ngyroe,memax,me,nparame,ncyclee,icycle,irksub,me1,zelectron1,ihybrid,&
      phit,dnet)
  elseif(species_name=="fast-ion")then
    call gkPushParticle(zfast,zfast0,wpfast,wtfast0,wtfast1,data1df,&
      jtfast0,jtfast1,diagfast,wzfast,rdtemf,kapanf,kapatf,meshnf,meshtf,&
      pfluxf,qfast,afast,fload,ngyrof,mfmax,mf,nparamf)
#ifdef _FRC
  elseif(species_name=="fast-electron")then
    call gkPushParticle(zfaste,zfaste0,wpfaste,wtfaste0,wtfaste1,&
     data1dfe,jtfaste0,jtfaste1,diagfaste,wzfaste,rdtemfe,kapanfe,kapatfe,&
     meshnfe,meshtfe,pfluxfe,qfaste,afaste,feload,ngyrofe,mfemax,mfe,nparamfe)
#else
  elseif(species_name=="fast-electron")then
    call gkPushParticle(zfaste,zfaste0,wpfaste,wtfaste0,wtfaste1,&
     data1dfe,jtfaste0,jtfaste1,diagfaste,wzfaste,rdtemfe,kapanfe,kapatfe,&
     meshnfe,meshtfe,pfluxfe,qfaste,afaste,feload,ngyrofe,mfemax,mfe,&
     nparamfe,ncyclefe,icycle,irksub,mfe1,zfaste1)
#endif
  else
    write(*,*)'push.F90: wrong choice: ',species_name
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine push
