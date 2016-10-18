subroutine mysnap
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use magnetic_island
  implicit none

integer i,j,isp,ij,nsnap
real(wp)pdum,dpx,dp2,eqmeshni(0:mpsi),eqmeshti(0:mpsi)
real,external::sprgpsi
character(len=64)cdum0,cdum2
!!$omp parallel do private(i,pdum,isp,dpx,dp2)
  do i=0,mpsi
     pdum=psimesh(i)
     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     dp2=dpx*dpx
     eqmeshti(i)=tipp(1,isp)+tipp(2,isp)*dpx+tipp(3,isp)*dp2
     eqmeshni(i)=nipp(1,isp)+nipp(2,isp)*dpx+nipp(3,isp)*dp2
 enddo


! open snapshot output file
  if(mype==0)then
     nsnap=mstepall+istep
     write(cdum0,'(i7.7,".out")')nsnap
  if(island==1)then
     cdum2='./densityn_history/'//'densityn'//trim(cdum0)
     open(6570,file=cdum2,status='replace')
     if(nsnap==8000)then
	open(6571,file='snap8000.out',status='replace')
	do i=0,mpsi
		do j=0,mtheta(i)
		ij=igrid(i)+j
if(nhybrid>0)then
write(6571,'(8f12.6)')sprgpsi(psimesh(i)),3.1415926*2*j/mtheta(i),eqmeshni(i),(densitye(1,ij)+1)*eqmeshni(i),eqmeshni(i)*eqmeshti(i)/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+eqmeshni(i)*(pressureepara(1,ij)))/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+eqmeshni(i)*(pressureeperp(1,ij)))/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+eqmeshni(i)*(pressureepara(1,ij)/3+2*pressureeperp(1,ij)/3))/eqmeshni(0)/eqmeshti(0)
else
		write(6571,'(8f12.6)')sprgpsi(psimesh(i)),3.1415926*2*j/mtheta(i),eqmeshni(i),(densityi(1,ij)+1)*eqmeshni(i),eqmeshni(i)*eqmeshti(i)/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+eqmeshni(i)*(pressureipara(1,ij)))/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+eqmeshni(i)*(pressureiperp(1,ij)))/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+eqmeshni(i)*(pressureipara(1,ij)/3+2*pressureiperp(1,ij)/3))/eqmeshni(0)/eqmeshti(0)     
endif
	enddo
	enddo
	close(6571)
	endif
do i=0, mpsi
        ij=igrid(i)
if(nhybrid>0)then
	write(6570,'(11f12.6)')sprgpsi(psimesh(i)),eqmeshni(i),(densityi(1,ij+mtheta(i)/2)+1)*eqmeshni(i),eqmeshni(i)*(densityi(1,ij)+1),eqmeshni(i)*(1+densitye(1,ij+mtheta(i)/2)),eqmeshni(i)*(1+densitye(1,ij)),eqmeshni(i)*eqmeshti(i)/eqmeshni(0)/eqmeshti(0),&
(pressureepara(1,ij)*eqmeshni(i)+eqmeshni(i)*eqmeshti(i))/eqmeshni(0)/eqmeshti(0),(pressureepara(1,ij+mtheta(i)/2)*eqmeshni(i)+eqmeshni(i)*eqmeshti(i))/eqmeshni(0)/eqmeshti(0),&
((2*pressureeperp(1,ij)+pressureepara(1,ij))*eqmeshni(i)/3+eqmeshni(i)*eqmeshti(i))/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+(2*pressureeperp(1,ij+mtheta(i)/2)+pressureepara(1,ij+mtheta(i)/2))*eqmeshni(i)/3)/eqmeshni(0)/eqmeshti(0)
else
	write(6570,'(9f12.8)')sprgpsi(psimesh(i)),eqmeshni(i),(densityi(1,ij+mtheta(i)/2)+1.0)*eqmeshni(i),eqmeshni(i)*(densityi(1,ij)+1.0),&
eqmeshni(i)*eqmeshti(i)/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+pressureipara(1,ij)*eqmeshni(i))/eqmeshni(0)/eqmeshti(0),&
(eqmeshni(i)*eqmeshti(i)+pressureipara(1,ij+mtheta(i)/2)*eqmeshni(i))/eqmeshni(0)/eqmeshti(0),(eqmeshni(i)*eqmeshti(i)+(pressureipara(1,ij)+2*pressureiperp(1,ij))/3*eqmeshni(i))/eqmeshni(0)/eqmeshti(0),&
(eqmeshni(i)*eqmeshti(i)+(pressureipara(1,ij+mtheta(i)/2)+2*pressureiperp(1,ij+mtheta(i)/2))/3*eqmeshni(i))/eqmeshni(0)/eqmeshti(0)
endif      
  enddo
     close(6570)
     endif
  endif
end subroutine mysnap
