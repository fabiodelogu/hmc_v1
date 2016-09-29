!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine Create_Meteo_GRID
!	geolocate rainfall matrix from original to basin grid
!***********************************************************************************  
 
subroutine Create_Meteo_GRID(idim,jdim,idimrain,jdimrain)

implicit none
include 'DeclarationH.f90'

integer idim,jdim,idimrain,jdimrain,m,sim_length
integer i,j,unit_grid,ios
character*50 info,text
character*400 rain_file  !file formato grid contenente i dati di precipitazione
real*8 rlat,rlon,nx,ny

REAL*8 a2dEx(idim,jdim)



a2iYLat=-9999
a2iXLon=-9999
!write(*,*)idim,jdim,idimrain,jdimrain
! georeferenziazione della matrice di pioggia
do i=1,idim
	do j=1,jdim
!FORALL(i=1:idim,j=1:jdim,a2dDem(i,j).ge.0.0)
	
		rlat=dXDemLat+(dDemPasLat*(i-1))
		rlon=dXDemLon+(dDemPasLon*(j-1))
		nx=nint((rlon-dXCNLon)/dCNPasLon)+1
		ny=nint((rlat-dXCNLat)/dCNPasLat)+1
		
		if (nx.gt.0.and.nx.le.jdimrain.and. &
				ny.gt.0.and.ny.le.idimrain) then
			a2iYLat(i,j)=ny
			a2iXLon(i,j)=nx
	    else
			write(*,*)nx,ny
			!Check fo approximation
			if(nx.lt.jdimrain+3)then
				nx=jdimrain
				a2iXLon(i,j)=nx
			endif
			if(ny.lt.idimrain+3)then
				ny=idimrain
				a2iYLat(i,j)=ny
			endif
		endif

!ENDFORALL
	enddo
enddo

!CALL Geoloc_meteo_file_GRID(idim,jdim,idimrain,jdimrain,a2dDem,a2dEx)
!open(unit=18,file='X.txt')
!do i=1,idim
!     write(18,'(2000(f8.1,1x))')(a2dXLon(i,j),j=1,jdim)
!end do
!close(18)
return
end subroutine