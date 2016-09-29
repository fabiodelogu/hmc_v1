!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine geoloc_rain_file_GRID
!	geolocate rainfall matrix from original to basin grid
!***********************************************************************************  
 
subroutine Geoloc_meteo_file_GRID(idim,jdim,idimrain,jdimrain,a2dMatrixO,a2dMatrixDem)

implicit none
include 'DeclarationH.f90'

integer idim,jdim,idimrain,jdimrain
integer i,j

real*8 a2dMatrixO(idimrain,jdimrain),a2dMatrixDem(idim,jdim)
	
! georeferenziazione della matrice meteo
FORALL (i=1:idim,j=1:jdim,a2dDem(i,j).ge.0.0)
	a2dMatrixDem(i,j)=a2dMatrixO(a2iYLat(i,j),a2iXLon(i,j))
ENDFORALL
return
end subroutine