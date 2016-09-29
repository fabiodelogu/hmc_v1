
!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine read_land_data
!	Legge il file sBasin.info.txt
!***********************************************************************************  
 
subroutine WriteStateMatrixBinary(iRows,iCols,xstring)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file
INTEGER*4 iLStr,iPathLenght,iSize
CHARACTER*12  xstring
CHARACTER*500 sVariable
CHARACTER*500 sPathState
REAL*8 tmp(iRows,iCols)

iPathLenght = iLStr(sPathResults)
700 FORMAT (1000(f9.2,1x))

sPathState=sPathResults(1:iPathLenght)//sBasin(1:iNameLenght)

!-----------------------------------------------------------------------------------
!   Writing matrix V
!-----------------------------------------------------------------------------------

sVariable='V'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dV,10000,iLinux)
!-----------------------------------------------------------------------------------
!   Writing matrix Retention
!-----------------------------------------------------------------------------------
sVariable='Ret'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dRetention,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix Hydro Level
!-----------------------------------------------------------------------------------

sVariable='Wl'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dHydro,100000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix Routing
!-----------------------------------------------------------------------------------
sVariable='Rou'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dRouting,100000,iLinux)


!-----------------------------------------------------------------------------------
!   Writing matrix Water Table Level
!-----------------------------------------------------------------------------------
tmp=0
WHERE(a2dDem.gt.0.0)
	tmp=(a2dDem-a2dVwt)*1000
ENDWHERE
sVariable='Vw'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,tmp,10000,iLinux)


!-----------------------------------------------------------------------------------
!   Writing matrix 3D Tmarked
!-----------------------------------------------------------------------------------
sVariable='Tmk'
iSize=size(a3dTMarked, DIM=3)
CALL Write3DMapBinary(xstring,iRows,iCols,iSize,sVariable,sPathState,a3dTMarked,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix 3D Temp24
!-----------------------------------------------------------------------------------
sVariable='T24'
iSize=size(a3dTemp24, DIM=3)
CALL Write3DMapBinary(xstring,iRows,iCols,iSize,sVariable,sPathState,a3dTemp24,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix Soil Temperature
!-----------------------------------------------------------------------------------
sVariable='Ts'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dTs,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix DeepFloweExfiltration
!-----------------------------------------------------------------------------------
sVariable='DFE'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dDeepFlow,10000,iLinux)


!-----------------------------------------------------------------------------------
!   Writing matrix SWE
!-----------------------------------------------------------------------------------
sVariable='SWE'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dSWE,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix SWE
!-----------------------------------------------------------------------------------
sVariable='Density'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dRoS,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix SWE
!-----------------------------------------------------------------------------------
sVariable='Albedo'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dAlbedo,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix SWE
!-----------------------------------------------------------------------------------
sVariable='Age'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dAge,10000,iLinux)

return
end subroutine