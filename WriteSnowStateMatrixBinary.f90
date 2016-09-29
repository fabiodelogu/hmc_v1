
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
 
subroutine WriteSnowStateMatrixBinary(iRows,iCols,xstring)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file
INTEGER*4 iLStr,iPathLenght,iSize
CHARACTER*12  xstring,sVariable
CHARACTER*500 sPathState
REAL*8 tmp(iRows,iCols)

iPathLenght = iLStr(sPathResults)
700 FORMAT (1000(f9.2,1x))

sPathState=sPathResults(1:iPathLenght)//sBasin(1:iNameLenght)

!-----------------------------------------------------------------------------------
!   Writing matrix SWE
!-----------------------------------------------------------------------------------

sVariable='SWE'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dSWE,10000,iLinux)
!-----------------------------------------------------------------------------------
!   Writing matrix Density
!-----------------------------------------------------------------------------------
sVariable='Density'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dRoS,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix Albedo
!-----------------------------------------------------------------------------------

sVariable='Albedo'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dAlbedo,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Writing matrix Routing
!-----------------------------------------------------------------------------------
sVariable='Age'
CALL WriteMeteoMapBinary(xstring,iRows,iCols,sVariable,sPathState,a2dAge,10000,iLinux)




return
end subroutine