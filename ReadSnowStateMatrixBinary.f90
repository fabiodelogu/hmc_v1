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
 
subroutine ReadSnowStateMatrixBinary(iRows,iCols,xstring)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file
INTEGER*4 iLStr,iPathLength,iVarLength,iSize
CHARACTER*12  xstring,sVariable
CHARACTER*500 sPathState,sFileState,sNameOut

iPathLength = iLStr(sPathResults)
700 FORMAT (1000(f9.2,1x))

sPathState=sPathResults(1:iPathLength)//sBasin(1:iNameLenght)

!-----------------------------------------------------------------------------------
!   Reading matrix SWE
!-----------------------------------------------------------------------------------

sVariable='SWE'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dSWE,10000,iLinux)
!-----------------------------------------------------------------------------------
!   Reading matrix Density
!-----------------------------------------------------------------------------------
sVariable='Density'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dRoS,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix Hydro Level
!-----------------------------------------------------------------------------------
sVariable='Albedo'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dAlbedo,10000,iLinux)
!-----------------------------------------------------------------------------------
!   Reading matrix Routing
!-----------------------------------------------------------------------------------
sVariable='Age'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dAge,10000,iLinux)



return
end subroutine