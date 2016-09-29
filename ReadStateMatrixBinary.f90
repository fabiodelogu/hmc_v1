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
 
subroutine ReadStateMatrixBinary(iRows,iCols,xstring)

implicit none
include 'DeclarationH.f90'

integer iRows,iCols
integer i,j
character*50 file
INTEGER*4 iLStr,iPathLength,iVarLength,iSize
CHARACTER*12  xstring
CHARACTER*500 sVariable
CHARACTER*500 sPathState,sFileState,sNameOut

iPathLength = iLStr(sPathResults)
700 FORMAT (1000(f9.2,1x))

sPathState=sPathResults(1:iPathLength)//sBasin(1:iNameLenght)

!-----------------------------------------------------------------------------------
!   Reading matrix V
!-----------------------------------------------------------------------------------

sVariable='V'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dV,10000,iLinux)
!-----------------------------------------------------------------------------------
!   Reading matrix Retention
!-----------------------------------------------------------------------------------
sVariable='Ret'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dRetention,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix Hydro Level
!-----------------------------------------------------------------------------------
sVariable='Wl'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dHydro,100000,iLinux)
!-----------------------------------------------------------------------------------
!   Reading matrix Routing
!-----------------------------------------------------------------------------------
sVariable='Rou'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dRouting,100000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix Water Table Level
!-----------------------------------------------------------------------------------
sVariable='Vw'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dVwt,10000,iLinux)


!-----------------------------------------------------------------------------------
!   Reading matrix 3D Tmarked
!-----------------------------------------------------------------------------------
sVariable='Tmk'
iSize=size(a3dTMarked, DIM=3)
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL Read3DMapBinary(iRows,iCols,iSize,sFileState,sNameOut,a3dTMarked,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix 3D Temp24
!-----------------------------------------------------------------------------------
sVariable='T24'
iSize=size(a3dTemp24, DIM=3)
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL Read3DMapBinary(iRows,iCols,iSize,sFileState,sNameOut,a3dTemp24,10000,iLinux)


!-----------------------------------------------------------------------------------
!   Reading matrix Ts
!-----------------------------------------------------------------------------------
sVariable='Ts'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dTs,10000,iLinux)


!-----------------------------------------------------------------------------------
!   Reading matrix DeepFloweExfiltration
!-----------------------------------------------------------------------------------
sVariable='DFE'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinary(iRows,iCols,sFileState,sNameOut,a2dDeepFlow,10000,iLinux)



!-----------------------------------------------------------------------------------
!   Reading matrix Snow Water Equivalent
!-----------------------------------------------------------------------------------
sVariable='SWE'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinarySkip(iRows,iCols,sFileState,sNameOut,a2dSWE,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix Density
!-----------------------------------------------------------------------------------
sVariable='Density'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinarySkip(iRows,iCols,sFileState,sNameOut,a2dRoS,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix Albedo
!-----------------------------------------------------------------------------------
sVariable='Albedo'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinarySkip(iRows,iCols,sFileState,sNameOut,a2dAlbedo,10000,iLinux)

!-----------------------------------------------------------------------------------
!   Reading matrix Age
!-----------------------------------------------------------------------------------
sVariable='Age'
iVarLength = iLStr(sVariable)
iPathLength = iLStr(sPathState)
sNameOut=sBasin(1:iNameLenght)//sVariable(1:iVarLength)//'_'//xstring
sFileState=sPathState(1:iPathLength)//sVariable(1:iVarLength)//'_'//xstring
CALL ReadMapBinarySkip(iRows,iCols,sFileState,sNameOut,a2dAge,10000,iLinux)

return
end subroutine