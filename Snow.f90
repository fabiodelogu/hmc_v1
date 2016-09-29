!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!--------------------------------------------------------------
!		Aggiorna la maschera e la copertura nevosa
!--------------------------------------------------------------

SUBROUTINE Snow(step,idim,jdim,iHours,Trif,sMc,sMcGlac,dem,T,P,R,L,a2dAlbedo,a2dCdN,iSwitchNetRad,iSwitchLW,SM,SWE)
         
IMPLICIT NONE

!INCLUDE 'declaration_sram.f90'
INTEGER i,j,idim,jdim,step,iHours,iSwitchNetRad,iSwitchLW

REAL*8 Trif,sMc(idim,jdim),sMcGlac(idim,jdim)
REAL*8 dem(idim,jdim),SWE(idim,jdim),T(idim,jdim),P(idim,jdim),SM(idim,jdim),R(idim,jdim),L(idim,jdim),f(idim,jdim)
REAL*8 a2dAlbedo(idim,jdim),a2dCdN(idim,jdim)

INTEGER NStep
REAL*8 ro_w
REAL*8 sigma,lambda_f

write(222,*)'Inizio sub. SNOW '


!--------------------------------------------------------------
! DEFINIZIONE FATTORE DI NUVOLOSITà


WHERE(dem.ge.0.0.and.P.ge.0.2)f=1.1
WHERE(dem.ge.0.0.and.P.ge.1.0)f=1.2
WHERE(dem.ge.0.0.and.P.le.0.2)f=0.9
	
write(222,*)'Determinato Fattore Nuvolosita '								
!--------------------------------------------------------------
! C'è NEVE E T<Trif

WHERE(dem.ge.0.0.and.T.lt.Trif.and.P.gt.0.0)
	SWE=SWE+P	
	P=0.0
ENDWHERE

write(222,*)'Aggiunto al manto nevoso il contributo di precipitazione '
!--------------------------------------------------------------
! FUSIONE DELLA NEVE se T>=1.0
  CALL Snow_melting_ibrido(step,idim,jdim,sMc,sMcGlac,Trif,f,dem,T,R,L,SWE,a2dAlbedo,a2dCdN,iSwitchNetRad,iSwitchLW,SM) !con fattore nuvolosita                    
!--------------------------------------------------------------

! Fuori dal GHIACCIAIO
WHERE(dem.ge.0.0.AND.SWE.gt.0.0.AND.SWE.le.SM.AND.a2dCdN.NE.51) !si è fusa tutta la neve
	SM=SWE
	SWE=0.0
	P=P+SM
ELSEWHERE(dem.ge.0.0.AND.SWE.gt.0.0.AND.a2dCdN.NE.51)
	SWE=SWE-SM
	P=P+SM
ENDWHERE

! Sul GHIACCIAIO si può sempre fondere
WHERE(dem.ge.0.0.AND.a2dCdN.EQ.51) 
	SWE=SWE-SM
	P=P+SM
ENDWHERE

WHERE(dem.ge.0.0.AND.SWE.LT.0.0) !
	SWE=0.0
ENDWHERE


!--------------------------------------------------------------
! PRECIPITAZIONE E SNOWMELTING CUMULATI
!SM_b = SUM(SUM(SM,DIM=1,MASK=dem.ge.0.0)) !DIM=1 columns
!RA   = SUM(SUM(P,DIM=1,MASK=dem.ge.0.0))         !DIM=1 columns 
!--------------------------------------------------------------------


RETURN
END SUBROUTINE

