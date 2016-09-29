
! DETERMINA SNOWMELTING CON IL METODO IBRIDO
!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

SUBROUTINE   Snow_melting_ibrido(step,idim,jdim,sMc,sMcGlac,Trif,f,dem,T,R,L,SWE,a2dAlbedo,a2dCdN,iSwitchNetRad,iSwitchLW,SM)                   

  IMPLICIT NONE

  INTEGER i,j,idim,jdim,step,NStep,iSwitchNetRad,iSwitchLW
  REAL*8 Trif,ro_w,sMc(idim,jdim),sMcGlac(idim,jdim)
  REAL*8 sigma,lambda_f
  REAL*8 f(idim,jdim),dem(idim,jdim),T(idim,jdim),R(idim,jdim),SWE(idim,jdim),SM(idim,jdim)
  REAL*8 a2dAlbedo(idim,jdim),L(idim,jdim),a2dCdN(idim,jdim)

  write(222,*)'Determinazione Melting '

  sigma=0.0000000049   	![MJ/(m^2 day-1 K^(-4))]
  ro_w=1000				! densità acqua [kg/m^3]

  NStep=60/step*24 !Step in un giorno
  IF(step.EQ.15)THEN
     NStep=(13*4)
  ENDIF
  IF(step.EQ.60)THEN
     NStep=24
  ENDIF
  IF(step.EQ.1440)THEN
     NStep=1
  ENDIF

  sigma=sigma/24 ![MJ/(m^2 h-1 K^(-4))]

  ![MJ/kg] calore latenete di fusione
  lambda_f=0.334

  !e1=-0.02+0.261*exp(-7.77*0.0001*T(i,j)**2)	        ! emissività netta tra atmosfera e suolo
  !L=(-sigma*f*e1*(T+273.2)**4)						! long wave radiation	
  !K=K*(1.0-albedo)*3600/1000000.0					! trasformazione in MJ/m^2/h e albedo 							
  !SM1=((1000.0*(K+L)/(ro_w*lambda_f)+sMc*(T)))/24.0	! SNOWMELTING

  !----Radiazione netta osservata
  IF(iSwitchNetRad==1)THEN
     write(222,*)'Uso Rad Netta Osservata '
     WHERE(dem.GE.0.0.AND.T.GE.Trif.AND.SWE.GT.0.0)
        R=R*step*60*10**(-6) ! trasformazione in MJ/m^2/h
        L=0
     ENDWHERE
  ENDIF

  !----Radiazione SW e LW Calcolate
  IF(iSwitchNetRad==0.AND.iSwitchLW==0)THEN
     write(222,*)'Calcolo LW e SW '
     WHERE(dem.GE.0.0.AND.T.GE.Trif.AND.SWE.GT.0.0)
        R=R*(1.0-a2dAlbedo)*step*60*10**(-6) ! trasformazione in MJ/m^2/h
        L=-sigma*f*(-0.02+0.261*dexp(-7.77*0.0001*T**2))*(T+273.2)**4
     ENDWHERE
  ENDIF
  !----Radiazione SW calcolata e LW osservata
  IF(iSwitchNetRad==0.AND.iSwitchLW==1)THEN
     write(222,*)'Calcolo LW e uso SW osservata '
     WHERE(dem.GE.0.0.AND.T.GE.Trif.AND.SWE.GT.0.0)
        R=R*(1.0-a2dAlbedo)*step*60*10**(-6)
        L=L*step*60*10**(-6) ! trasformazione in MJ/m^2/h 
     ENDWHERE
  ENDIF
  
  WHERE(dem.GE.0.0.AND.T.GE.Trif.AND.SWE.GT.0.0)
     SM=1000.0/(ro_w*lambda_f)*(R+L)+sMc*T
     SM=SM/float(NStep)
  ELSEWHERE
     SM=0.0
  ENDWHERE

  !Su ghiacciaio fondo sempre - albedo=0.9
  !WHERE(dem.GE.0.0.AND.T.GE.Trif.AND.SWE.EQ.0.0.AND.a2dCdN.EQ.51)
  !   R=R*(1.0-0.6)*step*60*10**(-6)
  !   L=-sigma*f*(-0.02+0.261*exp(-7.77*0.0001*T**2))*(T+273.2)**4
  !   SM=1000.0/(ro_w*lambda_f)*(R+L)+sMcGlac*T
  !   SM=SM/NStep 
  !ENDWHERE


  write(222,*)'Fine Determinazione Melting'

  ! emissività netta tra atmosfera e suolo: e1=-0.02+0.261*exp(-7.77*0.0001*T**2)		
  ! long wave radiation: L=(-sigma*f*(-0.02+0.261*exp(-7.77*0.0001*T**2))*(T+273.2)**4)		
  ! trasformazione in MJ/m^2/h e albedo=0.23

  WHERE(SM.lt.0.0)
     SM=0.0
  ENDWHERE

  RETURN
END SUBROUTINE Snow_melting_ibrido

