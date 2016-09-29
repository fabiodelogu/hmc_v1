!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************
!	Subroutine read_info_data
!	Legge il file con le info delle prese-rilasci: sBasinDamPreseInfo.txt
!***********************************************************************************  
 
subroutine ReadInfoPrese(iRows,iCols)

implicit none
include 'DeclarationH.f90'

INTEGER iRows,iCols
INTEGER i,iFlagInterp
CHARACTER*50 info,text
CHARACTER*12 xstring
INTEGER*4 iLStr,iPathLenght,inC,inCt,ii
REAL*8 dDintegr !Dt integrazione Hcdrift
REAL*8 dQr
!Opening info dighe
OPEN(20,file=sFileInfoPrese,status='old')
  
 
READ(20,*)
READ(20,*)
inCt=1
DO i=1,iNril
		READ(20,*)
		!---nome sezioni di rilascio
	    READ(20,*)a1sNameRilasci(i)
		!----coord i e j dell'opera di rilascio
		READ(20,*)a2dXYRilascio(i,2),a2dXYRilascio(i,1)
		a2dXYRilascio(i,2)=iRows-a2dXYRilascio(i,2)+1

		!---numero sezioni di presa
	    READ(20,*)inC
		
	 	a1dQmaxRil(i)=0.0	
		DO ii=inCt,inC+inCt-1
		 !----nome della centrale
	     READ(20,*)
		 a1sNamePrese(ii)=a1sNameRilasci(i)
		 !----tempo di corrivazione presa rilascio
	     READ(20,*)  a1dTcorrprese(ii)

		 !----coord i e j della sezione di presa
		 READ(20,*)a2dXYPresa(ii,2),a2dXYPresa(ii,1)
         a2dXYPresa(ii,2)=iRows-a2dXYPresa(ii,2)+1
		 !----Q max presa
	     READ(20,*)dQr
		 if(dQr.gt.0.0)then
			a1dQmaxRil(i)=a1dQmaxRil(i)+dQr
		 endif
		 !----Q min eco
	     READ(20,*)dQr
		 if(dQr.gt.0.0)then
			a1dQminEco(ii)=dQr
		 endif

		 !----Peso della presa rispetto al rilascio
	     READ(20,*)a1dPesoPresa(ii)
		 
		ENDDO
		inCt=inCt+inC
		
	    
ENDDO
CLOSE(20)



return
end subroutine