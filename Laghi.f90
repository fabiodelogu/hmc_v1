!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

subroutine Laghi(iRows,iCols)
include 'DeclarationH.f90'

integer iRows,iCols
integer i,ii,iii,j,jj,jjj,kkk,ttt
integer a,b,perc_tot

real*8 pi,DD,diff,fn,fMean,fNumPen
real*8 pdistance(iRows,iCols),LDD(iRows,iCols),mask_perc_tot(iRows,iCols)
real*8 diff_DD(iRows,iCols),pend(iRows,iCols),pend2(iRows,iCols),a2dTA(iRows,iCols)

!------------------------------------------------------------------


!-----------------------------------------------------------------
!     CALCOLA ZONE DEPRESSE E LE ALZA 
!     PUNTATORE MM(J,I)
!
!     - 0= VA BENE
!     - 1= INNALZATO
!-----------------------------------------------------------------

      i1=-1
      i2=1
      j1=-1
      j2=1
      ix=2
      jx=2
      kkk=0
 1000 kk=0
      do 11 i=ix,iCols-1
         do 10 j=jx,iRows-1
            if(a2dDem(j,i).lt.0) go to 10
            zmin=1.e20
            zx=a2dDem(j,i)
            
            do 2 iii=i1,i2
               do 2 jjj=j1,j2
                  if(iii.eq.0.and.jjj.eq.0) go to 22
                  ii=iii+i
                  jj=jjj+j
                  if(zx.gt.a2dDem(jj,ii)) go to 10
                  if(zmin.gt.a2dDem(jj,ii)) zmin=a2dDem(jj,ii)
 22               continue
 2          continue
            if(zx.le.zmin) then
               a2dDem(j,i)=zmin+.4
               kkk=kkk+1
               kk=1
               if (i.ne.2)ix=i-1
               if (j.ne.2)jx=j-1
               if(kkk/1000*1000.eq.kkk) then
                  write(*,'(i6,2i4,f8.2)')kkk,i,j,a2dDem(j,i)
               end if
               go to 1000
            end if
 10      continue
         jx=2
 11    continue
       
       write(*,*)'* ',kkk
       if(kk.eq.1) go to 1000
!      write(*,*)'End subroutine LAGO'



return
end subroutine