!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,M,NSTACK,indx(n)
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK),indxt,itemp
      REAL a,a2dTemp
	  do j=1,n
		indx(j)=j
	  enddo
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
		  indxt=indx(j)
          a=arr(indxt)
          do 11 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
			  indx(i+1)=indx(i)
            !arr(i+1)=arr(i)
11        continue
          i=l-1
2         indx(i+1)=indxt
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
		itemp=indx(k)
		indx(k)=indx(l+1)
		indx(l+1)=itemp
       
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
		  indx(l)=indx(ir)
		  indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
		  indx(l+1)=indx(ir)
		  indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
		  indx(l)=indx(l+1)
		  indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
		a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
		itemp=indx(i)
		indx(i)=indx(j)
		indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
!  (C) Copr. 1986-92 Numerical Recipes Software "a2dW#(1&.
