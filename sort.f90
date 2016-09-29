!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,a2dK,l,istack(NSTACK)
      REAL a,a2dTemp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        a2dK=(l+ir)/2
        a2dTemp=arr(a2dK)
        arr(a2dK)=arr(l+1)
        arr(l+1)=a2dTemp
        if(arr(l+1).gt.arr(ir))then
          a2dTemp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=a2dTemp
        endif
        if(arr(l).gt.arr(ir))then
          a2dTemp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=a2dTemp
        endif
        if(arr(l+1).gt.arr(l))then
          a2dTemp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=a2dTemp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        a2dTemp=arr(i)
        arr(i)=arr(j)
        arr(j)=a2dTemp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
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
