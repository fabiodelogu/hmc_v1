!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

!***********************************************************************************   
!	Subroutine lin_regr
!	Linear regression y = a + bx (flag = 1) , y = bx (flag = 0)
!
!	N = number of sample points
!	x(N),y(N) sample data
!	a,b = intersection, slope
!	r = regression coefficient
!***********************************************************************************

	subroutine LinRegr(N,x,y,a,b,r)
	Implicit none

	real*8 dNull
	PARAMETER(dNull=-9999)
    INTEGER N
	real*8 x(N),y(N)
	real*8 a, b, r 
	INTEGER flag 
	INTEGER i,j
	real*8 sx,sy,sxx,sxy,syy
	
	flag=1
	
	sx = 0.
	sy = 0.
	sxx = 0.
	sxy = 0.
	syy = 0.

	DO i = 1, N
		
		IF(y(i).ne.dNull)then
			sx = sx + x(i)
			sy = sy + y(i)
			sxx = sxx + x(i)*x(i)
			sxy = sxy + x(i)*y(i)
			syy = syy + y(i)*y(i)
		
		ELSE
	   
			sx = sx
			sy = sy 
			sxx = sxx 
			sxy = sxy 
			syy = syy
			N=N-1
		
		ENDIF
	END DO

	sx = sx / N
	sy = sy / N
	sxx = sxx / N
	sxy = sxy / N
	syy = syy / N

    b = (sxy-sx*sy)/(sxx-sx*sx)
    a = sy - b*sx
	r = (sxy-sx*sy)/sqrt(((sxx-sx*sx)*(syy-sy*sy)))

    IF (flag .eq. 0) Then
		b = sxy/sxx
		a = 0.
	END IF

	Return 
    
	END

!***********************************************************************************

