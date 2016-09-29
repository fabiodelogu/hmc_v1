!***********************************************************************************
!	EUPL
!	Code Version 1.0
!	Authors: Silvestro Francesco, Gabellani Simone, Delogu Fabio
!	This code is distributed following the European Union Public Licence:
!	https://joinup.ec.europa.eu/software/page/eupl/licence-eupl
!***********************************************************************************

SUBROUTINE DeterminaData_stepParam(xstring,xstring_up,iDtDate)

IMPLICIT NONE


INTEGER anno,iBanno,mese,giorno,ora,flag,min,iDtDate,iDiff,i
CHARACTER*12  xstring,xstring_up

READ( xstring(1:4), '(i4)' )  anno
READ( xstring(5:6), '(i2)' )  mese
READ( xstring(7:8), '(i2)' )  giorno
READ( xstring(9:10), '(i2)' ) ora
READ( xstring(11:12), '(i2)') min

flag=0
iDiff=60-iDtDate !Quantità da togliere per il cambio di ora
iBanno=int(anno/4)*4 !Anno bisestile più prossimo

IF(giorno.EQ.30.AND.ora.EQ.23.AND.min.EQ.iDiff.AND.flag.EQ.0)THEN
    IF(mese.EQ.4.OR.mese.EQ.6.OR.mese.EQ.9.OR.mese.EQ.11)THEN
			mese=mese+1
			giorno=1
			ora=0
			min=0
			flag=1
	ENDIF
ENDIF

IF(giorno.EQ.31.AND.ora.EQ.23.AND.min.EQ.iDiff.AND.flag.EQ.0)THEN
	IF(mese.EQ.1.OR.mese.EQ.3.OR.mese.EQ.5.OR.mese.EQ.7.OR.mese.EQ.8.OR.mese.EQ.10)THEN
		mese=mese+1
		giorno=1
		ora=0
		min=0
		flag=1
	ENDIF
ENDIF

IF(giorno.EQ.31.AND.ora.EQ.23.AND.min.EQ.iDiff.AND.flag.EQ.0)THEN
	IF(mese.EQ.12)THEN
		anno=anno+1
		mese=1
		giorno=1
		ora=0
		min=0
		flag=1
	ENDIF
ENDIF

IF(mese.EQ.2.AND.flag.EQ.0)THEN
	IF(giorno.EQ.29.AND.ora.EQ.23.AND.min.EQ.iDiff.AND.flag.EQ.0)THEN
		IF(anno.EQ.iBanno)THEN
			mese=mese+1
			giorno=1
			ora=0
			min=0
			flag=1
		ENDIF
	ENDIF
ENDIF


IF(mese.EQ.2.AND.flag.EQ.0)THEN
	 IF(giorno.EQ.28.AND.ora.EQ.23.AND.min.EQ.iDiff.AND.flag.EQ.0)THEN
		
			IF(anno.EQ.iBanno)THEN
				WRITE(*,*)'BISESTILE ',anno
				mese=mese
				giorno=giorno+1
				ora=0
				min=0
				flag=1
			ELSE
				mese=mese+1
				giorno=1
				ora=0
				min=0
				flag=1
			ENDIF

	ENDIF
ENDIF


IF(flag.EQ.0.AND.ora.EQ.23.AND.min.EQ.iDiff)THEN
	giorno=giorno+1
	ora=0
	min=0
	flag=1
ENDIF


IF(flag.EQ.0.AND.min.NE.iDiff)THEN
	min=min+iDtDate
	ora=ora
	flag=1
ENDIF

IF(flag.EQ.0.AND.ora.NE.23.AND.min.EQ.iDiff)THEN
	ora=ora+1
	min=0
ENDIF


WRITE( xstring_up(1:4), '(i4)' )  anno

IF(mese.GT.9)THEN
	WRITE( xstring_up(5:6), '(i2)' )mese
	ELSE
	WRITE( xstring_up(5:5), '(i1)' )0
	WRITE( xstring_up(6:6), '(i1)' )mese
ENDIF
IF(giorno.GT.9)THEN
	WRITE( xstring_up(7:8), '(i2)' )  giorno
	ELSE
	WRITE( xstring_up(7:7), '(i1)' )  0
	WRITE( xstring_up(8:8), '(i1)' )  giorno
ENDIF
IF(ora.GT.9)THEN
	WRITE( xstring_up(9:10), '(i2)' ) ora
	ELSE
	WRITE( xstring_up(9:9), '(i1)' ) 0
	WRITE( xstring_up(10:10), '(i1)' ) ora
ENDIF


IF(min.GT.0)THEN
	WRITE( xstring_up(11:12), '(i2)' ) min
	ELSE
	WRITE( xstring_up(11:11), '(i1)' ) 0
	WRITE( xstring_up(12:12), '(i1)' ) min
ENDIF


RETURN
END SUBROUTINE
