!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for saving data to the output file				*
!*																	*
!********************************************************************
!*																	*
!* The mesh is shown below											*
!*																	*
!* (bottom left corner)												*
!*																	*
!* v----o----v----o  y(2)											*
!*      |         |													*
!* P    u    P    u													*
!*      |         |													*
!* v----o----v----o  y(1)											*
!*	    |         |													*
!* P    u    P    u													*
!*																	*
!*     x(1)     x(2)												*
!********************************************************************
! Modified by Jessica Sanders 16 juillet 2008 to remove temp data
! And again 15 Aout 2008 to remove volume variable "F", and add third
! phase
SUBROUTINE SAVEDATA(Nx,Ny,Lx,Ly,time,x,y,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v,file_count,int_data)

	IMPLICIT NONE
	
	! local integers
	INTEGER :: size,ios,i,j,k,myoutput

	! integers
	INTEGER, INTENT(IN):: Nx,Ny,Lx,Ly,file_count 
	
	! real variables
	real(kind=8) :: deltaTavg
	REAL(kind=8), INTENT(IN) :: time

	! real arrays
	REAL(kind=8), DIMENSION(Nx+1), INTENT(IN) :: x, int_data
	REAL(kind=8), DIMENSION(Ny+1), INTENT(IN) :: y
	REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS,H,P
	REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: s_phiLS,s_H
	REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: w_phiLS,w_H
	REAL(kind=8), DIMENSION(Nx+1,Ny+2), INTENT(IN) :: u
	REAL(kind=8), DIMENSION(Nx+2,Ny+1), INTENT(IN) :: v
	
	! characters
	CHARACTER (LEN=20) :: filename,string1,string2
	CHARACTER (LEN=5) :: string
 
	
	! creating the filename and printing data to scheme
	CALL to_string(file_count, string, size)
	
	if (file_count < 10) then
		string1 = "file00"  ! using two zeros as a place setter, file00x
	elseif (file_count >= 10 .and. file_count < 100) then
		string1 = "file0"  ! using two zeros as a place setter, file0x
	else
		string1 = "file"
	endif
	string2 = string(:size)   !I am extracting the text number
	filename = TRIM(string1)//TRIM(string2)//".dat"
	
	
	! writing to the screen
	PRINT *, " "
	PRINT *, "-----"
	PRINT *, "saving: ", filename
	PRINT *, "time =", time, " (sec)"
	PRINT *, " "

	! CREATING DATA FILE
	OPEN (UNIT=1,FILE=filename,STATUS="REPLACE", IOSTAT=ios)
	
	! TECHPLOT OUTPUT FILE
	WRITE (UNIT=1,FMT=50) 'VARIABLES = "X", "Y", "PHILS", "S_PHILS", "H", "S_H", "P", "U", "V"'
	WRITE (UNIT=1,FMT=60) 'ZONE I=',Nx+1,'J=',Ny+1,'F=POINT'
	! note: the I ZONE corresponds to j indice and the J ZONE corresponds to the i indice
	DO j=1,Ny+1
		DO i=1,Nx+1
			WRITE (UNIT=1,FMT=70)	x(i), y(j), &
									0.25*(phiLS(i,j)+phiLS(i,j+1)+phiLS(i+1,j)+phiLS(i+1,j+1)), &
									0.25*(s_phiLS(i,j)+s_phiLS(i,j+1)+s_phiLS(i+1,j)+s_phiLS(i+1,j+1)), &
									0.25*(w_phiLS(i,j)+w_phiLS(i,j+1)+w_phiLS(i+1,j)+w_phiLS(i+1,j+1)), &
									0.25*(H(i,j)+H(i,j+1)+H(i+1,j)+H(i+1,j+1)), &
									0.25*(s_H(i,j)+s_H(i,j+1)+s_H(i+1,j)+s_H(i+1,j+1)), &
									0.25*(w_H(i,j)+w_H(i,j+1)+w_H(i+1,j)+w_H(i+1,j+1)), &
									0.25*(P(i,j)+P(i,j+1)+P(i+1,j)+P(i+1,j+1)), &
									0.5*(u(i,j)+u(i,j+1)), 0.5*(v(i,j)+v(i+1,j)), &
                                                                        int_data(i)
		ENDDO
	END DO

	CLOSE (UNIT=1)

	! writing the wall temperature gradient to file
	!!if (file_count == 0) then
	!!	OPEN (UNIT=2,FILE="Nusselt.dat",STATUS="REPLACE",IOSTAT=ios)
	!!	WRITE (UNIT=2, FMT=50) 'This file stores the average temperature gradient at the wall.'
	!!	WRITE (UNIT=2, FMT=50) 'The left column is Time and the right column is dT/dy.'
	!!else
	!!	OPEN (UNIT=2,FILE="Nusselt.dat",STATUS="OLD",POSITION="APPEND",IOSTAT=ios)
	!!endif

	!!deltaTavg = 0.0
	!!do i=1,Nx+2
	!!	deltaTavg = deltaTavg + (T(i,2)-T(i,1))/(y(2)-y(1))
	!!enddo
	!!deltaTavg = deltaTavg/(Nx+2)
	!!WRITE (UNIT=2, FMT=70) time, deltaTavg
		
	
	!!CLOSE (UNIT=2)
		

	!FORMATS
50  FORMAT (A)
60  FORMAT (A, I3, 1X, A, I3, 1X ,A)
70	FORMAT (9F15.4)  !9F10.4

80  FORMAT (A)
90	FORMAT (F10.4)
100 FORMAT (I8)

END SUBROUTINE SAVEDATA



!*************************************************************
!   Converts the integer stored in number into an ascii 
!   string.  The string is returned in string.  The number of 
!   digits is returned in size.
SUBROUTINE to_string(number, string, size)
	integer number
	character *(*) string
	integer size

	character*100 temp
	integer local
	integer last_digit
	integer i

	local = number
        i = 0

!   strip digits off starting with least significant
!   do-while loop
100	    last_digit = mod(local,10)
	    local = local/10
	    i = i + 1
	    temp(i:i) = char(last_digit + ichar('0'))
	if (local.ne.0) go to 100
	
	size = i

!  reverse digits
	do 200 i = 1, size
	    string(size-i+1:size-i+1) = temp(i:i)
200	continue
!
	RETURN
END SUBROUTINE to_string
!   to_string
!
!
