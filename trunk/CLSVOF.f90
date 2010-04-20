!********************************************************************
!*																	*
!*							   CLSVOF								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a subroutine for tracking the volume using the CLSVOF	*
!* method.															*
!*																	*	
!********************************************************************
SUBROUTINE CLSVOF(Nx,Ny,u,v,uINT,vINT,phi,F,A,B,C,deltaX,deltaY,deltaT,&
				  xwhole,ywhole,xdirgo,ydirgo,tol,rho_one,rho_two)

	implicit none
	
	!----------------------------------- HEADER -----------------------------------
	! integers
	integer, intent(in) :: Nx,Ny,xdirgo,ydirgo

	! reals
	real(kind=8), intent(in) :: deltaX,deltaY,deltaT,tol,rho_one,rho_two

	! real arrays
	real(kind=8), dimension(Nx+1,Ny+2), intent(in) :: u,uINT
	real(kind=8), dimension(Nx+2,Ny+1), intent(in) :: v,vINT
	real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi,F,A,B,C
	real(kind=8), dimension(Nx+2), intent(in) :: xwhole
	real(kind=8), dimension(Ny+2), intent(in) :: ywhole

	! local integers
	integer :: i,j,iter,scheme,dummyint

	! local reals
	real(kind=8) :: Res,uchar,deltaArtT,mynorm,s_nn,s_ss,s_ee,s_ww

	! local real arrays
	real(kind=8), dimension(Nx+2,Ny+2) :: phi_tilda,phi_carrot,phi_dummy,&
										  F_tilda,H,Density
	
	! externals
	real(kind=8), external :: advectVOF,advectval,finaladvectval,LSHEAVY,NORMDIFFPHI
	!------------------------------------------------------------------------------
	

	!rho_one = 0.0
	!rho_two = 1.0


	! ------------- ADVECT DIRECTION 1 -------------
	! Advect phi in direction 1
    do i=2,Nx+1
        do j=2,Ny+1
			if (i==2 .or. i==Nx+1 .or. j==2 .or. j==Ny+1) then
				s_nn = 0.0
				s_ss = 0.0
				s_ee = 0.0
				s_ww = 0.0
				scheme = 1 ! do first-order donor cell
			else
				s_nn = phi(i,j+2)
				s_ss = phi(i,j-2)
				s_ee = phi(i+2,j)
				s_ww = phi(i-2,j)
				scheme = 1 ! do higher-order scheme		************ the higher order scheme assume div(V)=0
			endif
            phi_tilda(i,j) = advectval(phi(i,j),phi(i,j+1),phi(i,j-1),phi(i+1,j),phi(i-1,j), &
									   s_nn,s_ss,s_ee,s_ww,scheme, &
                                       uINT(i,j),uINT(i-1,j),vINT(i,j),vINT(i,j-1), &
                                       deltaT,deltaX,deltaY,xdirgo,ydirgo)
		enddo
    enddo
	! apply boundary conditions
    phi_tilda(:,1)=phi_tilda(:,2)
    phi_tilda(:,Ny+2)=phi_tilda(:,Ny+1)
    phi_tilda(1,:)=phi_tilda(2,:)
    phi_tilda(Nx+2,:)=phi_tilda(Nx+1,:)

	! advect density in direction 1
	do i=2,Nx+1
        do j=2,Ny+1
            Density(i,j) = advectVOF(A(i,j),B(i,j),C(i,j),A(i,j+1),B(i,j+1),C(i,j+1),A(i,j-1),B(i,j-1),C(i,j-1),&
									 A(i+1,j),B(i+1,j),C(i+1,j),A(i-1,j),B(i-1,j),C(i-1,j),F(i,j), &
                                     phi(i,j),phi(i,j+1),phi(i,j-1),phi(i+1,j),phi(i-1,j), &
                                     u(i,j),u(i-1,j),v(i,j),v(i,j-1),deltaT,deltaX,deltaY,&
									 xdirgo,ydirgo,xwhole(i),ywhole(j),rho_one,rho_two)
        enddo
    enddo
	! apply boundary conditions
    Density(:,1)=Density(:,2)
    Density(:,Ny+2)=Density(:,Ny+1)
    Density(1,:)=Density(2,:)
    Density(Nx+2,:)=Density(Nx+1,:)
	
    ! get volume fractions
    do i=1,Nx+2
        do j=1,Ny+2
			F_tilda(i,j) = (Density(i,j)-rho_one)/(rho_two-rho_one)

			! clip vol fracs
			if (F_tilda(i,j)>1.0) then
				F_tilda(i,j)=1.0
			elseif (F_tilda(i,j)<0.0) then
				F_tilda(i,j)=0.0
			endif
        enddo
    enddo

	! find the interface that satisfies F_tilda
    do i=2,Nx+1
        do j=2,Ny+1          
            call ReconInt(A(i,j),B(i,j),C(i,j),Res,F_tilda(i,j),F_tilda(i+1,j+1),F_tilda(i,j+1),F_tilda(i-1,j+1), &
                          F_tilda(i+1,j),F_tilda(i-1,j),F_tilda(i+1,j-1),F_tilda(i,j-1),F_tilda(i-1,j-1), &
                          phi_tilda(i,j),phi_tilda(i+1,j),phi_tilda(i-1,j),phi_tilda(i,j+1),phi_tilda(i,j-1), &
                          deltaX,deltaY,xwhole(i),ywhole(j),tol)
		enddo
    enddo


	! construct a signed distance function from the interface
    do i=2,Nx+1
        do j=2,Ny+1
            if (F_tilda(i,j)<=1.0-tol .and. F_tilda(i,j)>=tol) then
				if ( abs(phi_tilda(i,j))<=1.2*deltaX ) then
					phi_tilda(i,j) = A(i,j)
				endif
            endif
			H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,phi_tilda)
        enddo
    enddo
    deltaArtT = deltaX/2.0
    do iter = 1,7		!****************************
        do i=2,Nx+1
            do j=2,Ny+1
                uchar = 2.0*(H(i,j)-0.5)
				myNorm = NORMDIFFPHI(Nx,Ny,i,j,deltaX,deltaY,uchar,phi_tilda)
                phi_dummy(i,j) = phi_tilda(i,j) - deltaArtT*uchar*(myNorm-1.0)
            enddo
        enddo
        phi_tilda=phi_dummy
		phi_tilda(:,1)=phi_tilda(:,2)
		phi_tilda(:,Ny+2)=phi_tilda(:,Ny+1)
		phi_tilda(1,:)=phi_tilda(2,:)
		phi_tilda(Nx+2,:)=phi_tilda(Nx+1,:)
    enddo 



	! ---------------------- ADVECT DIRECTION 2 ----------------------

    ! Advect phi in direction 2
    do i=2,nx+1
        do j=2,ny+1
			if (i==2 .or. i==Nx+1 .or. j==2 .or. j==Ny+1) then
				s_nn = 0.0
				s_ss = 0.0
				s_ee = 0.0
				s_ww = 0.0
				scheme = 1  ! do first-order donor cell
			else
				s_nn = phi_tilda(i,j+2)
				s_ss = phi_tilda(i,j-2)
				s_ee = phi_tilda(i+2,j)
				s_ww = phi_tilda(i-2,j)
				scheme = 1  ! do higher-order scheme *****************************
			endif
            phi(i,j) = advectval(phi_tilda(i,j),phi_tilda(i,j+1),phi_tilda(i,j-1),&
								 phi_tilda(i+1,j),phi_tilda(i-1,j), &
								 s_nn,s_ss,s_ee,s_ww,scheme, &
                                 uINT(i,j),uINT(i-1,j),vINT(i,j),vINT(i,j-1), &
                                 deltaT,deltaX,deltaY,xdirgo,ydirgo)
        enddo
    enddo
	! apply boundary conditions
    phi(:,1)=phi(:,2)
    phi(:,Ny+2)=phi(:,Ny+1)
    phi(1,:)=phi(2,:)
    phi(Nx+2,:)=phi(Nx+1,:)

	! Advect density in direction 2
    do i=2,Nx+1
        do j=2,Ny+1
              Density(i,j) = advectVOF(A(i,j),B(i,j),C(i,j),A(i,j+1),B(i,j+1),C(i,j+1),A(i,j-1),B(i,j-1),C(i,j-1),&
									   A(i+1,j),B(i+1,j),C(i+1,j),A(i-1,j),B(i-1,j),C(i-1,j),F_tilda(i,j), &
                                       phi_tilda(i,j),phi_tilda(i,j+1),phi_tilda(i,j-1),phi_tilda(i+1,j),phi_tilda(i-1,j),& 
									   u(i,j),u(i-1,j),v(i,j),v(i,j-1),deltaT,deltaX,deltaY,&
									   xdirgo,ydirgo,xwhole(i),ywhole(j),rho_one,rho_two)
        enddo
    enddo
    ! clipping vol fracs
    do i=2,Nx+1
        do j=2,Ny+1
			F(i,j) = (Density(i,j)-rho_one)/(rho_two-rho_one)

			! clip vol fracs
			if (F(i,j)>1.0) then
				F(i,j)=1.0
			elseif (F(i,j)<0.0) then
				F(i,j)=0.0
			endif
        enddo
    enddo
	
	! apply boundary conditions
    F(:,1)=F(:,2)
    F(:,Ny+2)=F(:,Ny+1)
    F(1,:)=F(2,:)
    F(Nx+2,:)=F(Nx+1,:)

   
    
    ! find phiRec(A,B,C) that satisfy F
    do i=2,Nx+1
        do j=2,Ny+1
            call ReconInt(A(i,j),B(i,j),C(i,j),Res,F(i,j),F(i+1,j+1),F(i,j+1),F(i-1,j+1), &
                          F(i+1,j),F(i-1,j),F(i+1,j-1),F(i,j-1),F(i-1,j-1), &            
                          phi(i,j),phi(i+1,j),phi(i-1,j),phi(i,j+1),phi(i,j-1), &
                          deltaX,deltaY,xwhole(i),ywhole(j),tol)
		enddo
    enddo


	! construct a signed distance function from the interface
    do i=2,Nx+1
        do j=2,Ny+1
            if (F(i,j)<=1.0-tol .and. F(i,j)>=tol) then
				if ( abs(phi(i,j))<=1.2*deltaX ) then
					phi(i,j) = A(i,j)
				endif
            endif
        	H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,phi)
        enddo
    enddo
    deltaArtT = deltaX/2.0
    do iter = 1,7		!**************************************
        do i=2,Nx+1
            do j=2,Ny+1
                uchar = 2.0*(H(i,j)-0.5)
                myNorm = NORMDIFFPHI(Nx,Ny,i,j,deltaX,deltaY,uchar,phi)
                phi_dummy(i,j) = phi(i,j) - deltaArtT*uchar*(myNorm-1.0)
            enddo
        enddo
        phi=phi_dummy
		phi(:,1)=phi(:,2)
		phi(:,Ny+2)=phi(:,Ny+1)
		phi(1,:)=phi(2,:)
		phi(Nx+2,:)=phi(Nx+1,:)
    enddo    



ENDSUBROUTINE
!********************************************************************
!*																	*
!*							   advectVOF							*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a function for advecting the volume fraction in a		*
!* specified direction.												*
!*																	*	
!********************************************************************
REAL(kind=8) FUNCTION advectVOF(A_c,B_c,C_c,A_n,B_n,C_n,A_s,B_s,C_s, &
								A_e,B_e,C_e,A_w,B_w,C_w, &
								F_c,phi_c,phi_n,phi_s,phi_e,phi_w, &
								u_e,u_w,v_n,v_s,deltaT,deltaX,deltaY,&
								xdirgo,ydirgo,x_c,y_c,rho_one,rho_two)

	implicit none

	!----------------------------------- HEADER -----------------------------------
	! integers
	integer, intent(in) :: xdirgo,ydirgo
	
	! reals
	real(kind=8), intent(in) :: A_c,B_c,C_c,A_n,B_n,C_n,A_s,B_s,C_s, &
								A_e,B_e,C_e,A_w,B_w,C_w, &
								F_c,phi_c,phi_e,phi_w,phi_n,phi_s, &
								u_e,u_w,v_n,v_s,deltaT,deltaX,deltaY,x_c,y_c,&
								rho_one,rho_two

	! local reals
	real(kind=8) :: Foneface,Ftwoface,G_one,G_two,vel_one,vel_two,delta,&
					rho_c,rho_e,rho_w,rho_n,rho_s
	!------------------------------------------------------------------------------

	if (xdirgo==1) then
		vel_one = u_e		! ueast
		vel_two = u_w		! uwest
        ! find the volume fluxes     
		call volfracflux(Foneface,Ftwoface,x_c,y_c,A_c,B_c,C_c,A_n,B_n,C_n,A_s,B_s,C_s, &
						 A_e,B_e,C_e,A_w,B_w,C_w, &
					 	 u_e,u_w,v_n,v_s,deltaT,deltaX,deltaY,xdirgo,ydirgo)
		rho_e = rho_one + Foneface*(rho_two-rho_one)
		rho_w = rho_one + Ftwoface*(rho_two-rho_one)
		G_one = vel_one*rho_e
		G_two = vel_two*rho_w
		delta = deltaX
	elseif (ydirgo==1) then
		vel_one = v_n	! vnorth
		vel_two = v_s	! vsouth
		! find the volume fluxes
		call volfracflux(Foneface,Ftwoface,x_c,y_c,A_c,B_c,C_c,A_n,B_n,C_n,A_s,B_s,C_s, &
						 A_e,B_e,C_e,A_w,B_w,C_w, &
						 u_e,u_w,v_n,v_s,deltaT,deltaX,deltaY,xdirgo,ydirgo)
		rho_n = rho_one + Foneface*(rho_two-rho_one)
		rho_s = rho_one + Ftwoface*(rho_two-rho_one)
		G_one = vel_one*rho_n
		G_two = vel_two*rho_s
		delta = deltaY
	endif
	   	
	rho_c = rho_one + F_c*(rho_two-rho_one)	! the density of the cell

	! the density in this cell at the next time level
	advectVOF = rho_c - deltaT/delta*(G_one-G_two)

ENDFUNCTION
!********************************************************************
!*																	*
!*							   advectval							*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a function for advecting a scalar quantity in a 			*
!* specified direction.												*
!*																	*	
!********************************************************************
REAL(kind=8) FUNCTION advectval(s_c,s_n,s_s,s_e,s_w,s_nn,s_ss,s_ee,s_ww,&
								scheme,u_e,u_w,v_n,v_s,&
								deltaT,deltaX,deltaY,xdirgo,ydirgo)
	
	implicit none

	!----------------------------------- HEADER -----------------------------------
	! integers
	integer, intent(in) :: xdirgo,ydirgo,scheme

	! reals
	real(kind=8), intent(in) :: s_c,s_n,s_s,s_e,s_w,s_nn,s_ss,s_ee,s_ww,&
								u_e,u_w,v_n,v_s, &
								deltaT,deltaX,deltaY
	
	! local integers
	real(kind=8) :: vel_one,vel_two,G_one,G_two,delta,sface
	!------------------------------------------------------------------------------


	if (xdirgo==1) then
		vel_one = u_e		!ueast
		vel_two = u_w		!uwest
        
		if (scheme==1) then
			! first-order upwinding     
			if (u_e>0.0) then
				G_one = s_c*vel_one		! Geast
			else
				G_one = s_e*vel_one		! Geast
			endif

			if (u_w>0.0) then
				G_two = s_w*vel_two		! Gwest
			else
				G_two = s_c*vel_two		! Gwest
			endif
		else
			! second-order method
			if (u_e>0.0) then
				sface = s_c + deltaX/2.0*(1.0-vel_one*deltaT/deltaX)*(s_e-s_w)/deltaX
				G_one = sface*vel_one	! Geast
			else
				sface = s_e - deltaX/2.0*(1.0+vel_one*deltaT/deltaX)*(s_ee-s_c)/deltaX
				G_one = sface*vel_one	! Geast
			endif

			if (u_w>0.0) then
				sface = s_w + deltaX/2.0*(1.0-vel_two*deltaT/deltaX)*(s_c-s_ww)/deltaX
				G_two = sface*vel_two	! Gwest
			else
				sface = s_c - deltaX/2.0*(1.0+vel_two*deltaT/deltaX)*(s_e-s_w)/deltaX
				G_two = sface*vel_two	! Gwest
			endif
		endif

		delta = deltaX
	elseif (ydirgo==1) then
		vel_one = v_n	! vnorth
		vel_two = v_s	! vsouth
        
		if (scheme==1) then         
			! first-order upwinding   
			if (v_n>0.0) then
				G_one = s_c*vel_one		! Gnorth
			else
				G_one = s_n*vel_one		! Gnorth
			endif

			if (v_s>0.0) then
				G_two = s_s*vel_two		! Gsouth 
			else
				G_two = s_c*vel_two		! Gsouth
			endif
		else
			! second-order method
			! first-order upwinding   
			if (v_n>0.0) then
				sface = s_c + deltaY/2.0*(1.0-vel_one*deltaT/deltaY)*(s_n-s_s)/deltaY
				G_one = sface*vel_one	! Gnorth
			else
				sface = s_n - deltaY/2.0*(1.0+vel_one*deltaT/deltaY)*(s_nn-s_c)/deltaY
				G_one = sface*vel_one	! Gnorth
			endif


			if (v_s>0.0) then
				sface = s_s + deltaY/2.0*(1.0-vel_two*deltaT/deltaY)*(s_c-s_ss)/deltaY
				G_two = sface*vel_two	! Gsouth 
			else
				sface = s_c - deltaY/2.0*(1.0+vel_two*deltaT/deltaY)*(s_n-s_s)/deltaY
				G_two = sface*vel_two	! Gsouth
			endif
		endif

		delta = deltaY
	endif            
	advectval = s_c - deltaT/delta*(G_one-G_two)

ENDFUNCTION

!********************************************************************
!*																	*
!*							ReconInt								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a subroutine for reconstructing the interface			*
!*																	*	
!********************************************************************
SUBROUTINE ReconInt(A,B,C,Res,F_c,F_ne,F_n,F_nw,F_e,F_w,F_se,F_s,F_sw, & 
                    phi_c,phi_e,phi_w,phi_n,phi_s,deltaX,deltaY,x_c,y_c,tol)

	implicit none

	!----------------------------------- HEADER -----------------------------------
	! reals
	real(kind=8), intent(inout) :: A,B,C,Res
	real(kind=8), intent(in) :: F_c,F_ne,F_n,F_nw,F_e,F_w,F_se,F_s,F_sw, &
								phi_c,phi_e,phi_w,phi_n,phi_s,deltaX,deltaY,x_c,y_c,tol 
	! local integers
	integer :: N,i,iter

	! local reals
	real(kind=8) :: alpha_1,alpha_2,beta_1,beta_2,mag,x1,x2,xnew,deltaA,Amin,Amax,Res1,Res2,&
					lx,ly,Adummy,ResOld,L,psi,theta
	
	! parameter
	parameter(N=4)

	! local arrays
	real(kind=8), allocatable :: Res_A(:), Aguess(:)

	! externals
	real(kind=8), external :: ResFunc
	
	allocate( Res_A(N) )
	allocate( Aguess(N) )
	!------------------------------------------------------------------------------

	
	! Find the A value that satisfies VOF
	alpha_1 = x_c - deltaX/2.0
	alpha_2 = x_c + deltaX/2.0
	beta_1 = y_c - deltaY/2.0
	beta_2 = y_c + deltaY/2.0

	if (F_c<1.0 .and. F_c>0.0) then ! cell contains interface so construct interface
    
		! construct interface
		A = phi_c
		B = (phi_e-phi_w)/(2.0*deltaX)	! normal x-direction
		C = (phi_n-phi_s)/(2.0*deltaY)	! normal y-direction
		! F_ne,F_n,F_nw,F_e,F_c,F_w,F_se,F_s,F_sw
		!B = -(F_ne + 2.0*F_e + F_se - F_nw - 2.0*F_w - F_sw)/deltaX
		!C = -(F_ne + 2.0*F_n + F_nw - F_se - 2.0*F_s - F_sw)/deltaY
                
		mag = sqrt(B**2.0 + C**2.0)
		A = A/mag
		B = B/mag
		C = C/mag
                
                
		x1 = A
		Res = ResFunc(x1,B,C,F_c,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)


		! the exact Amin and Amax
		L = 0.5*sqrt(deltaX**2.0 + deltaY**2.0)
		psi = atan(deltaY/deltaX)
		theta = atan(abs(C/B))
		Amax = 1.0001*L*abs(cos(psi-theta))	! I need extra coef to overcome round-off errors
		Amin = -Amax

	
		deltaA = (Amax-Amin)/real(N-1)
		Aguess(1) = Amin
		Aguess(N) = Amax
		do i=2,N-1
			Aguess(i) = Aguess(i-1) + deltaA
		enddo

		do i=1,N
			res_A(i) = ResFunc(Aguess(i),B,C,F_c,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)
		enddo
		do i=1,N-1
			if (res_A(i)*res_A(i+1)<=0.0) then
				x1 = Aguess(i)
				x2 = Aguess(i+1)
			endif
		enddo

		Res1 = ResFunc(Amin,B,C,F_c,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)
		Res2 = ResFunc(Amax,B,C,F_c,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)
		if (Res1*Res2>0.0) then
			print *, Res1, Res2
			pause
			if (abs(Res)>tol) then
				if (abs(Res1)<abs(Res2)) then
					A = Amin
					Res = Res1
				else
					A = Amax
					Res = Res2
				endif
			endif
		endif
		
		if (abs(Res)>tol) then
			call zbrent(A,Res,x1,x2,tol,B,C,F_c,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c) 
		endif
          
	else ! cell does not contain interface so do not construct interface
		A = phi_c
		B = (phi_e-phi_w)/(2.0*deltaX)
		C = (phi_n-phi_s)/(2.0*deltaY)  
		Res = 0.0

	endif ! end of if statement

ENDSUBROUTINE
!********************************************************************
!*																	*
!*							   ResFunc								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a function for calculating the volume fraction residual	*
!*																	*	
!********************************************************************
REAL(kind=8) FUNCTION ResFunc(A,B,C,F,alpha_1,alpha_2,beta_1,beta_2,&
							  deltaX,deltaY,xcenter,ycenter)
	
	implicit none

	!----------------------------------- HEADER -----------------------------------
	! reals
	real(kind=8), intent(in) :: A,B,C,F,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,xcenter,ycenter
	
	! local reals
	real(kind=8) :: volfrac

	! externals
	real(kind=8), external :: volfraction
	!------------------------------------------------------------------------------


	volfrac = volfraction(A,B,C,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,xcenter,ycenter)
	ResFunc = volfrac-F

ENDFUNCTION
!********************************************************************
!*																	*
!*							volfracflux								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a subroutine for calculating the volume fraction flux	*
!* from a cell														*
!*																	*	
!********************************************************************
SUBROUTINE volfracflux(Foneface,Ftwoface,x_c,y_c,A_c,B_c,C_c,A_n,B_n,C_n,A_s,B_s,C_s, &
					   A_e,B_e,C_e,A_w,B_w,C_w, &
					   u_e,u_w,v_n,v_s,deltaT,deltaX,deltaY,xdirgo,ydirgo)

	implicit none

	!----------------------------------- HEADER -----------------------------------
	! integers
	integer, intent(in) :: xdirgo,ydirgo
	
	! reals
	real(kind=8), intent(in) :: x_c,y_c,A_c,B_c,C_c,A_n,B_n,C_n,A_s,B_s,C_s, &
								A_e,B_e,C_e,A_w,B_w,C_w, &
								u_e,u_w,v_n,v_s,deltaT,deltaX,deltaY
	real(kind=8), intent(out) :: Foneface,Ftwoface

	! local reals
	real(kind=8) :: alpha_1,alpha_2,beta_1,beta_2,G_one,G_two,vel_one,vel_two,delta,&
					volfrac,vol,A,B,C,x,y

	! externals
	real(kind=8), external :: volfraction
	!------------------------------------------------------------------------------



	if (xdirgo==1) then
		beta_1 = y_c-deltaY/2.0
		beta_2 = y_c+deltaY/2.0
    
		! ----- EAST FACE -----
		vel_one = u_e		! ueast
		if (vel_one>0.0) then
			! flux is going out east face
			alpha_1 = x_c+deltaX/2.0 - vel_one*deltaT
			alpha_2 = x_c+deltaX/2.0

			A = A_c
			B = B_c
			C = C_c
			x = x_c
			y = y_c
		elseif (vel_one<0.0) then
			! flux is comming in east face
			alpha_1 = x_c+deltaX/2.0
			alpha_2 = x_c+deltaX/2.0 - vel_one*deltaT

			A = A_e !A_e
			B = B_e !B_e
			C = C_e !C_e
			x = x_c+deltaX
			y = y_c
		else
			! no flux
			alpha_1 = x_c
			alpha_2 = x_c
		endif          
    
		if (abs(alpha_2-alpha_1)/=0.0) then
			volfrac = volfraction(A,B,C,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x,y) ! volume fraction
			Foneface = volfrac ! the volume fraction entering/leaving
		else
			Foneface = 0.0		! velocity is zero so flux is zero
		endif
    
		! ----- WEST FACE -----
		vel_two = u_w	! uwest
		if (vel_two>0.0) then
			! flux is going in west face
			alpha_1 = x_c-deltaX/2.0 - vel_two*deltaT
			alpha_2 = x_c-deltaX/2.0

			A = A_w ! A_w
			B = B_w ! B_w
			C = C_w ! C_w
			x = x_c-deltaX
			y = y_c
		elseif (vel_two<0.0) then
			! flux is going out west face
			alpha_1 = x_c-deltaX/2.0
			alpha_2 = x_c-deltaX/2.0 - vel_two*deltaT

			A = A_c
			B = B_c
			C = C_c
			x = x_c
			y = y_c
		else
			! no flux
			alpha_1 = x_c
			alpha_2 = x_c
		endif          
    
		if (abs(alpha_2-alpha_1)/=0.0) then
			volfrac = volfraction(A,B,C,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x,y) ! volume fraction
			Ftwoface = volfrac ! the volume fraction entering/leaving
		else
			Ftwoface = 0.0	! velocity is zero so flux is zero
		endif
    
	elseif (ydirgo==1) then
    
		alpha_1 = x_c-deltaX/2.0
		alpha_2 = x_c+deltaX/2.0
    
		! NORTH FACE
		vel_one = v_n		! vnorth
    
		if (vel_one>0.0) then
			! flux is going out north face
			beta_1 = y_c+deltaY/2.0-vel_one*deltaT
			beta_2 = y_c+deltaY/2.0
			
			A = A_c
			B = B_c
			C = C_c
			x = x_c
			y = y_c
		elseif (vel_one<0.0) then
			! flux is going in north face
			beta_1 = y_c+deltaY/2.0
			beta_2 = y_c+deltaY/2.0-vel_one*deltaT

			A = A_n !A_n
			B = B_n !B_n
			C = C_n !C_n
			x = x_c
			y = y_c+deltaY
		else
			! no flux
			beta_1 = y_c
			beta_2 = y_c
		endif

		if (abs(beta_2-beta_1)/=0.0) then
			volfrac = volfraction(A,B,C,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x,y) ! volume fraction
			Foneface = volfrac ! the volume fraction entering/leaving
		else
			Foneface = 0.0	! velocity is zero so flux is zero
		endif 
    
		! SOUTH FACE
		vel_two = v_s		! vsouth
		if (vel_two>0.0) then
			! flux is going in south face
			beta_1 = y_c-deltaY/2.0-vel_two*deltaT
			beta_2 = y_c-deltaY/2.0
			
			A = A_s !A_s
			B = B_s !B_s
			C = C_s !C_s
			x = x_c
			y = y_c-deltaY
		elseif (vel_two<0.0) then
			! flux is going out south face
			beta_1 = y_c-deltaY/2.0
			beta_2 = y_c-deltaY/2.0-vel_two*deltaT

			A = A_c
			B = B_c
			C = C_c
			x = x_c
			y = y_c
		else
			! no flux
			beta_1 = y_c
			beta_2 = y_c
		endif

		if (abs(beta_2-beta_1)/=0) then
			volfrac = volfraction(A,B,C,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x,y) 
			Ftwoface = volfrac ! the volume fraction entering/leaving
		else
			Ftwoface = 0.0	! velocity is zero so flux is zero
		endif 
 
	endif  


ENDSUBROUTINE
!********************************************************************
!*																	*
!*							volfraction								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a function for calculating the volume fraction in a cell *
!*																	*	
!********************************************************************
REAL(kind=8) FUNCTION volfraction(A,B,C,alpha_1,alpha_2,beta_1,beta_2,&
								  deltaX,deltaY,x_c,y_c)
                             
	implicit none

	!----------------------------------- HEADER -----------------------------------
	! reals
	real(kind=8), intent(in) :: A,B,C,alpha_1,alpha_2,beta_1,beta_2,&
								deltaX,deltaY,x_c,y_c

	! local integers
	integer :: k,v

	! local reals
	real(kind=8) :: distance,delta_beta,delta_alpha,&
					phibar_sw,phibar_se,phibar_ne,phibar_nw,&
					x_s,y_e,x_n,y_w,Vol,dummy,deltaX0,deltaY0
					
	! local real arrays
	real(kind=8), dimension(8) :: x,y
	!------------------------------------------------------------------------------


	distance = 1.0
	volfraction = 0.0

	dummy = 0.0

	if (abs(A)<=distance*max(deltaX,deltaY))  then
		phibar_sw = A + B*(alpha_1-x_c) + C*(beta_1-y_c)
		phibar_se = A + B*(alpha_2-x_c) + C*(beta_1-y_c)
		phibar_ne = A + B*(alpha_2-x_c) + C*(beta_2-y_c)
		phibar_nw = A + B*(alpha_1-x_c) + C*(beta_2-y_c)

		x_s = -A/(B+1.e-16) - C/(B+1.e-16)*(beta_1-y_c) + x_c
		y_e = -A/(C+1.e-16) - B/(C+1.e-16)*(alpha_2-x_c) + y_c
		x_n = -A/(B+1.e-16) - C/(B+1.e-16)*(beta_2-y_c) + x_c
		y_w = -A/(C+1.e-16) - B/(C+1.e-16)*(alpha_1-x_c) + y_c


		k=1
		! sw corner
		if (phibar_sw<=dummy) then
			x(k)=alpha_1
			y(k)=beta_1
			k=k+1
		endif
		! south face
		if (x_s>alpha_1 .and. x_s<alpha_2) then
			x(k)=x_s
			y(k)=beta_1
			k=k+1
		endif
		! se corner
		if (phibar_se<=dummy) then
			x(k) = alpha_2
			y(k) = beta_1
			k=k+1
		endif
		! east face
		if (y_e>beta_1 .and. y_e<beta_2) then
			x(k) = alpha_2
			y(k) = y_e
			k=k+1
		endif
		! ne corner
		if (phibar_ne<=dummy) then
			x(k) = alpha_2
			y(k) = beta_2
			k=k+1
		endif
		! north face
		if (x_n>alpha_1 .and. x_n<alpha_2) then
			x(k) = x_n
			y(k) = beta_2
			k=k+1
		endif
		! nw corner
		if  (phibar_nw<=dummy) then
	        x(k) = alpha_1
		    y(k) = beta_2
			k=k+1
		endif
		! west face
		if (y_w>beta_1 .and. y_w<beta_2) then
			x(k) = alpha_1
			y(k) = y_w
			k=k+1
		endif
		x(k) = x(1)
		y(k) = y(1)
    
	    Vol = 0.0
		do v=1,k-1
			Vol = Vol + 1.0/2.0*(x(v)*y(v+1) - x(v+1)*y(v))
	    enddo
		volfraction = Vol/( (alpha_2-alpha_1)*(beta_2-beta_1) )
		
		! in the event the interface is not in the cell
		if (phibar_sw<=0.0 .and. phibar_se<=0.0 .and. phibar_ne<=0.0 .and. phibar_nw<=0.0) then
			volfraction = 1.0
		elseif (phibar_sw>=0.0 .and. phibar_se>=0.0 .and. phibar_ne>=0.0 .and. phibar_nw>=0.0) then
			volfraction = 0.0
		endif


	elseif (A>distance*max(deltaX,deltaY)) then
		volfraction = 0.0
	else
		volfraction = 1.0 
    
	endif ! end of if statement

ENDFUNCTION
!********************************************************************
!*																	*
!*							  zbrent								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for calculating the route of a nonlinear		*
!* function, FCN=0													*
!*																	*	
!********************************************************************
SUBROUTINE zbrent(root,Res,x1,x2,tol,B_coef,C_coef,F,&
				  alpha_1,alpha_2,beta_1,beta_2,&
				  deltaX,deltaY,x_c,y_c)

	implicit none
	
	!----------------------------------- HEADER -----------------------------------
	! reals
	real(kind=8), intent(in) :: x1,x2,tol,B_coef,C_coef,F,alpha_1,alpha_2,beta_1,beta_2,&
								deltaX,deltaY,x_c,y_c
	real(kind=8), intent(out) :: root,Res

	! local integers
	integer :: itermax,iter

	! local reals
	real(kind=8) :: EPS,a,b,c,d,e,fa,fb,fc,xm,s,q,r,tol1,p

	! externals
	real(kind=8), external :: ResFunc 
	!------------------------------------------------------------------------------


	itermax = 40
	EPS = 3.e-8

	a=x1
	b=x2
	fa = ResFunc(a,B_coef,C_coef,F,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)
	fb = ResFunc(b,B_coef,C_coef,F,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)

	if ( (fa>0.0 .and. fb>0.0) .or. (fa<0.0 .and. fb<0.0) ) then
		print *, "Error - root must be bracted"
		iter=0
	endif

	c=b
	fc=fb

	do iter=1,itermax
		if ( (fb>0.0 .and. fc>0.0) .or. (fb<0.0 .and. fc<0.0) ) then    ! Rename a,b,c, and adjust bounding interval d
			c=a
			fc=fa
			d=b-a
			e=d
		endif
		if (abs(fc) < abs(fb)) then
			a=b
			b=c
			c=a
			fa=fb
			fb=fc
			fc=fa
		endif
		tol1=2.0*EPS*abs(b) + 0.5*tol    ! Convergence check
		xm=0.5*(c-b)
		if (abs(xm)<=tol1 .or. fb==0) then
			root=b
			Res=xm
			return
		endif
		if ( abs(e)>= tol1 .and. abs(fa)>abs(fb) ) then
			s=fb/fa        ! Attempt inverse quadratic interpolation
			if (a==c) then
				p=2*xm*s
				q=1-s
			else
				q=fa/fc
				r=fb/fc
				p=s*(2*xm*q*(q-r) - (b-a)*(r-1))
				q=(q-1)*(r-1)*(s-1)
			endif
			if (p>0) then 
				q=-q       ! Check whether in bounds
			endif
			p=abs(p)
			if ( 2*p<min(3*xm*q - abs(tol1*q), abs(e*q)) ) then
				e=d        ! Accept interpolation
				d=p/q
			else
				d=xm       ! Interpolation failed, use bisection
				e=d
			endif
		else               ! Bounds decreasing too slowly, use bisection
			d=xm
			e=d
		endif
		a=b                ! move last best guess to a
		fa=fb
		if( abs(d)>tol1 )  then ! Evaluate new trial root
			b=b+d
		else
			b=b+sign(tol1,xm)
		endif
		fb=ResFunc(b,B_coef,C_coef,F,alpha_1,alpha_2,beta_1,beta_2,deltaX,deltaY,x_c,y_c)
	enddo
    root=b	
	Res=xm

ENDSUBROUTINE

