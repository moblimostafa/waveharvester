!********************************************************************
!*                                                                  *
!*                                                          MPXLIB  *
!*                                                                  *
!********************************************************************
!* Author: Nathaniel R. Morgan (2004)                               *
!*                                                                  *
!* MPXLIB is a (M)ulti(P)hase dynamics(X) (LIB)ary for modeling     *
!* 2-phase flow problems.  MPXLIB solves the full transient         *
!* Navier-Stokes equations and has the ability to capture the       *
!* interface between two fluids of different properties using the   *
!* Level Set (LS) method.                                           *    
!*                                                                  *
!*                                                                  *
!* The staggered grid is shown below:                               *
!*                                                                  *
!* (bottom left corner)                                             *
!*                                                                  *
!* v----o----v----o  y(2)                                           *
!*      |         |                                                 *
!* P    u    P    u                                                 *
!*      |         |                                                 *
!* v----o----v----o  y(1)                                           *
!*            |         |                                           *
!* P    u    P    u                                                 *
!*                                                                  *
!*     x(1)     x(2)                                                *
!*                                                                  *
!********************************************************************

! NOTES:
! Modified by Jessica Sanders Summer, 2008
! Modified by Curtis Lee, Summer 2010
PROGRAM MPXLIB

        IMPLICIT NONE
        !
        ! INTEGERS 
        ! JDS 7/17/08 Removing temperature variables
        integer :: Nx,Ny,i,j,step,r,l,k                             ! loop integers
        integer :: ios                                              ! input file
        integer :: file_count,count,dumpcount                       ! counting integers
        integer :: iter,itermax,iterprint                           ! iteration integers, not used
        integer :: MBDCND,N,M,NBDCND,IDIMF,IERROR                   ! for fish pack solver
        integer :: solvertype,precond
        integer :: IC
        integer :: BCnorth,BCsouth,BCeast,BCwest, &
                           BCPnorth,BCPsouth,BCPeast,BCPwest, &
                           BCphisouth,BCphinorth,BCphiwest,&
                           BCphieast                                ! BC types - wall,flow,etc
        integer :: rb_shape                                         ! Shape of the rigid body
                                                                    ! 1 = circle
                                                                    ! 2 = ellipse
                                                                    ! 3 = square
        integer :: int_type                                         ! Type of internal component(s)
                                                                    ! 0 = none
                                                                    ! 1 = torsional disc
                                                                    ! 2 = 2D trans
                                                                    ! 3 = 1D trans
        integer :: RB_lock                                          ! Restrain motion of R
                                                                    ! 0 = No restraints
                                                                    ! 1 = Vertical motion only
        integer :: aom_counter                                      ! For detecting RB internal instability
        !
        !
        ! LOGICAL AND CHARACTER
        logical :: startup,trackphase,immersedbc                    ! booleans        
        logical :: rigidbody
        logical :: wedge                                            ! Wedge
        logical :: inside
        character(20) :: filename,variables,startfile               ! strings
        !
        !
        ! REALS
        real(kind=8) :: pi                                          ! pi
        real(kind=8) :: Lx,Ly                                       ! domain size
        real(kind=8) :: umax,vmax,dtvisc,dtadv,&
                        dtg,dtsigma,dtinternal                      ! time step variables
        real(kind=8) :: deltaT,deltaX,deltaY,deltaTau               ! incremenit
        real(kind=8) :: time,time_max, &
                        time_print,time_count                       ! time variables
        real(kind=8) :: stability
        real(kind=8) :: rho_one,rho_two,rho_three,mu_one,&        
                        mu_two,mu_three,sigma                       ! properties
        real(kind=8) :: kappa_n,kappa_s,kappa_w,kappa_e             ! interface geometry
        real(kind=8) :: gx,gy,Fsigma                                ! body forces
        real(kind=8) :: Tolerance                                   ! iteration variables
        real(kind=8) :: unorth,usouth,uwest,ueast,&
                        vnorth,vsouth,vwest,veast                                
        real(kind=8) :: Pnorth,Psouth,Pwest,Peast                   ! dirclet pressure BC's
        real(kind=8) :: contC,contA,contR,contM                     ! contact line BC's
        real(kind=8) :: u_E,u_W,u_N,u_S,&
                        u_ne,u_se,u_nw,u_sw,u_O, &
                        v_E,v_W,v_N,v_S, &
                        v_ne,v_nw,v_se,v_sw,v_O                     ! velocities
        real(kind=8) :: u_xe,u_xw,u_yn,u_ys,&
                        v_xn,v_xs,u_xO,u_yO                         ! velocity derivatives
        real(kind=8) :: u_ye,u_yw,v_yn,v_ys,&
                        v_xe,v_xw,v_xO,v_yO                         ! velocity derivatives
        real(kind=8) :: mu_n,mu_s,mu_e,mu_w, &                        
                        rho_O,rho_n,rho_s,rho_e,rho_w               ! properties
        real(kind=8) :: u_plus,v_plus                               ! upwinding terms
        real(kind=8) :: u_minus,v_minus                             ! upwinding terms
        real(kind=8) :: dudx_plus,dudx_minus                        ! upwinding terms
        real(kind=8) :: dudy_plus,dudy_minus                        ! upwinding terms
        real(kind=8) :: dvdx_plus,dvdx_minus                        ! upwinding terms
        real(kind=8) :: dvdy_plus,dvdy_minus                        ! upwinding terms
        real(kind=8) :: Visc                                        ! viscous or diffusion term
        real(kind=8) :: DiffX, DiffY                                ! derivatives in the respective directions
        real(kind=8) :: A,B,C,D,PERTRB,ELMBDA                       ! for fishpac solver
        real(kind=8) :: an,as,aw,ae,aO,omega                        ! for SOR
        real(kind=8) :: Pavg,PavgOld,PRes,maxPRes                   ! pressure variables (SOR)
        real(kind=8) :: deltau,deltav                               ! change in u and v (immersed boundary)                                                        
        real(kind=8) :: xcenter,ycenter                             ! dummy variables for ICs
        real(kind=8) :: a_phiLS,b_phiLS,c_phiLS,d_phiLS,&
                        e_phiLS,f_phiLS,v_phiLS
        !
        real(kind=8) :: Coef(4)
        real(kind=8) :: normal(2)
        real(kind=8) :: total_energy                                ! Total system energy
        real(kind=8) :: lin_momx                                    ! Momentum in x-direction
        real(kind=8) :: lin_momy                                    ! Momentum in y-direction
        real(kind=8) :: net_mom                                     ! Net momentum

        ! Rigid Body Variables
        real(kind=8) :: a_rbLS, b_rbLS, c_rbLS, d_rbLS              ! Rigid body geometry
        
        real(kind=8) :: rb1_mass,big_J                              ! Rigid body mass
        real(kind=8) :: rb_mass                                     ! TOTAL mass (RB + Internals)
        
        real(kind=8) :: n_rb_dx,n_rb_dy                             ! Rigid body displacement        
        real(kind=8) :: n1_rb_dx,n1_rb_dy
        real(kind=8) :: rx,ry
        
        real(kind=8) :: n_w_dy                                      ! Wedge velocities
        real(kind=8) :: n1_w_dy        
        real(kind=8) :: n12_w_vy
        real(kind=8) :: n1_w_vy        

        
        real(kind=8) :: n_rb_vx,n_rb_vy                             ! Rigid body velocities
        real(kind=8) :: n12_rb_vx,n12_rb_vy
        real(kind=8) :: n1_rb_vx,n1_rb_vy
        
        real(kind=8) :: n_rb_ax,n_rb_ay                             ! Rigid body acceleration        
        real(kind=8) :: n1_rb_ax,n1_rb_ay

        real(kind=8) :: big_theta       
        real(kind=8) :: th_change,th_domain,th_old,th_counter       ! Change in RB's angle (from original value)
        real(kind=8) :: n_rb_vom,n12_rb_vom                         ! Rotational velocities
        real(kind=8) :: n1_rb_vom
        real(kind=8) :: n_rb_aom,n1_rb_aom                          ! Rotational acceleration
        real(kind=8) :: xc,yc                                       ! Center of mass (Actually, geometric center... -CL)
        real(kind=8) :: u_12,v_12                                   ! Intermediate variables for setting
                                                                    ! fluid velocity to RB velocity
        real(kind=8) :: adj_factor

        ! Torsional Disc Variables               
        real(kind=8) :: small_theta_n, small_theta_n1        
        real(kind=8) :: n_rb2_vom,n12_rb2_vom                       ! Rotational velocities
        real(kind=8) :: n1_rb2_vom
        real(kind=8) :: n_rb2_aom,n1_rb2_aom                        ! Rotational acceleration
        
        real(kind=8) :: delta_theta                                 ! Angle difference between disc and RB
        real(kind=8) :: delta_vom                                   ! Angular velocity difference
        real(kind=8) :: int_torque                                  ! Internal torque 

        real(kind=8) :: disk_mass,small_J                            ! Mass/moment of disc
        real(kind=8) :: k12                                         ! Torsional spring constant
        real(kind=8) :: c12                                         ! Torsional damping constant
        real(kind=8) :: aom_limit                                   ! Used in detection of instability
      
        ! Slider Variables
        real(kind=8) :: n_rb2_dx,n_rb2_dy                           ! Displacement
        real(kind=8) :: n1_rb2_dx,n1_rb2_dy                         ! "px" = "pseudo x" - buoy's ref frame
        
        real(kind=8) :: n_rb2_vx,n_rb2_vy                           ! Velocities
        real(kind=8) :: n12_rb2_vx,n12_rb2_vy
        real(kind=8) :: n1_rb2_vx,n1_rb2_vy
        
        real(kind=8) :: n_rb2_ax,n_rb2_ay                           ! Acceleration        
        real(kind=8) :: n1_rb2_ax,n1_rb2_ay
                 
        real(kind=8) :: delta_dpx, delta_dpy                        ! Differences in displacement and velocity
        real(kind=8) :: delta_vpx,delta_vpy                        
        real(kind=8) :: int_forcex,int_forcey
        real(kind=8) :: int_forcepx,int_forcepy

        real(kind=8) :: slider_mass, slider_J                       ! Mass, MoI, and spring/damping constants
        real(kind=8) :: buoy_J,old_big_J,diff_J
        real(kind=8) :: scgx,scgy                                   ! Slider Center of Mass
        real(kind=8) :: kx12            
        real(kind=8) :: cx12                                     
        real(kind=8) :: ky12
        real(kind=8) :: cy12

        real(kind=8) :: max_e,max_w,max_n,max_s                     ! Spring constant adapts to keep mass in buoy
        real(kind=8) :: kx_adj,ky_adj                               ! Adjusted spring constants
        real(kind=8) :: ext_px,ext_py                               ! external forces in buoy reference frame


        ! Integral terms for the RHS of rigid body equations
        real(kind=8) :: n_xint_conv, n_xint_grav, n_xint_inert      ! For the first step the convective and gravity terms
                                                                    ! For the correction step, the inertial terms
        real(kind=8) :: n_yint_conv, n_yint_grav, n_yint_inert      ! For the first step the convective and gravity terms
                                                                    ! For the correction step the inertial term
        ! Integral terms for the RHS of rigid body equations
        real(kind=8) :: n1_xint_conv, n1_xint_grav                  ! For the first step the convective and gravity terms
        real(kind=8) :: n1_xint_inert                               ! For the correction step, the inertial terms
        real(kind=8) :: n1_yint_conv, n1_yint_grav                  ! For the first step the convective and gravity terms
        real(kind=8) :: n1_yint_inert                               ! For the correction step, the intertial terms
        
        ! Integral terms for calculating fluid forces on wedge
        
        real(kind=8) :: n_w_conv, n_w_grav                  
        real(kind=8) :: n_w_inert                     
        real(kind=8) :: n1_w_conv, n1_w_grav                
        real(kind=8) :: n1_w_inert
        real(kind=8) :: wedge_force
        real(kind=8) :: wedge_work        
        real(kind=8) :: int_work

        ! Integral terms for the RHS of rigid body equations
        real(kind=8) :: r_int_conv, r_int_grav                      ! For the first step the convective and gravity terms
        real(kind=8) :: r_int_inert                                 ! For the correction step, the inertial terms
        real(kind=8) :: r1_int_conv, r1_int_grav                    ! For the first step the convective and gravity terms
        real(kind=8) :: r1_int_inert                                ! For the correction step, the inertial terms
                                                                                                                
        real(kind=8) :: ustar_int, vstar_int                        ! n+1 convective dummy variables for integral claculations
        
        real(kind=8) :: aellipse, bellipse                          ! ellipse variables
        real(kind=8) :: xo, yo, th, xn, yn, tho
        real(kind=8) :: c1,c2        
        
        real(kind=8) :: x1,x2,x3,x4                                 ! square variables
        real(kind=8) :: y1,y2,y3,y4
        real(kind=8) :: nx1,nx2,nx3,nx4                             ! square variables
        real(kind=8) :: ny1,ny2,ny3,ny4
        real(kind=8) :: cgx,cgy,ncgx,ncgy                           ! center of mass (CL)
        real(kind=8) :: netcgx,netcgy                               ! net center of mass (buoy + internals)

        real(kind=8) :: cost,sint
        real(kind=8) :: Det,a_det,b_det,c_det,temp        
        real(kind=8) :: c_xo,c_yo,c_tho
        
        !
        ! ALLCOCATABLE ARRAYS
        real(kind=8), ALLOCATABLE :: x(:),y(:)                       ! coordinates
        real(kind=8), ALLOCATABLE :: phi(:,:),phiLS(:,:)             ! interface variables
        real(kind=8), allocatable :: phiLSn(:,:)
        !!
        !! For the third phase...(solid)
        real(kind=8), ALLOCATABLE :: s_phiLS(:,:)                    ! interface variables
        real(kind=8), allocatable :: s_phiLSn(:,:)
        real(kind=8), ALLOCATABLE :: s_H(:,:),s_HOLD(:,:)            ! heavieside field
        !!
        real(kind=8), allocatable :: w_phiLS(:,:)
        real(kind=8), ALLOCATABLE :: w_H(:,:)                        ! heavieside field
        !!
        real(kind=8), ALLOCATABLE :: H(:,:),HOLD(:,:)                ! heavieside field
        real(kind=8), ALLOCATABLE :: u(:,:),v(:,:),&
                                     P(:,:),Pold(:,:),&
                                     uINT(:,:),vINT(:,:)             ! variables
        real(kind=8), ALLOCATABLE :: ustar(:,:)                      ! temporary u-comp velocity
        real(kind=8), ALLOCATABLE :: vstar(:,:)                      ! temporary v-comp velocity
        real(kind=8), ALLOCATABLE :: u_old(:,:)                      ! An array to hold a variable
        real(kind=8), ALLOCATABLE :: v_old(:,:)                      ! An array to hold a variable        
        real(kind=8), ALLOCATABLE :: phiLSstar(:,:)                  ! temporary level set field
        real(kind=8), ALLOCATABLE :: s_phiLSstar(:,:)                ! temporary level set field        
        real(kind=8), ALLOCATABLE :: PPERHS(:,:)                     ! Poission right-hand-side
        real(kind=8), ALLOCATABLE :: BDA(:),BDB(:)                   ! BC arrays for fishpac solver
        real(kind=8), ALLOCATABLE :: BDC(:),BDD(:)                   ! BC arrays for fishpac solver
        real(kind=8), dimension(10000) :: WORK                       ! work array for fishpac solution
        real(kind=8), ALLOCATABLE :: xwhole(:),ywhole(:)             ! whole node locations
        
        real(kind=8), ALLOCATABLE :: n_rbdelom(:,:)                  ! Rotation tensor
        real(kind=8), ALLOCATABLE :: n1_rbdelom(:,:)                 ! Rotation tensor
        real(kind=8), ALLOCATABLE :: exp_big_theta(:,:) 
        
        real(kind=8), ALLOCATABLE :: x_v(:),y_v(:),z_v(:)            ! square vairables
        real(kind=8), ALLOCATABLE :: nx_v(:),ny_v(:)
        real(kind=8), ALLOCATABLE :: pt(:),rot_mat(:,:),t(:,:)
        !
        ! For the internals
        real(kind=8), ALLOCATABLE :: int_data(:)                    ! Internal data - used for graphic
        !
        !
        ! EXTERNAL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: LSHEAVY
        REAL(kind=8), EXTERNAL :: LSCURVATURE
        REAL(kind=8), EXTERNAL :: EXPINT, ERF
        REAL(kind=8), EXTERNAL :: NORM
        !
        !
        !
        ! //////////////////////////// INPUTS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\        
        filename = 'inputfile_new.txt'                                   ! The name of the input file
        
        pi = 3.1415926535897932384626433832795

        ! OPENING DATA FILE
        OPEN (UNIT=2,FILE=filename,STATUS="OLD",IOSTAT=ios)
        
        ! READ THE INTPUT FILE
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        !
        ! ** GEOMETRY **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) Nx                                      ! Nx corresponds to the number of interior x whole nodes
        READ (UNIT=2,FMT=60) Ny                                      ! Ny corresponds to the number of interior y whole nodes
        READ (UNIT=2,FMT=70) Lx                                      ! total x length 
        READ (UNIT=2,FMT=70) Ly                                      ! total y length
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        ! ** SOLVERS **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) solvertype                              ! solver for pressure poisson equation
                                                                     ! (1=FFT,2=PCG,3=SOR)
                                                                     ! NOTE: FFT only works for single-phase flow and square grids
        READ (UNIT=2,FMT=60) itermax                                 ! maximum iterations for pressure convergence -> int
        READ (UNIT=2,FMT=70) tolerance                               ! tolerance to stop pressure iteration
        READ (UNIT=2,FMT=60) iterprint                               ! print pressure residual every iterprint interations -> int
        READ (UNIT=2,FMT=60) precond                                 ! preconditioner (PCG solver), 0=none, 1=diagonal -> int
        READ (UNIT=2,FMT=70) omega                                   ! SOR weight, omega=1 corresponds to GS
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        ! ** TIME **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) time_print                              ! print every time_print seconds  -> real(8)
        READ (UNIT=2,FMT=70) time_max                                ! max time -> real(8)
        READ (UNIT=2,FMT=70) stability                               ! coefficient to reduce time step such that code is stable.
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        ! ** INTERFACE **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=80) trackphase                              ! track the interface
        READ (UNIT=2,FMT=80) immersedbc                              ! use immersedbc
        READ (UNIT=2,FMT=80) rigidbody                               ! add rigid body
        READ (UNIT=2,FMT=80) wedge                                   ! use wedge
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        !** PROPERTIES ** 
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) rb_shape                                ! shape of the rigid body
                                                                     !   1 = circle
                                                                     !   2 = ellipse
        READ (UNIT=2,FMT=60) int_type                                ! type of internal component(s)
                                                                     !   0 = none
                                                                     !   1 = disc
                                                                     !   2 = 2D trans
                                                                     !   3 = 1D trans
        READ (UNIT=2,FMT=60) RB_lock                                 ! restrict motion of RB
                                                                     !   0 = No restraints
                                                                     !   1 = Vertical motion only
        READ (UNIT=2,FMT=70) rb1_mass                                ! mass of the rigid body
        READ (UNIT=2,FMT=70) buoy_J                                  ! rotational moment of inertia of the rigid body
        READ (UNIT=2,FMT=70) rho_one                                 ! density of first fluid, H=1, phiLS>0        (liquid)
        READ (UNIT=2,FMT=70) rho_two                                 ! density of second fluid, H=0, phiLS<0        (vapor)
        READ (UNIT=2,FMT=70) mu_one                                  ! viscosity of first fluid
        READ (UNIT=2,FMT=70) mu_two                                  ! viscosity of second fluid
        READ (UNIT=2,FMT=70) sigma                                   ! surface tension        
        READ (UNIT=2,FMT=70) gx                                      ! gravity constant (positive is down)
        READ (UNIT=2,FMT=70) gy                                      ! gravity constant (positive is down)
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        !** BOUNDARY TYPES **        
        ! velocity BC types - (1) dirclet wall, (2) neumman wall (also symmetry),   
        !                     (3) free surface, (4) pressure, (5) flow, (6) contact line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) BCnorth 
        READ (UNIT=2,FMT=60) BCsouth  
        READ (UNIT=2,FMT=60) BCeast
        READ (UNIT=2,FMT=60) BCwest
        !
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        !
        !** BOUNDARY CONDITIONS **
        ! velocity boundary conditions [applicable for BC types (1) and (4)]
        READ (UNIT=2,FMT=50) variables                               ! skip text line        
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) unorth 
        READ (UNIT=2,FMT=70) usouth 
        READ (UNIT=2,FMT=70) ueast
        READ (UNIT=2,FMT=70) uwest
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !        
        READ (UNIT=2,FMT=70) vnorth 
        READ (UNIT=2,FMT=70) vsouth 
        READ (UNIT=2,FMT=70) veast 
        READ (UNIT=2,FMT=70) vwest
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        ! pressure dirclet BC's
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) Pnorth 
        READ (UNIT=2,FMT=70) Psouth 
        READ (UNIT=2,FMT=70) Peast
        READ (UNIT=2,FMT=70) Pwest 
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        ! contact line boundary conditions 
        ! Ucl = contC*(angle - contA)**contM if angle(phi) > contA
        ! Ucl = -contC*(contR - angle)**contM if angle(phi) < contR
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) contC                                   ! constant        
        READ (UNIT=2,FMT=70) contA                                   ! static advancing contact angle
        READ (UNIT=2,FMT=70) contR                                   ! static receding contact angle
        READ (UNIT=2,FMT=70) contM                                   ! mobility exponent
        !        
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR FREE SURFACE**        
        ! phiLS = a*sqrt( (x-b)^2 +(y-c)^2 ) + d + ex + fy
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) a_phiLS
        READ (UNIT=2,FMT=70) b_phiLS
        READ (UNIT=2,FMT=70) c_phiLS
        READ (UNIT=2,FMT=70) d_phiLS
        READ (UNIT=2,FMT=70) e_phiLS
        READ (UNIT=2,FMT=70) f_phiLS
        READ (UNIT=2,FMT=70) v_phiLS
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR CIRCULAR RIGID BODY**        
        ! s_phiLS = a*sqrt( (x-b)^2 +(y-c)^2 ) + d
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) a_rbLS
        READ (UNIT=2,FMT=70) b_rbLS
        READ (UNIT=2,FMT=70) c_rbLS
        READ (UNIT=2,FMT=70) d_rbLS
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR ELLIPTICAL RIGID BODY**        
        ! s_phiLS = 
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) aellipse
        READ (UNIT=2,FMT=70) bellipse
        READ (UNIT=2,FMT=70) xo
        READ (UNIT=2,FMT=70) yo
        READ (UNIT=2,FMT=70) tho        
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR SQUARE RIGID BODY**        
        ! s_phiLS = 
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) c_xo
        READ (UNIT=2,FMT=70) c_yo
        READ (UNIT=2,FMT=70) x1
        READ (UNIT=2,FMT=70) y1
        READ (UNIT=2,FMT=70) x2
        READ (UNIT=2,FMT=70) y2
        READ (UNIT=2,FMT=70) x3
        READ (UNIT=2,FMT=70) y3
        READ (UNIT=2,FMT=70) x4
        READ (UNIT=2,FMT=70) y4
        READ (UNIT=2,FMT=70) c_tho
        READ (UNIT=2,FMT=70) cgx
        READ (UNIT=2,FMT=70) cgy
        !        
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** TORSIONAL DISC PARAMETERS (int_type = 1) ***        
        READ (UNIT=2,FMT=50) variables   
        READ (UNIT=2,FMT=70) disk_mass                                ! mass of the torsional disc
        READ (UNIT=2,FMT=70) small_J                                 ! rotational moment of inertia of the disc
        READ (UNIT=2,FMT=70) k12                                     ! torsional spring constant
        READ (UNIT=2,FMT=70) c12                                     ! torsional damping constant
        !        
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** SLIDER PARAMETERS (int_type = 2 OR 3) ***
        ! int_type = 2 IGNORES slider_J
        ! int_type = 3 IGNORES kx12 and cx12       
        READ (UNIT=2,FMT=50) variables   
        READ (UNIT=2,FMT=70) slider_mass                             ! mass of the slider
        READ (UNIT=2,FMT=70) slider_J                                ! rotational moment of inertia of the slider
        READ (UNIT=2,FMT=70) kx12                                    ! slider spring constant (x-direction)
        READ (UNIT=2,FMT=70) cx12                                    ! slider damping constant (x-direction)
        READ (UNIT=2,FMT=70) ky12                                    ! slider spring constant (y-direction)
        READ (UNIT=2,FMT=70) cy12                                    ! slider damping constant (y-direction)
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        CLOSE (UNIT=2)                                               ! close the input file
        !
        !
        !
        ! ///////////////////// INITIALIZE VARIABLES \\\\\\\\\\\\\\\\\\\\\\
        ! allocating the size of each arrays        
        ALLOCATE (x(Nx+1))
        ALLOCATE (y(Ny+1))
        ALLOCATE (xwhole(Nx+2))
        ALLOCATE (ywhole(Ny+2))
        !        
        ALLOCATE (phi(Nx+2,Ny+2))
        ALLOCATE (phiLS(Nx+2,Ny+2))
        ALLOCATE (phiLSn(Nx+2,Ny+2))
        AllOCATE (H(Nx+2,Ny+2))
        AllOCATE (HOLD(Nx+2,Ny+2))
        
        ALLOCATE (s_phiLS(Nx+2,Ny+2))
        ALLOCATE (s_phiLSn(Nx+2,Ny+2))
        AllOCATE (s_H(Nx+2,Ny+2))
        AllOCATE (s_HOLD(Nx+2,Ny+2))
        !
        ALLOCATE (w_phiLS(Nx+2,Ny+2))
        AllOCATE (w_H(Nx+2,Ny+2))
        !
        AllOCATE (u(Nx+1,Ny+2))
        AllOCATE (v(Nx+2,Ny+1))
        AllOCATE (u_old(Nx+1,Ny+2))
        AllOCATE (v_old(Nx+2,Ny+1))
        AllOCATE (P(Nx+2,Ny+2))
        ALLOCATE (Pold(Nx+2,Ny+2))
        AllOCATE (uINT(Nx+1,Ny+2))
        AllOCATE (vINT(Nx+2,Ny+1))
        
        !For the square
        AllOCATE (x_v(3))
        AllOCATE (y_v(3))
        AllOCATE (nx_v(3))
        AllOCATE (ny_v(3))
        AllOCATE (z_v(3))
        AllOCATE (pt(2))
        ALLOCATE (rot_mat(2,2))
        AllOCATE (t(2,3))
        !
        AllOCATE (ustar(Nx+1,Ny+2))
        AllOCATE (vstar(Nx+2,Ny+1)) 
        AllOCATE (phiLSstar(Nx+2,Ny+2))
        AllOCATE (s_phiLSstar(Nx+2,Ny+2))
        
        ALLOCATE (n_rbdelom(2,2))
        ALLOCATE (n1_rbdelom(2,2))
        ALLOCATE (exp_big_theta(2,2))
        !        
        ! PPE RHS
        AllOCATE (PPERHS(Nx,Ny))
        !
        ! BC arrays for fishpac solver
        AllOCATE (BDA(Ny))
        AllOCATE (BDB(Ny))   
        AllOCATE (BDC(Nx))
        AllOCATE (BDD(Nx))
        !
        ! Internal Data
        ALLOCATE (int_data(Nx+1))
        !
        !
        ! time and counter initialization
        time = 0.0                                                   ! time is zero at start
        file_count = 0                                               ! print file counter 
        time_count = 0.0                                             ! counter for when to print
        dumpcount = 0                                                ! dump file counter
        
        !
        ! RK4 coefficients
        Coef(1) = 1.0/4.0
        Coef(2) = 1.0/3.0
        Coef(3) = 1.0/2.0
        Coef(4) = 1.0

        ! Initialize rigid body variables
        
        n_rb_dx = 0.0
        n_rb_dy = 0.0
        n1_rb_dx = 0.0
        n1_rb_dy = 0.0
        n_w_dy = 0.0
        n1_w_dy = 0.0
        n12_w_vy = 0.0
        n1_w_vy = 0.0
        
        n_rb_vx = 0.0
        n_rb_vy = 0.0
        n12_rb_vx = 0.0
        n12_rb_vy = 0.0
        n1_rb_vx = 0.0
        n1_rb_vy = 0.0

        n_rb_ax = 0.0
        n_rb_ay = 0.0
        n1_rb_ax = 0.0
        n1_rb_ay = 0.0
        
        ! Rotation variables
        
        c1 = 0.0
        c2 = 0.0
        
        rx = 0.0
        ry = 0.0
        
        n_rbdelom(1,1) = 1.0
        n_rbdelom(1,2) = 0.0
        n_rbdelom(2,1) = 0.0
        n_rbdelom(2,2) = 1.0
        
        n1_rbdelom(1,1) = 1.0
        n1_rbdelom(1,2) = 0.0
        n1_rbdelom(2,1) = 0.0
        n1_rbdelom(2,2) = 1.0
        
        ! The exponential map
        exp_big_theta(1,1) = 0.0
        exp_big_theta(1,2) = 0.0
        exp_big_theta(2,1) = 0.0
        exp_big_theta(2,2) = 0.0
        
        big_theta = 0.0
        
        n_rb_vom = 0.0
        n12_rb_vom = 0.0
        n1_rb_vom = 0.0
        
        n_rb_aom = 0.0
        n1_rb_aom = 0.0

        th_change = 0.0
        th_domain = 0.0
        th_old = 0.0
        th_counter = 0.0

           ! Internal Parameters
 
           int_data = 0.0                                ! Internal data for postprocessing
           int_data(1) = int_type

           ! RB/Internals components
           
           dtinternal = 10.0

           ncgx = cgx
           ncgy = cgy
           netcgx = cgx
           netcgy = cgy

           if (int_type == 0) then
             
             rb_mass = rb1_mass
             big_J = buoy_J

           elseif (int_type == 1) then
             
             rb_mass = rb1_mass + disk_mass
             big_J = buoy_J
             
             if ((k12 /= 0.0) .and. (c12 /= 0.0)) then
             dtinternal = min( sqrt(min(big_J,small_J)/k12), sqrt(min(big_J,small_J)/c12) )
             endif

             aom_counter = 0
             
             if (k12 /= 0.0) then
                aom_limit = 400/k12
             else
                aom_limit = 500
             endif
           
           elseif (int_type == 2) then

             rb_mass = rb1_mass
             big_J = buoy_J
             old_big_J = big_J
             scgx = cgx                                
             scgy = cgy
             max_w = x1 - cgx
             max_e = x2 - cgx
             max_n = y1 - cgy
             max_s = y3 - cgy
             
             if ((kx12 /= 0.0) .and. (ky12 /= 0.0)) then
             dtinternal = 0.5 * min( sqrt(min(rb1_mass,slider_mass)/kx12), &
                                     sqrt(min(rb1_mass,slider_mass)/ky12) )
             endif

           elseif (int_type == 3) then

             rb_mass = rb1_mass
             big_J = buoy_J + slider_J                ! Initially
             old_big_J = big_J
             scgx = cgx                                
             scgy = cgy
             max_w = x1 - cgx
             max_e = x2 - cgx
             max_n = y1 - cgy
             max_s = y3 - cgy


           end if

           ! Torsional Disc
           small_theta_n = 0.0
           small_theta_n1 = 0.0
 
           n_rb2_vom = 0.0
           n12_rb2_vom = 0.0
           n1_rb2_vom = 0.0
        
           n_rb2_aom = 0.0
           n1_rb2_aom = 0.0
           delta_theta = 0.0
           delta_vom = 0.0
           int_torque = 0.0
           int_work = 0.0

           ! Slider
           n_rb2_dx = 0.0                            ! Displacement        
           n_rb2_dy = 0.0
           n1_rb2_dx = 0.0
           n1_rb2_dy = 0.0
        
           n_rb2_vx = 0.0                            ! Velocities
           n_rb2_vy = 0.0                          
           n12_rb2_vx = 0.0
           n12_rb2_vy = 0.0
           n1_rb2_vx = 0.0
           n1_rb2_vy = 0.0
        
           n_rb2_ax = 0.0                            ! Acceleration  
           n_rb2_ay = 0.0
           n1_rb2_ax = 0.0
           n1_rb2_ay = 0.0
       
           delta_dpx = 0.0                           ! Differences in displacement and velocity
           delta_dpy = 0.0                     
           delta_vpx = 0.0
           delta_vpy = 0.0

           diff_J = 0.0
           int_forcex = 0.0                          ! Internal forces                  
           int_forcey = 0.0
           int_forcepx = 0.0
           int_forcepy = 0.0
           
           ext_px = 0.0                              ! External forces
           ext_py = 0.0

        ! Initialize intergral variables two different directions
        n_xint_conv  = 0.0
        n_xint_grav  = 0.0
        n_xint_inert = 0.0
        
        n_yint_conv  = 0.0
        n_yint_grav  = 0.0
        n_yint_inert = 0.0
        
        ! Initialize intergral variables two different directions
        n1_xint_conv  = 0.0
        n1_xint_grav  = 0.0
        n1_xint_inert = 0.0
        
        n1_yint_conv  = 0.0
        n1_yint_grav  = 0.0
        n1_yint_inert = 0.0

        ! Initialize integral variables for wedge (in one direction)
        
        n_w_conv  = 0.0
        n_w_grav  = 0.0
        n_w_inert = 0.0
        
        n1_w_conv  = 0.0
        n1_w_grav  = 0.0
        n1_w_inert = 0.0

        wedge_force = 0.0
        wedge_work = 0.0

        ! Initialize rotation integral variables
        r_int_conv  = 0.0
        r_int_grav  = 0.0
        r_int_inert = 0.0
        
        r1_int_conv  = 0.0
        r1_int_grav  = 0.0
        r1_int_inert = 0.0
        
        ustar_int = 0.0
        vstar_int = 0.0
        
        xn = 0.0
        yn = 0.0
        th = 0.0
        xc = 0.0
        yc = 0.0

        ! 1 = neumman, 2 = dirclet, 3 = extrapolation
        BCphisouth = 1        ! 3
        BCphinorth = 1        ! 3
        BCphiwest = 1        ! 1
        BCphieast = 1        ! 1

        ! ////////////////////// COORDINATE GENERATION \\\\\\\\\\\\\\\\\\\\\\
        deltaX = Lx/(real(Nx))
        deltaY = Ly/(real(Ny))
        x(1) = 0.0
        y(1) = 0.0
        x(Nx+1) = Lx
        y(Ny+1) = Ly

        do i=2,Nx
          x(i) = deltaX*real(i-1)
        enddo

        do j=2,Ny
          y(j) = deltaY*real(j-1)
        enddo

        ! whole node locations
        do i=1,Nx+1
          xwhole(i) = x(i)-0.5*deltaX
        enddo

        xwhole(Nx+2)=x(Nx+1)+0.5*deltaX

        do j=1,Ny+1
          ywhole(j) = y(j)-0.5*deltaY
        enddo

        ywhole(Ny+2)=y(Ny+1)+0.5*deltaY

        
        ! ///////////////// INITIAL CONDITIONS \\\\\\\\\\\\\\\\\\\

        ! ---- IMMERSED BOUNDARY ----        
        
        s_phiLS = 0.0                                                ! Rigid Body
        phiLS = 0.0                                                  ! Free surface
        phi = 0.0                                                    ! Immersed Boundary
 
         
        if (immersedbc .EQV. .TRUE.) then
          do i=1,Nx+2
            do j=1,Ny+2
              ! phi>0 is solid, phi<0 is fluid
              phi(i,j) = -sqrt( (xwhole(i)-0.5)**2.0 +&
                                (ywhole(j)-0.5)**2.0 ) + 0.2 
            enddo
          enddo
        endif
                
        ! ---- PRESSURE        ----        
        P = 0.0                                                      ! initializing array to zero

        ! ---- INITIAL LEVEL SET FIELDS  ----
        ! The free surface
        if (trackphase .EQV. .TRUE.) then
          do i=1,Nx+2
            do j=1,Ny+2

              phiLS(i,j) = a_phiLS*sqrt( (xwhole(i)-b_phiLS)**2.0+&
                                         (ywhole(j)-c_phiLS)**2.0 )+&
                                         d_phiLS+&
                                         e_phiLS*xwhole(i)+&
                                         f_phiLS*ywhole(j)
                                                        
              !if (ywhole(j) < (v_phiLS + 0.05)) then
              !  phiLS(i,j) = -ywhole(j) + v_phiLS
              !endif

            enddo
          enddo
        endif

        ! The wedge
        if (wedge .EQV. .TRUE.) then
          do i=1,Nx+2
            do j=1,Ny+2        
        
              if (xwhole(i) < 0.1) then
                  w_phiLS(i,j) = -(-ywhole(j) + 0.65)
              else
                  w_phiLS(i,j) = -(-ywhole(j) + xwhole(i) + 0.55)
              endif
                        
            enddo
          enddo
        endif

        ! The Rigid Body
        ! Circle
        if ((rigidbody .EQV. .TRUE.) .AND. (rb_shape == 1)) then
          do i=1,Nx+2
            do j=1,Ny+2
              
              th = 0.0
              th_old = th
              xc = b_rbLS
              yc = c_rbLS

              s_phiLS(i,j) = a_rbLS*sqrt( (xwhole(i)-b_rbLS)**2.0+&
                                          (ywhole(j)-c_rbLS)**2.0 )+&
                                          d_rbLS
                                                        
            enddo
          enddo
                
        ! Ellipse
        elseif ((rigidbody .EQV. .TRUE.) .AND. (rb_shape == 2)) then
        
          xn = xo
          yn = yo
          xc = xn
          yc = yn
          th = pi/tho
          th_old = th
          print*, "tho"
          print*, tho
          print*, "th"
          print*, th
          do i=1,Nx+2
            do j=1,Ny+2
                                
              s_phiLS(i,j) = -100*&
                             ( bellipse**2*( xwhole(i)*cos(th)+&
                                             ywhole(j)*sin(th)-&
                                             (xo*cos(th) + yo*sin(th)))**2 + &
                               aellipse**2*(-xwhole(i)*sin(th) + ywhole(j)*cos(th) - &
                                           (-xo*sin(th) + yo*cos(th)))**2 - &
                               aellipse**2*bellipse**2)

            enddo
          enddo
        
        ! Square
        elseif ((rigidbody .EQV. .TRUE.) .AND. (rb_shape == 3)) then
        
          xc = c_xo  ! Geometric Center
          yc = c_yo
                
          !th = pi/c_tho
          th = c_tho
          th_old = th
          cost = cos(c_tho)
          sint = sin(c_tho)
          rot_mat(1,1) = cost
          rot_mat(1,2) = -sint
          rot_mat(2,1) = sint
          rot_mat(2,2) = cost
                
          !print*, "theta"
          !print*, c_tho
          !print*, "xc"
          !print*, xc
          !print*, "yc"
          !print*, yc
          !print*, "x1"
          !print*, x1
          !print*, "y1"
          !print*, y1
          !print*, "x2"
          !print*, x2
          !print*, "y2"
          !print*, y2
          !print*, "x3"
          !print*, x3
          !print*, "y3"
          !print*, y3
          !print*, "x4"
          !print*, x4
          !print*, "y4"
          !print*, y4
                        
          ! x1-4,y-4 are ordered counterclockwise
                
          do r=1,4
                
            if (r == 1) then
                        
              ! Region 1
              x_v(1) = xc;
              x_v(2) = x1;
              x_v(3) = x2;
                                
              y_v(1) = yc;
              y_v(2) = y1;
              y_v(3) = y2;
                                
            elseif (r == 2) then
                        
              ! Region 2
              x_v(1) = xc;
              x_v(2) = x2;
              x_v(3) = x3;
                                
              y_v(1) = yc;
              y_v(2) = y2;
              y_v(3) = y3;
                                
            elseif (r == 3) then
                        
              ! Region 3
              x_v(1) = xc;
              x_v(2) = x3;
              x_v(3) = x4;
                                
              y_v(1) = yc;
              y_v(2) = y3;
              y_v(3) = y4;
                                
            elseif (r == 4) then
                        
              ! Region 4
              x_v(1) = xc;
              x_v(2) = x4;
              x_v(3) = x1;
                                
              y_v(1) = yc;
              y_v(2) = y4;
              y_v(3) = y1;
                                
            endif

            nx_v(1) = rot_mat(1,1)*(xc-cgx) + rot_mat(1,2)*(yc-cgy)
            ny_v(1) = rot_mat(2,1)*(xc-cgx) + rot_mat(2,2)*(yc-cgy)                        
            nx_v(2) = rot_mat(1,1)*(x_v(2)-cgx) + rot_mat(1,2)*(y_v(2)-cgy)
            ny_v(2) = rot_mat(2,1)*(x_v(2)-cgx) + rot_mat(2,2)*(y_v(2)-cgy)
            nx_v(3) = rot_mat(1,1)*(x_v(3)-cgx) + rot_mat(1,2)*(y_v(3)-cgy)
            ny_v(3) = rot_mat(2,1)*(x_v(3)-cgx) + rot_mat(2,2)*(y_v(3)-cgy)

            x_v(1) = nx_v(1) + cgx
            y_v(1) = ny_v(1) + cgy                        
            x_v(2) = nx_v(2) + cgx
            y_v(2) = ny_v(2) + cgy
            x_v(3) = nx_v(3) + cgx
            y_v(3) = ny_v(3) + cgy
                                
            z_v(1) = 1.0;
            z_v(2) = 0.0;
            z_v(3) = 0.0;
                                
            Det = x_v(1)*(y_v(2)*z_v(3)-y_v(3)*z_v(2)) - &
                      y_v(1)*(x_v(2)*z_v(3)-x_v(3)*z_v(2)) + &
                      z_v(1)*(x_v(2)*y_v(3)-x_v(3)*y_v(2))
                                          
            a_det = y_v(2)*z_v(3)-y_v(3)*z_v(2) - &
                            y_v(1)*(z_v(3)-z_v(2)) + &
                            z_v(1)*(y_v(3)-y_v(2))
                                          
            b_det = x_v(1)*(z_v(3)-z_v(2)) - &
                            x_v(2)*z_v(3)-x_v(3)*z_v(2) + &
                            z_v(1)*(x_v(2)-x_v(3))
                                          
            c_det = x_v(1)*(y_v(2)-y_v(3)) - &
                            y_v(1)*(x_v(2)-x_v(3)) + &
                            x_v(2)*y_v(3)-x_v(3)*y_v(2)
                                                
            a_rbLS = a_det/Det
            b_rbLS = b_det/Det
            c_rbLS = c_det/Det
                                
            do i=1,Nx+2
              do j=1,Ny+2
                        
                ! Check on the location of the grid point.
                ! Is it in the correct region?
                                
                pt(1) = xwhole(i)
                pt(2) = ywhole(j)
                                
                t(1,1) = x_v(1);
                t(2,1) = y_v(1);
                t(1,2) = x_v(2);
                t(2,2) = y_v(2);
                t(1,3) = x_v(3);
                t(2,3) = y_v(3);
                                
                do m = 1,3,2
                                
                  k = mod(m,3) + 1
                                        
                  temp = ( pt(1) - t(1,m) ) * ( t(2,k) - t(2,m) ) - &
                         ( pt(2) - t(2,m) ) * ( t(1,k) - t(1,m) )
                                                   
                  if (0.0 < temp) then
                    inside = .FALSE.
                    exit
                  endif
                  
                  inside = .TRUE.
                enddo
                                
                if (inside .EQV. .TRUE.) then
                                
                  s_phiLS(i,j) = (-a_rbLS*xwhole(i)-&
                                   b_rbLS*ywhole(j) + 1)/c_rbLS
                                                                   
                endif
              enddo                
            enddo
          enddo
        endif
                                                     

        ! initialize the new LS field to a signed distance function
        CALL LSREIN(Nx,Ny,deltaX,deltaY,phiLS,H,&
                    BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)

        !CALL LSREIN(Nx,Ny,deltaX,deltaY,s_phiLS,s_H,&
        !           BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)                                
                                

        ! ---- HEAVIESIDE FUNCTION ----
        if (trackphase .EQV. .TRUE.) then
          do i=1,Nx+2
            do j=1,Ny+2
              ! Heavieside function
              H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,phiLS)
            enddo
          enddo
        else
          H = 1.0                ! outside interface=1, inside interface=0
        endif        

        ! ---- HEAVIESIDE FUNCTION SOLID ----
        if (rigidbody .EQV. .TRUE.) then
          do i=1,Nx+2
            do j=1,Ny+2
              ! Heavieside function
              s_H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,s_phiLS)
              !s_H(i,j) = 0.0
            enddo
          enddo
        else
          s_H = 1.0                ! outside interface=1, inside interface=0
        endif        

        ! The wedge
        if (wedge .EQV. .TRUE.) then
          do i=1,Nx+2
            do j=1,Ny+2
              ! Heavieside function
              w_H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,w_phiLS)
              !s_H(i,j) = 0.0
            enddo
          enddo
        endif
        
        ! ---- FLUID VELOCITY ----
        u = 0.0                
        v = 0.0
        u_old = 0.0
        v_old = 0.0
        deltaT = min(1.0/( max(mu_one/rho_one,mu_two/rho_two)*&
                         (2.0/(deltaX**2.0) + 2.0/(deltaY**2.0)) ),dtinternal)

        ! Calculate initial work done by wedge
        if (wedge .EQV. .TRUE.) then 
          
          do i=2,Nx+1
            do j=2,Ny     
              
              ! density at the half node
              rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                      0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))

              if (0.5*(w_H(i,j)+w_H(i,j+1)) > 0) then                                      
                
                vstar_int = 0.0
                
                n_w_inert = n1_w_inert + &
                                0.5*(w_H(i,j)+w_H(i,j+1))*&
                                rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                deltaX*deltaY
                                                                                        
                n_w_conv = n1_w_conv + &
                               0.5*(w_H(i,j)+w_H(i,j+1))*&
                               rho_O*vstar_int*deltaX*deltaY
                                
                n_w_grav = n1_w_grav + 0.5*&
                               (w_H(i,j)+w_H(i,j+1))*rho_O*deltaX*deltaY*gy

              endif
              
              wedge_force = n_w_inert + n_w_conv + n_w_grav
              wedge_work = wedge_work - (wedge_force * n1_w_vy * deltaT)

            enddo
          enddo
        endif


        ! --- INITIAL SET OF RIGID BODY FORCES

        total_energy = 0.0
        lin_momx = 0.0
        lin_momy = 0.0
        
        ! Fluid Energy Contribution
        if (rigidbody .EQV. .FALSE.) then
          do i=2,Nx+1
            do j=2,Ny+1
  
              if (i < Nx+1) then
                rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))+&
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                 
                total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*deltaX*deltaY  ! kinetic (in x)
                lin_momx = lin_momx + rho_O*u(i,j)*deltaX*deltaY
              endif
              
              if (j < Ny+1) then
                rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))                          

                total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*deltaX*deltaY ! kinetic (in y)
                lin_momy = lin_momy + rho_O*v(i,j)*deltaX*deltaY
              endif

              if ((i < Nx+1) .AND. (j < Ny+1)) then 
                rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                total_energy = total_energy + &
                               ( -P(i,j) + rho_O*gy*(j-1)*deltaY ) * &
                               deltaX*deltaY
              endif

            enddo
          enddo
        endif

        if (rigidbody .EQV. .TRUE.) then
          ! Collect some integrals
                
          do i=2,Nx
            do j=2,Ny+1
                        
              ! density at the half node
              rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))+&
                      0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                 
              total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*deltaX*deltaY  ! kinetic (in x)
              lin_momx = lin_momx + rho_O*u(i,j)*deltaX*deltaY

              if (0.5*(s_H(i,j)+s_H(i+1,j)) > 0) then
                                                
                ! dummy variable for the convective opperator
                ustar_int= 0.0
                                                
                n1_xint_inert = n1_xint_inert+&
                                0.5*(s_H(i,j)+s_H(i+1,j))*&
                                rho_O/DeltaT*(u(i,j) - u_old(i,j))*&
                                deltaX*deltaY
                                                
                n1_xint_conv = n1_xint_conv+&
                               0.5*(s_H(i,j)+s_H(i+1,j))*&
                               rho_O*ustar_int*deltaX*deltaY
                                
                n1_xint_grav = n1_xint_grav+&
                               0.5*(s_H(i,j)+s_H(i+1,j))*&
                               rho_O*deltaX*deltaY*gx
                                                
                ! Rotations
                ! rs is the s-distance from the location at i,j to the centroid of the object.
                ry = ywhole(j) - netcgy
                                                
                r1_int_inert = r1_int_inert - &
                               ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                               rho_O/deltaT*(u(i,j) - u_old(i,j))*&
                               deltaX*deltaY
                                                
                r1_int_conv = r1_int_conv - &
                              ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                              rho_O*ustar_int*deltaX*deltaY
                                
                r1_int_grav = r1_int_grav - &
                              ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                              rho_O*deltaX*deltaY*gx
                                
              endif

            enddo
          enddo


          do i=2,Nx+1
            do j=2,Ny
                        
              ! density at the half node
              rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                      0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))                          

              total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*deltaX*deltaY ! kinetic (in y)
              lin_momy = lin_momy + rho_O*v(i,j)*deltaX*deltaY

              if (0.5*(s_H(i,j)+s_H(i,j+1)) > 0) then
                                        
                ! dummy variable for the convective operator
                vstar_int = 0.0
                                                        
                n_yint_inert = n1_yint_inert +&
                               0.5*(s_H(i,j)+s_H(i,j+1))*&
                               rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                               deltaX*deltaY
                                                                                        
                n_yint_conv = n1_yint_conv + &
                              0.5*(s_H(i,j)+s_H(i,j+1))*&
                              rho_O*vstar_int*deltaX*deltaY
                                
                n_yint_grav = n1_yint_grav + &
                              0.5*(s_H(i,j)+s_H(i,j+1))*&
                              rho_O*deltaX*deltaY*gy
                                                
                ! Rotations
                ! rs is the s-distance from the location at i,j to the centroid of the object.
                rx = xwhole(i) - netcgx
                                                
                r_int_inert = r1_int_inert + &
                              rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                              rho_O/deltaT*(v(i,j) - v_old(i,j))*&
                              deltaX*deltaY
                                                                                        
                r_int_conv = r1_int_conv + &
                             rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                             rho_O*vstar_int*deltaX*deltaY
                                
                r_int_grav = r1_int_grav + &
                             rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                             rho_O*deltaX*deltaY*gy
              endif                                  
                
            enddo
          enddo

          ! Pressure and gravitational contributions to total energy
          do i=2,Nx
            do j=2,Ny     
              
              ! density at the whole node
              rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

              total_energy = total_energy + &
                             ( -P(i,j) + rho_O*gy*(j-1)*deltaY ) * &
                             deltaX*deltaY

            enddo
          enddo

          ! Initial internal forces          
          if (int_type == 1) then                               ! Torque exerted between disc and RB
              
              small_theta_n = th
              delta_theta = small_theta_n - th
              delta_vom = n_rb2_vom - n_rb_vom
              int_torque = k12*delta_theta + c12*delta_vom 
              int_data(2) = cgx                                 ! Save data for post-processing
              int_data(3) = cgy
              int_data(4) = small_theta_n
              
              ! Acceleration
              n_rb2_aom = -int_torque/small_J

              ! Energy
              int_work = int_work + c12*(delta_vom**2)*deltaT
              total_energy = total_energy + &
                             0.5*small_J*(n_rb2_vom**2) + &
                             0.5*k12*(delta_theta**2)
     
          elseif (int_type == 2) then

              ! Coordinate Transformation
              delta_vpx =  cos(th)*(n_rb2_vx - n_rb_vx) + sin(th)*(n_rb2_vy - n_rb_vy)
              delta_vpy = -sin(th)*(n_rb2_vx - n_rb_vx) + cos(th)*(n_rb2_vy - n_rb_vy)

              delta_dpx =  cos(th)*(n_rb2_dx - n_rb_dx) + sin(th)*(n_rb2_dy - n_rb_dy)
              delta_dpy = -sin(th)*(n_rb2_dx - n_rb_dx) + cos(th)*(n_rb2_dy - n_rb_dy)

              if ((delta_dpx >= 0.8*max_e) .or. (delta_dpx <= 0.8*max_w)) then
                      kx_adj = 100*kx12
              else
                      kx_adj = (.64*kx12*max_e*max_w) / ((delta_dpx - 0.8*max_e)*(delta_dpx - 0.8*max_w))
              endif

              if ((delta_dpy >= 0.8*max_n) .or. (delta_dpy <= 0.8*max_s)) then
                      ky_adj = 100*ky12
              else
                      ky_adj = (.64*ky12*max_n*max_s) / ((delta_dpy - 0.8*max_n)*(delta_dpy - 0.8*max_s))

              endif

              int_forcepx = kx_adj*delta_dpx + cx12*delta_vpx
              int_forcepy = ky_adj*delta_dpy + cy12*delta_vpy

              ! Inverse Transformation
              int_forcex = cos(th)*int_forcepx - sin(th)*int_forcepy
              int_forcey = sin(th)*int_forcepx + cos(th)*int_forcepy
              
              ! Center of Mass update
              scgx = cgx + n_rb2_dx
              scgy = cgy + n_rb2_dy
              int_data(2) = cgx
              int_data(3) = cgy
              int_data(4) = scgx
              int_data(5) = scgy

              ! Acceleration
              n_rb2_ax = - int_forcex/slider_mass
              n_rb2_ay = - int_forcey/slider_mass - gy

              ! Energy
              int_work = int_work + deltaT * (cx12*(delta_vpx**2) + cy12*(delta_vpy**2))
              total_energy = total_energy + &
                             0.5*slider_mass*((n_rb2_vx**2) + (n_rb2_vy**2)) + &
                             0.5*kx12*(delta_dpx**2) + &
                             0.5*ky12*(delta_dpy**2)
              lin_momx = lin_momx + slider_mass*n_rb2_vx
              lin_momy = lin_momy + slider_mass*n_rb2_vy

          elseif (int_type == 3) then

              ! Coordinate Transformation

              ! delta_vpy = -sin(th)*(n_rb2_vx - n_rb_vx) + cos(th)*(n_rb2_vy - n_rb_vy)           
              ! delta_dpy = -sin(th)*(n_rb2_dx - n_rb_dx) + cos(th)*(n_rb2_dy - n_rb_dy)

              if ((delta_dpy >= max_n) .or. (delta_dpy <= max_s)) then
                      ky_adj = 10000
              else
                      ky_adj = (ky12*max_n*max_s) / ((delta_dpy - max_n)*(delta_dpy - max_s))
              endif

              int_forcepy = ky_adj*delta_dpy + cy12*delta_vpy + cos(th)*slider_mass*gy
              
              ! External Forces
              ext_px =  cos(th)*(n_xint_inert + n_xint_conv + n_xint_grav) + &
                        sin(th)*(n_yint_inert + n_yint_conv + n_yint_grav + sin(th)*slider_mass*gy)
              ext_py = -sin(th)*(n_xint_inert + n_xint_conv + n_xint_grav) + &
                        cos(th)*(n_yint_inert + n_yint_conv + n_yint_grav + sin(th)*slider_mass*gy)

              ! Inverse Transformation
              int_forcex = - sin(th)*int_forcepy
              int_forcey =   cos(th)*int_forcepy
              
              ! Center of Mass / Moment of Inertia update
              n_rb2_dx = - sin(th)*delta_dpy
              n_rb2_dy =   cos(th)*delta_dpy
              scgx = cgx + n_rb2_dx
              scgy = cgy + n_rb2_dy
              netcgx = (rb1_mass*cgx + slider_mass*scgx)/(rb1_mass + slider_mass)
              netcgy = (rb1_mass*cgy + slider_mass*scgy)/(rb1_mass + slider_mass)
              big_J = buoy_J + rb1_mass*((netcgx - cgx)**2 + (netcgy - cgy)**2) + &
                      slider_J + slider_mass*((netcgx - scgx)**2 + (netcgy - scgy)**2)
              diff_J = (big_J - old_big_J)/deltaT
              int_data(2) = cgx
              int_data(3) = cgy
              int_data(4) = scgx
              int_data(5) = scgy

              ! Buoy Acceleration
              n_rb_ax = cos(th)*(ext_px/(rb1_mass+slider_mass)) - sin(th)*(ext_py/rb1_mass) + &
                        int_forcex/rb1_mass - gx
              n_rb_ay = sin(th)*(ext_px/(rb1_mass+slider_mass)) + cos(th)*(ext_py/rb1_mass) + &
                        int_forcey/rb1_mass - gy
              n_rb_aom = (r_int_inert + r_int_conv + r_int_grav + int_torque - &
                         n_rb_vom*diff_J)/big_J

              ! Internal Acceleration
              ! n_rb2_ax = - int_forcex/slider_mass
              n_rb2_ay = - int_forcepy/slider_mass

              ! Energy
              int_work = int_work + deltaT * (cx12*(delta_vpx**2) + cy12*(delta_vpy**2))
              total_energy = total_energy + &
                             0.5*slider_mass*((n_rb2_vx**2) + (n_rb2_vy**2)) + &
                             0.5*kx12*(delta_dpx**2) + &
                             0.5*ky12*(delta_dpy**2)

          endif
      
          ! First set of rigid body accelerations
          if (int_type /= 3) then
          n_rb_ax = (n_xint_inert + n_xint_conv + n_xint_grav + int_forcex)/rb_mass - gx        
          n_rb_ay = (n_yint_inert + n_yint_conv + n_yint_grav + int_forcey)/rb_mass - gy        
          n_rb_aom = (r_int_inert + r_int_conv + r_int_grav + int_torque - &
                      n_rb_vom*diff_J)/big_J
          endif


          ! Restrict motion
          if (RB_lock == 1) then
             n_rb_ax = 0.0
             n_rb_aom = 0.0
             n1_rb_ax = 0.0
             n1_rb_aom = 0.0
          endif

          ! Energy Update
          total_energy = total_energy + &
                         rb_mass*gy*cgy + &
                         0.5*rb_mass*((n_rb_vx**2)+(n_rb_vy**2)) + &
                         0.5*big_J*(n_rb_vom**2)
          lin_momx = lin_momx + rb_mass*n_rb_vx
          lin_momy = lin_momy + rb_mass*n_rb_vy
          

        endif

        ! ---- SAVE INITIAL DATA ----
        call SAVEDATA(Nx,Ny,Lx,Ly,time,x,y,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v,file_count,int_data)
        
        ! //////////////////////// SOLVING EQUATIONS \\\\\\\\\\\\\\\\\\\\\\\
        do  ! time loop
        
        !Calculate net momentum
        net_mom = sqrt((lin_momx**2) + (lin_momy**2))


          ! CREATING DATA FILE
          OPEN (UNIT=8,FILE='centers.dat', POSITION ='append', IOSTAT=ios)
          WRITE (UNIT=8,FMT=90) n_rb_dx, n_rb_dy, n_rb_vx, n_rb_vy,&
                                n_yint_grav, n_yint_conv, n_yint_inert
          CLOSE (UNIT=8)
                
          OPEN (UNIT=7,FILE='rotations.dat', POSITION ='append', IOSTAT=ios)
          WRITE (UNIT=7,FMT=90) deltaT, th, net_mom, n_rb_aom,&
                                r_int_grav, r_int_conv, r_int_inert, total_energy
          CLOSE (UNIT=7)
          
          OPEN (UNIT=7,FILE='time_step.dat', POSITION ='append', IOSTAT=ios)
          WRITE (UNIT=7,FMT=90) deltaT, dtvisc, dtadv, dtg, &
                                dtsigma, dtinternal, stability
          CLOSE (UNIT=7)

          if ((rigidbody .EQV. .TRUE.) .and. (int_type == 1)) then
              OPEN (UNIT=7,FILE='disk.dat', POSITION ='append', IOSTAT=ios)
              WRITE (UNIT=7,FMT=90) deltaT, th, small_theta_n, total_energy,& 
                                    wedge_work, int_work, n_rb_aom, int_torque
              CLOSE (UNIT=7)

          elseif ((rigidbody .EQV. .TRUE.) .and. ((int_type == 2) .or. (int_type == 3))) then
              OPEN (UNIT=7,FILE='trans.dat', POSITION ='append', IOSTAT=ios)
              WRITE (UNIT=7,FMT=90) deltaT, delta_dpx, delta_dpy, total_energy,&
                                    wedge_work, int_work, int_forcex, int_forcey
              CLOSE (UNIT=7)
          
          endif

          !---- FIND MAX VELOCITY ----
          ! finding the stable timestep
          umax = max(unorth,usouth)
          vmax = max(vwest,veast)
          do i=1,Nx+1
            do j=1,Ny+2
              umax = max(abs(u(i,j)),umax) ! largest u velocity
            enddo
          enddo
          do i=1,Nx+2 
            do j=1,Ny+1
              vmax = max(abs(v(i,j)),vmax) ! largest velocity
            enddo
          enddo                

          !--- NEXT TIME STEP ---  
          ! if solution blows up or wrinkles form, then reduce stability parameter
          dtvisc=1.0/( max(mu_one/rho_one,mu_two/rho_two)*&
                       (2.0/(deltaX**2.0) + 2.0/(deltaY**2.0)) )
          dtadv = 1.0/(umax/deltaX + vmax/deltaY)
          dtg=sqrt( min(deltaX,deltaY)/max(gx,gy,1.0e-16) )
          dtsigma=sqrt( min(rho_one,rho_two)*min(deltaX,deltaY)**3.0/&
                        (4.0*pi*sigma+1.0e-16) )

          if ((rigidbody .EQV. .TRUE.) .and. (int_type == 2) .and. &
              (kx_adj /= 0.0) .and. (ky_adj /= 0.0)) then

             dtinternal = 0.5 * min( sqrt(min(rb1_mass,slider_mass)/kx_adj), &
                                     sqrt(min(rb1_mass,slider_mass)/ky_adj) )
          endif            

          deltaT=min(dtvisc,dtadv,dtg,dtsigma,dtinternal)
          deltaT=stability*deltaT  ! stability must be less than 0.9
          deltaT = min(deltaT,time_print)

!          print *, " "
!          print *,  "time =", time, " deltaT =", deltaT
!          print *, " "
!          print *, "/"
!          print *, " "
!          print *, " Viscous dt =", dtvisc
!          print *, " "
!          print *, "/"        
!          print *, " "
!          print *, " Advective dt =", dtadv
!          print *, " "
!          print *, "/"
!          print *, " Sigma dt =", dtsigma
!          print *, " "
!          print *, "/"
!          print *, " "
!          print *, " umax =", umax
!          print *, " "
!          print *, "/"        
!          print *, " "
!          print *, " vmax =", vmax
!          print *, " "
!          print *, "/"

          ! Calculate Rigid Body motion
          if (rigidbody .EQV. .TRUE.) then
            
            ! Check for motion restrictions      
            if (RB_lock /= 1) then

              ! Compute Relative rotation vector 
              big_theta = deltaT*n_rb_vom + 0.5*deltaT**2*n_rb_aom
                
              c1 = (sin(abs(big_theta)))/abs(big_theta)
              c2 = 2*(sin(0.5*abs(big_theta))**2)/(big_theta**2)        
                
              ! The exponential map
              exp_big_theta(1,1) = 1.0 - c2*big_theta**2
              exp_big_theta(1,2) = -c1*big_theta
              exp_big_theta(2,1) =  c1*big_theta
              exp_big_theta(2,2) = 1.0 - c2*big_theta**2
                
              ! Rotational velocity at (n+1/2)
              n12_rb_vom = n_rb_vom + 0.5*deltaT/big_J*&
                                       (r_int_inert + r_int_conv + r_int_grav)
            endif

            ! Calculate the rigid body velocity at v(n+1/2)
            n12_rb_vx = n_rb_vx + 0.5*deltaT*n_rb_ax
            n12_rb_vy = n_rb_vy + 0.5*deltaT*n_rb_ay
            
            ! Calculate d(n+1)                
            n1_rb_dx = n_rb_dx + deltaT*n12_rb_vx 
            n1_rb_dy = n_rb_dy + deltaT*n12_rb_vy
            
            ! Check for motion restrictions   
            if (RB_lock /= 1) then 
              n1_rbdelom(1,1) = n_rbdelom(1,1)*exp_big_theta(1,1) +&
                                n_rbdelom(1,2)*exp_big_theta(2,1)
              n1_rbdelom(1,2) = n_rbdelom(1,1)*exp_big_theta(1,2) +&
                                n_rbdelom(1,2)*exp_big_theta(2,2)
              n1_rbdelom(2,1) = n_rbdelom(2,1)*exp_big_theta(1,1) +&
                                n_rbdelom(2,2)*exp_big_theta(2,1)
              n1_rbdelom(2,2) = n_rbdelom(2,1)*exp_big_theta(1,2) +&
                                n_rbdelom(2,2)*exp_big_theta(2,2)
            endif


            if (int_type == 1) then

              ! Internal rotational velocity at (n+1/2)
              n12_rb2_vom = n_rb2_vom + 0.5*deltaT*n_rb2_aom

            elseif (int_type == 2) then
              
              ! Internal translational velocity at (n+1/2)  
              n12_rb2_vx = n_rb2_vx + 0.5*deltaT*n_rb2_ax
              n12_rb2_vy = n_rb2_vy + 0.5*deltaT*n_rb2_ay

              ! Calculate internal displacement d(n+1)
              n1_rb2_dx = n_rb2_dx + deltaT*n12_rb2_vx 
              n1_rb2_dy = n_rb2_dy + deltaT*n12_rb2_vy        
            
              ! Calculate new CoM and MoI           
              ncgx = cgx + n1_rb_dx
              ncgy = cgy + n1_rb_dy
              scgx = cgx + n1_rb2_dx
              scgy = cgy + n1_rb2_dy
              int_data(2) = ncgx
              int_data(3) = ncgy
              int_data(4) = scgx
              int_data(5) = scgy 

            elseif (int_type == 3) then
              
              ! Internal translational velocity at (n+1/2)  
              ! n12_rb2_vx = n_rb2_vx + 0.5*deltaT*n_rb2_ax
              ! n12_rb2_vy = n_rb2_vy + 0.5*deltaT*n_rb2_ay
              delta_vpy = delta_vpy + 0.5*n_rb2_ay

              ! Calculate internal displacement d(n+1)
              ! n1_rb2_dx = n_rb2_dx + deltaT*n12_rb2_vx 
              ! n1_rb2_dy = n_rb2_dy + deltaT*n12_rb2_vy        
              delta_dpy = delta_dpy + deltaT*delta_vpy

              ! Adjust
              ! adj_factor = deltaT*(n12_rb_vx*cos(th)-n12_rb_vy*sin(th))
              ! n1_rb2_dx = n1_rb2_dx + (adj_factor*cos(th))
              ! n1_rb2_dy = n1_rb2_dy - (adj_factor*sin(th))

              ! Calculate new CoM and MoI
              n1_rb2_dx = - sin(th)*delta_dpy
              n1_rb2_dy =   cos(th)*delta_dpy
              ncgx = cgx + n1_rb_dx
              ncgy = cgy + n1_rb_dy
              scgx = ncgx + n1_rb2_dx
              scgy = ncgy + n1_rb2_dy
              netcgx = (rb1_mass*ncgx + slider_mass*scgx)/(rb1_mass + slider_mass)
              netcgy = (rb1_mass*ncgy + slider_mass*scgy)/(rb1_mass + slider_mass)
              old_big_J = big_J
              big_J = buoy_J + rb1_mass*((netcgx - ncgx)**2 + (netcgy - ncgy)**2) + &
                      slider_J + slider_mass*((netcgx - scgx)**2 + (netcgy - scgy)**2)
              diff_J = (big_J - old_big_J)/deltaT
              int_data(2) = ncgx
              int_data(3) = ncgy
              int_data(4) = scgx
              int_data(5) = scgy 
            
            endif     

            if (int_type /= 3) then
              netcgx = ncgx
              netcgy = ncgy
            endif

                      
          endif
                
          ! ////////////////////////// VELOCITY \\\\\\\\\\\\\\\\\\\\\\\\\\\
          ! apply advection BC's
          call BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                          unorth,usouth,ueast,uwest,&
                          vnorth,vsouth,veast,vwest,&
                          deltaX,deltaY,phiLS,contC,contA,contR,contM,1)
                                                
          ! ---- FORCE THE FLUID WITH THE RIGID BODY ---- !
                           
          if (rigidbody .EQV. .TRUE.) then          
            
            ! Enforce the rigid body boundary condition on fluid velocity                 
            do i=1,Nx+1
              do j=1,Ny+1
                                                
                if ( 0.5*(s_phiLS(i+1,j) + s_phiLS(i,j)) >= 0) then

                  ry = ywhole(j) - netcgy

                  ! new velocity, n+1/2 time level
                  u(i,j) = n12_rb_vx - n12_rb_vom*ry


                endif
              enddo
            enddo

            do i=1,Nx+2
              do j=1,Ny+1
                                
                if ( 0.5*(s_phiLS(i,j) + s_phiLS(i,j+1)) >= 0) then

                  rx = xwhole(i) - netcgx

                  ! new velocity, n+1/2 time level
                  v(i,j) = n12_rb_vy + n12_rb_vom*rx

                endif
              enddo
            enddo
          endif

          if (wedge .EQV. .TRUE.) then
            do i=1,Nx+1
              do j=1,Ny+1
                                                
                if ( 0.5*(w_phiLS(i+1,j) + w_phiLS(i,j)) >= 0) then

                ! new velocity, n+1/2 time level
                u(i,j) = 0.0

                endif
              enddo
            enddo
 

            do i=1,Nx+2
              do j=1,Ny+1
                                
                if ( 0.5*(w_phiLS(i,j) + w_phiLS(i,j+1)) >= 0) then

                ! new velocity, n+1/2 time level
                v(i,j) = n1_w_vy


                endif
              enddo
            enddo
          endif

                
          ! ---- CONVECTIVE TERMS WITH NO UPWINDING ---- !
          ! ---- Unstable for very convective problems ---- !
          ! calculating the convection term for the x-momentum equation
          !do j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
          !        do i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
          !
          !                ! dummy variable for the convective opperator
          !                ustar(i,j)=((u(i+1,j)+u(i,j))**2.0-(u(i,j)+u(i-1,j))**2.0)/deltaX*0.25 + &
          !                                   ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - &
          !                                   (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))/deltaY*0.25
          !                ! i corresponds to half nodes
          !                ! j corresponds to whole nodes        
          !        enddo
          !enddo
                
          ! calculating the convection term for the y-momentum equation
          !do j=2,Ny                ! j=1 and j=Ny+1 are boundary nodes for v
          !        do i=2,Nx+1        ! i=1 and i=Nx+2 are boundary nodes for v
          !
          !                ! dummy variable for the convective opperator
          !                vstar(i,j)=((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - &
          !                                   (u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/deltaX*0.25 + &
          !                                   ((v(i,j+1)+v(i,j))**2.0-(v(i,j)+v(i,j-1))**2.0)/deltaY*0.25
          !                ! i corresponds to whole nodes
          !                ! j corresponds to half nodes
          !
          !        enddo
          !enddo

          ! ---- CONVECTIVE TERMS WITH FIRST ORDER UPWINDING ---- !
          ! ---- Diffusive for high velocity gradients ---- !
          ! First order upwinded convective terms
          ! calculating the convection term for the x-momentum equation
          !do j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
          !        do i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
          !
          !                u_plus = 0.5*(u(i,j) + abs(u(i,j)))
          !                u_minus = 0.5*(u(i,j) - abs(u(i,j)))
          !
          !                v_plus = 0.5*(v(i,j) + abs(v(i,j)))
          !                v_minus = 0.5*(v(i,j) - abs(v(i,j)))
          !                
          !                dudx_minus = (u(i,j) - u(i-1,j))/deltaX
          !                dudx_plus  = (u(i+1,j) - u(i,j))/deltaX
          !
          !                dudy_minus = (u(i,j) - u(i,j-1))/deltaY
          !                dudy_plus  = (u(i,j+1) - u(i,j))/deltaY
          !
          !                ! dummy variable for the convective opperator
          !                ustar(i,j) = u_plus*dudx_minus + u_minus*dudx_plus + &
          !                             v_plus*dudy_minus + v_minus*dudy_plus
          !                
          !                ! i corresponds to half nodes
          !                ! j corresponds to whole nodes        
          !        enddo
          !enddo

          ! calculating the convection term for the y-momentum equation
          !do j=2,Ny                ! j=1 and j=Ny+1 are boundary nodes for v
          !        do i=2,Nx+1        ! i=1 and i=Nx+2 are boundary nodes for v
          !
          !
          !                u_plus = 0.5*(u(i,j) + abs(u(i,j)))
          !                u_minus = 0.5*(u(i,j) - abs(u(i,j)))
          !
          !                v_plus = 0.5*(v(i,j) + abs(v(i,j)))
          !                v_minus = 0.5*(v(i,j) - abs(v(i,j)))
          !                
          !                dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
          !                dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
          !
          !                dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
          !                dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
          !
          !                ! dummy variable for the convective opperator
          !                vstar(i,j) = u_plus*dvdx_minus + u_minus*dvdx_plus + &
          !                             v_plus*dvdy_minus + v_minus*dvdy_plus
          !                
          !                ! i corresponds to half nodes
          !                ! j corresponds to whole nodes        
          !        enddo
          !enddo

          ! ---- CONVECTIVE TERMS WITH SECOND ORDER UPWINDING ---- !
          ! First order near the boundaries
          ! Diffusive for high velocity gradients
          ! calculating the convection term for the x-momentum equation
          do j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
            do i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
                
              u_plus = 0.5*(u(i,j) + abs(u(i,j)))
              u_minus = 0.5*(u(i,j) - abs(u(i,j)))

              v_plus = 0.5*(v(i,j) + abs(v(i,j)))
              v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                                                        
              if (i == 2) then
                dudx_minus = (u(i,j) - u(i-1,j))/deltaX
                dudx_plus  = (u(i+1,j) - u(i,j))/deltaX
              elseif (i == Nx) then
                dudx_minus = (u(i,j) - u(i-1,j))/deltaX
                dudx_plus  = (u(i+1,j) - u(i,j))/deltaX
              else
                dudx_minus = (3*u(i,j) - 4*u(i-1,j) + u(i-2,j))/(2*deltaX)
                dudx_plus  = (-u(i+2,j) + 4*u(i+1,j) - 3*u(i,j))/(2*deltaX)
              endif                                   

              if (j == 2) then
                dudy_minus = (u(i,j) - u(i,j-1))/deltaY
                dudy_plus  = (u(i,j+1) - u(i,j))/deltaY
              elseif (j == Ny+1) then
                dudy_minus = (u(i,j) - u(i,j-1))/deltaY
                dudy_plus  = (u(i,j+1) - u(i,j))/deltaY
              else
                dudy_minus = (3*u(i,j) - 4*u(i,j-1) + u(i,j-2))/(2*deltaY)
                dudy_plus  = (-u(i,j+2) + 4*u(i,j+1) - 3*u(i,j))/(2*deltaY)
              endif
 
              ! dummy variable for the convective opperator
              ustar(i,j) = u_plus*dudx_minus + u_minus*dudx_plus + &
                           v_plus*dudy_minus + v_minus*dudy_plus
                                
              ! i corresponds to half nodes
              ! j corresponds to whole nodes        
            enddo
          enddo

          ! calculating the convection term for the y-momentum equation
          do j=2,Ny                ! j=1 and j=Ny+1 are boundary nodes for v
            do i=2,Nx+1        ! i=1 and i=Nx+2 are boundary nodes for v
                
                
              u_plus = 0.5*(u(i,j) + abs(u(i,j)))
              u_minus = 0.5*(u(i,j) - abs(u(i,j)))

              v_plus = 0.5*(v(i,j) + abs(v(i,j)))
              v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                               
              if (i == 2) then
                dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
                dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
              elseif (i == Nx+1) then
                dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
                dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
              else
                dvdx_minus = (3*v(i,j) - 4*v(i-1,j) + v(i-2,j))/(2*deltaX)
                dvdx_plus  = (-v(i+2,j) + 4*v(i+1,j) - 3*v(i,j))/(2*deltaX)
              endif                                   

              if (j == 2) then
                dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
                dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
              elseif (j == Ny) then
                dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
                dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
              else
                dvdy_minus = (3*v(i,j) - 4*v(i,j-1) + v(i,j-2))/(2*deltaY)
                dvdy_plus  = (-v(i,j+2) + 4*v(i,j+1) - 3*v(i,j))/(2*deltaY)
              endif

              ! dummy variable for the convective opperator
              vstar(i,j) = u_plus*dvdx_minus + u_minus*dvdx_plus + &
                           v_plus*dvdy_minus + v_minus*dvdy_plus
                                
              ! i corresponds to half nodes
              ! j corresponds to whole nodes        
            enddo
          enddo
                
!          print *, "finished with convection"
!          print *, "/"
        
          ! apply viscous BC's
          call BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                          unorth,usouth,ueast,uwest,&
                          vnorth,vsouth,veast,vwest,&
                          deltaX,deltaY,phiLS,contC,contA,contR,contM,2)

          ! finding the intermediate velocity field for the x-momentum equation
          !do i = 2,Nx
          !        do j = 2,Ny+1
          do j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
            do i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
                        
              ! finding the viscous term
              ! dummy variables
              u_E = u(i+1,j)                ! note, i corresponds to half nodes
              u_O = u(i,j)                ! "
              u_W = u(i-1,j)                ! "
              u_N = u(i,j+1)                ! "
              u_S = u(i,j-1)                ! "

              v_ne = v(i+1,j)                ! note, j corresponds to half nodes
              v_nw = v(i,j)                ! "
              v_se = v(i+1,j-1)        ! "
              v_sw = v(i,j-1)                ! "

              ! u derivatives at the respective locations, n,s,e,w
              u_xe = (u_E - u_O)/deltaX
              u_xw = (u_O - u_W)/deltaX
              u_yn = (u_N - u_O)/deltaY
              u_ys = (u_O - u_S)/deltaY

              ! v derivatives derivatives at the respective locations, n,s
              v_xn = (v_ne - v_nw)/deltaX
              v_xs = (v_se - v_sw)/deltaX

              ! viscosity at the respective locations, n,s,e,w
              mu_n = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + &
                           mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j)) + &
                           mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1)) + &
                           mu_one*H(i+1,j+1) + mu_two*(1.0-H(i+1,j+1)))
              mu_s = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + & 
                           mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j)) + &
                           mu_one*H(i,j-1) + mu_two*(1.0-H(i,j-1)) + &
                           mu_one*H(i+1,j-1) + mu_two*(1.0-H(i+1,j-1)))
              mu_e = mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j))
              mu_w = mu_one*H(i,j) + mu_two*(1.0-H(i,j))                                                         

                                
              ! the viscous term
              Visc = (2.0*mu_e*u_xe - 2.0*mu_w*u_xw)/deltaX + &
                     (mu_n*(u_yn + v_xn) - mu_s*(u_ys + v_xs))/deltaY                        
                     ! note, i corresponds to ihalf
                                
              ! density at the half node
              rho_O = 0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j))) + &
                      0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))
                                                
              !print *, i
              !print *, j
              !print *, "rho_O"
              !print *, rho_O
                                
              ! surface tension - only for fluid/vapor interface
              !if (trackphase .EQV. .TRUE.) then 
              !        ! get the curvature 
              !        kappa_e = LSCURVATURE(Nx,Ny,i+1,j,deltaX,deltaY,phiLS)
              !        kappa_w = LSCURVATURE(Nx,Ny,i,j,deltaX,deltaY,phiLS)
              !        ! body force
              !        Fsigma = sigma*0.5*(kappa_e+kappa_w)*(H(i+1,j)-H(i,j))/deltaX
              !else
                      Fsigma = 0.0
              !endif
                                
              ! collect integral terms, if the node is underneath
              ! the rigid body
                                

              ! ustar on left is intermediate velocity, ustar on right is the conv term
              ustar(i,j) = u(i,j) + deltaT*(-ustar(i,j) + Visc/rho_O - gx - Fsigma/rho_O)        
                                        
            enddo
          enddo ! end of calculating x-momentum's viscous term

          ! finding the intermediate velocity field for the y-momentum equation
          do j=2,Ny                        ! j=1 and j=Ny+1 are boundary nodes for v
            do i=2,Nx+1                ! i=1 and i=Nx+2 are boundary nodes for v
                
                
              ! dummy variables
              u_ne = u(i,j+1)                ! note, i corresponds to ihalf
              u_se = u(i,j)                ! "
              u_nw = u(i-1,j+1)        ! "
              u_sw = u(i-1,j)                ! "

              v_E = v(i+1,j)                ! note, j corresponds to jhalf
              v_O = v(i,j)                ! "
              v_W = v(i-1,j)                ! "
              v_N = v(i,j+1)                ! "
              v_S = v(i,j-1)                ! "

              ! u derivatives at the respective locations, e,w
              u_ye = (u_ne - u_se)/deltaY
              u_yw = (u_nw - u_sw)/deltaY

              ! v derivatives at the respective locations, n,s,e,w
              v_yn = (v_N - v_O)/deltaY
              v_ys = (v_O - v_S)/deltaY
              v_xe = (v_E - v_O)/deltaX
              v_xw = (v_O - v_W)/deltaX
                        
              ! viscosity at the respective locations, n,s,e,w
              mu_n = mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1))
              mu_s = mu_one*H(i,j) + mu_two*(1.0-H(i,j))
              mu_e = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + & 
                           mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j)) + &
                           mu_one*H(i+1,j+1) + mu_two*(1.0-H(i+1,j+1)) + &
                           mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1)))
              mu_w = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + &
                           mu_one*H(i-1,j) + mu_two*(1.0-H(i-1,j)) + &
                           mu_one*H(i-1,j+1) + mu_two*(1.0-H(i-1,j+1)) + &
                           mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1)) )                                                         

              Visc = (mu_e*(u_ye + v_xe) - mu_w*(u_yw + v_xw))/deltaX + &
                     (2.0*mu_n*v_yn - 2.0*mu_s*v_ys)/deltaY                                        
                     ! note, j corresponds to jhalf

              ! density at the half node
              rho_O = 0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j))) + &
                      0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
                                
              !print *, "rho_O"
              !print *, rho_O

              ! surface tension
              !if (trackphase .EQV. .TRUE.) then
              !        ! get the curvature
              !        kappa_n = LSCURVATURE(Nx,Ny,i,j+1,deltaX,deltaY,phiLS)
              !        kappa_s = LSCURVATURE(Nx,Ny,i,j,deltaX,deltaY,phiLS)
              !        ! body force
              !        Fsigma = sigma*0.5*(kappa_n+kappa_s)*(H(i,j+1)-H(i,j))/deltaY
              !else
                      Fsigma = 0.0
              !endif
                                
              !print *, "kappa_n"
              !print *, kappa_n
              !print *, "kappa_s"
              !print *, kappa_s

              !rho_O = rho_one
              !beta_O = beta_one
                                
              !print *, "convective vstarij"
              !print *, vstar(i,j)
                                
              ! vstar on left is intermediate velocity, vstar on right is the conv term
              vstar(i,j) = v(i,j) + deltaT*(-vstar(i,j) + Visc/rho_O - gy - Fsigma/rho_O)        
                                
              !print *, "vij"
              !print *, v(i,j)
                                
              !print *, "deltaT"
              !print *, deltaT
                                
              !print *, "Visc"
              !print *, Visc
                                
              !print *, "Fsigma"
              !print *, Fsigma
                                
              !print *, "intermediate vstarij"
              !print *, vstar(i,j)
                                        
            enddo
          enddo
                
          !print *, "intermediate velocity"
          !print *, "/"

          ! ///////////////////////// NEW PRESSURE \\\\\\\\\\\\\\\\\\\\\\\\\\\\
          if (solvertype==1) then ! solver 1 only works for single phase
                        
            ! solving for pressure at n+1 implicitly so all the variables on the RHS of
            ! the PPE must be at n+1 or in the case of velocity they must be at star

            ! this is the RHS of the pressure poission equation (PPE)
            PPERHS(1:Nx,1:Ny) = rho_one* &
                                ((ustar(2:Nx+1,2:Ny+1) - ustar(1:Nx,2:Ny+1))/deltaX + &
                                 (vstar(2:Nx+1,2:Ny+1) - vstar(2:Nx+1,1:Ny))/deltaY)/deltaT
                                                                
                                                                
            A=0.0                ! X(0)
            B=Lx                ! X(XL)
            M=Nx
            MBDCND=3        ! 3=DERIVATIVE SPECIFIED AT BOTH ENDS
            BDA=0.0
            BDB=0.0
            C=0.0                ! Y(0)
            D=Ly                ! Y(YL)
            N=Ny
            NBDCND=3        ! 3=DERIVATIVE SPECIFIED AT BOTH ENDS, 2=reference p=0 at y=0
            BDC=0.0
            BDD=0.0
            ELMBDA=0.0        ! LAMBDA COEFFICIENT IN HELMHOLTZ EQ. (ALWAYS ZERO FOR THIS PPE)
            IDIMF=Nx

            ! FISHPAK ROUTINES FOR PPE SOLUTION
            CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD, &
                         ELMBDA,PPERHS,IDIMF,PERTRB,IERROR,WORK)

            IF (IERROR.NE.0) THEN
               WRITE(*,*) 'ERROR IN FISHPACK IN HSTCRT, IERROR= ',IERROR
               PAUSE
            ENDIF

            IF (WORK(1).GE.100000.0) THEN
              WRITE(*,*) 'WORK ARRAY IN FISHPACK UNDERDIMENSIONED'
              PAUSE
            ENDIF

            P(2:nx+1,2:ny+1) = PPERHS(1:nx,1:ny)
                        
            ! Neumann BC's for pressure
            P(1:nx+2,1) = P(1:nx+2,2)
            P(1:nx+2,ny+2) = P(1:nx+2,ny+1)
            P(1,1:ny+2) = P(2,1:ny+2)
            P(nx+2,1:ny+2) = P(nx+1,1:ny+2)
            ! ---- END OF PPE SOLUTION ----
                
          endif
                
          if (solvertype==2) then
                        
            !--- Specifing the pressure BC type ---
            ! north
            if (BCnorth==3 .or. BCnorth==4) then
              ! Dirclet Pressure BC
              BCPnorth = 2
            else
              ! Neumann Pressure BC
              BCPnorth = 1
            endif
            ! south
            if (BCsouth==3 .or. BCsouth==4) then
              ! Dirclet Pressure BC
              BCPsouth = 2
            else
              ! Neumann Pressure BC
              BCPsouth = 1
            endif 
            ! east
            if (BCeast==3 .or. BCeast==4) then
              ! Dirclet Pressure BC
              BCPeast = 2
            else
              ! Neumann Pressure BC
              BCPeast = 1
            endif
            ! west
            if (BCwest==3 .or. BCwest==4) then
              ! Dirclet Pressure BC
              BCPwest = 2
            else
              ! Neumann Pressure BC
              BCPwest = 1
            endif
                        
                        
            ! this is the RHS of the pressure poission equation (PPE) 
            do i=1,Nx  ! looping over all interior nodes
              do j=1,Ny
                                        
                ! ---- PHASE TRANSITION ----
                ! get the normal vector
                call LSNORMAL(Nx,Ny,i+1,j+1,deltaX,deltaY,phiLS,normal)
                                        
                !print *, i
                !print *, j

                !---- PPE RHS ----
                deltau = ustar(i+1,j+1)-ustar(i,j+1)
                deltav = vstar(i+1,j+1)-vstar(i+1,j)
                !print *, "delta_u"
                !print *, deltau
                !print *, "delta_v"
                !print *, deltav
                ! JDS 6/17/2008 Removed masstrans term from PPERHS
                PPERHS(i,j) = ( deltau/deltaX + deltav/deltaY )/deltaT
                ! note, i in ustar corresponds to ihalf
                ! note, j in vstar corresponds to jhalf
                ! note, delstar is at whole nodes
                !print *, "deltaX"
                !print *, deltaX
                !print *, "deltaY"
                !print *, deltaY
                !print *, "deltaT"
                !print *, deltaT
                !print *, "PPERHSij"
                !print *, PPERHS(i,j)
                !call ArrtoVec(Nx,Ny,PPERHS,junk)
                !junk2 = norm(Nx,Ny,junk)
                !print *, "norm of PPERHS"
                !print *, junk2
              enddo
            enddo ! end of the PPERHS                        
                        
            Pold = P
            call PCG(Nx,Ny,deltaX,deltaY,P(2:Nx+1,2:Ny+1),PPERHS,H,s_H,&
                     rho_one,rho_two,rho_three,&
                     tolerance,itermax,iterprint,precond,&
                     BCPnorth,BCPsouth,BCPeast,BCPwest,Pnorth,Psouth,Peast,Pwest)
                        

            ! Neumann BC's for pressure
            if (BCPsouth==1) then
              P(1:nx+2,1) = P(1:nx+2,2)
            endif
            if (BCPnorth==1) then
              P(1:nx+2,ny+2) = P(1:nx+2,ny+1)
            endif
            if (BCPwest==1) then
              P(1,1:ny+2) = P(2,1:ny+2)
            endif
            if (BCPeast==1) then
              P(nx+2,1:ny+2) = P(nx+1,1:ny+2)
            endif
                        
            ! Dirclet BC's for pressure
            if (BCPsouth==2) then
              P(1:nx+2,1) = Psouth
            endif
            if (BCPnorth==2) then
              P(1:nx+2,ny+2) = Pnorth
            endif
            if (BCPwest==2) then
              P(1,1:ny+2) = Pwest
            endif
            if (BCPeast==2) then
              P(nx+2,1:ny+2) = Peast
            endif
            ! ---- END OF PPE SOLUTION ----
          endif


                
          ! ---- SOR -----
          if (solvertype==3) then 
                        
            ! this is the RHS of the pressure poission equation (PPE) rho(2:Nx+1,2:Ny+1)
            PPERHS(1:Nx,1:Ny) = ((ustar(2:Nx+1,2:Ny+1) - ustar(1:Nx,2:Ny+1))/deltaX + &
                                 (vstar(2:Nx+1,2:Ny+1) - vstar(2:Nx+1,1:Ny))/deltaY)/deltaT
                                 ! note, i in ustar corresponds to ihalf
                                 ! note, j in vstar corresponds to jhalf
                                 ! note, delstar is at whole nodes
                        
            do iter = 1,itermax

              ! ---- APPLY BC's ----
              ! Neumann BC's for pressure
              P(1:nx+2,1) = P(1:nx+2,2)
              P(1:nx+2,ny+2) = P(1:nx+2,ny+1)
              P(1,1:ny+2) = P(2,1:ny+2)
              P(nx+2,1:ny+2) = P(nx+1,1:ny+2)

              ! corners (not used in derivative stencil)
              P(1,1) = 0.5*(P(2,1)+P(1,2))
              P(1,Ny+2)=0.5*(P(2,Ny+2)+P(1,Ny-1))
              P(Nx+2,1)=0.5*(P(Nx+1,1)+P(Nx+2,2))
              P(Nx+2,Ny+2)=0.5*(P(Nx+1,Ny+2)+P(Nx+2,Ny+1))
                                
                                 
              P(floor(real(Nx)/2.0)+1,1) = 0.0

              Pavg = 0.0        ! initializing back to zero
              maxPRes = 0.0 ! initializing back to zero
              do i=1,Nx
                do j=1,Ny
                  ! density 
                  rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
                          0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
                  rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
                          0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
                  rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
                          0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
                  rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
                          0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))
                                                
                  aw = 1.0/(rho_w*deltaX**2.0) !1.0/deltaX**2.0
                  ae = 1.0/(rho_e*deltaX**2.0) !1.0/deltaX**2.0
                  an = 1.0/(rho_n*deltaY**2.0) !1.0/deltaX**2.0
                  as = 1.0/(rho_s*deltaY**2.0) !1.0/deltaX**2.0
                  aO = -aw-ae-an-as                         !-4.0/deltaX**2.0

                  PRes = PPERHS(i,j)-ae*P(i+2,j+1)-aw*P(i,j+1) - &
                         an*P(i+1,j+2)-as*P(i+1,j)-ao*P(i+1,j+1)
                  P(i+1,j+1) = P(i+1,j+1) + omega/aO*PRes
                                                
                  maxPRes = max(abs(PRes),maxPRes)  ! largest error
                  Pavg = P(i+1,j+1) + Pavg        ! summing up the pressure
                enddo
              enddo
              Pavg = Pavg/(Nx*Ny)

              ! for Neumann BC's
              !do i=2,Nx+2
              !        do j=2,Ny+2
              !                P(i,j) = P(i,j)-Pavg  ! keeping the solution bounded
              !        enddo
              !enddo

              ! ---- check for convergence ----
              if (iter>1) then
                ! convergence criteria
                if (maxPRes<=tolerance) then
                  exit
                endif
              endif

              if (mod(iter,iterprint)==0) then
                print *, 'Pressure Residual =', maxPRes
              endif

              if (iter+1==itermax) then
                print *, 'ERROR-> Pressure Failed to Converge!'
              endif
                                        
              PavgOld = Pavg

            enddo  ! end of iteration
                
            ! ---- END OF PPE SOLUTION ----
          endif
                
          !print *, "PPE finished"
          !print *, "/"
                


          ! /////////////////////// NEW VELOCITY \\\\\\\\\\\\\\\\\\\\\\\\\\
          do i=2,Nx
            do j=2,Ny+1
              ! density at the half node
              rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j))) + &
                      0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))

              ! Hold old variables
              u_old(i,j) = u(i,j)

              ! new velocity, n+1 time level
              u(i,j) = ustar(i,j) - deltaT/rho_O*(P(i+1,j)-P(i,j))/deltaX  ! note, i corresponds to ihalf
            enddo
          enddo        
        
          do i=2,Nx+1
            do j=2,Ny
              ! density at the half node
              rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1))) + &
                      0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                
              ! Hold old variable
              v_old(i,j) = v(i,j)

              ! new velocity, n+1 time level
              v(i,j) = vstar(i,j) - deltaT/rho_O*(P(i,j+1)-P(i,j))/deltaY  ! note, j corresponds to jhalf
            enddo
          enddo
                
          !/////////////////////// NEW LEVEL-SET FIELD \\\\\\\\\\\\\\\\\\\\\\\\\\
          if (trackphase .EQV. .TRUE.) then
            ! All variables,(phiLS, u, v, P, etc) are at time level n, 
            ! but after the new level-set field is calculated the properties 
            ! will be at the n+1 time level.
                        
            uINT = u
            vINT = v

            !---- USING RK4 ----
            ! PDE => d(phiLS)/d(t) = -V*grad(phiLS)
            phiLSn = phiLS  ! storing the nth level values
            DO step=1,4
              ! getting the RHS                
              DO i=2,Nx+1
                DO j=2,Ny+1
                  ! calculate spatial derivatives
                  if (i<4 .or. j<4) then
                    CALL ENODIFF(Nx,Ny,i,j,deltaX,deltaY,u,v,phiLS,DiffX,DiffY)
                  elseif (i>Nx-2 .or. j>Ny-2) then
                    CALL ENODIFF(Nx,Ny,i,j,deltaX,deltaY,u,v,phiLS,DiffX,DiffY)
                  else
                    !CALL ENODIFF(Nx,Ny,i,j,deltaX,deltaY,u,v,phiLS,DiffX,DiffY)
                                                        
                    ! use WENO on interior (up to 5th order accurate)
                    CALL WENODIFF(Nx,Ny,i,j,deltaX,deltaY,u,v,phiLS,DiffX,DiffY)
                  endif
                                                
                                                
                  ! ---- PHASE TRANSITION ----
                  ! get the normal vector
                  !call LSNORMAL(Nx,Ny,i,j,deltaX,deltaY,phiLS,normal)
                                                
                  ! density at the O node
                  !rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                  !  JDS 7/16/2008 Took umflux and vmflux out of equation                                                
                  ! phiLSstar is a dummy variable for advection
                  phiLSstar(i,j) = -( 0.5*(u(i,j)+u(i-1,j)))*DiffX - &
                                    ( 0.5*(v(i,j)+v(i,j-1)))*DiffY                        

                ENDDO
              ENDDO
                                
              ! get next value
              DO i=2,Nx+1
                DO j=2,Ny+1
                  ! note, Coef is positive
                  phiLS(i,j) = phiLSn(i,j) + Coef(step)*deltaT*phiLSstar(i,j)                                
                ENDDO
              ENDDO

              ! apply boundary conditions (extrapolation)
              CALL BCLEVELSET(Nx,Ny,deltaX,deltaY,ywhole,phiLS,&
                              BCphisouth,BCphinorth,BCphiwest,BCphieast,d_phiLS)
            ENDDO
            ! phiLS is now at n+1
                                
            !---- REINITIALIZATION ----
            ! reinitialize the new LS field to a signed distance function
            !CALL LSREIN(Nx,Ny,deltaX,deltaY,phiLS,H,&
            !            BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)
                        
            !---------------------END OF STANDARD LSM---------------------

            ! ---- NEW HEAVIESIDE FCN ----
            ! finding the new heavieside function at each node
            do i=1,Nx+2
              do j=1,Ny+2
                ! storing H at level n for PPE, T, and the v(n+1) correction step
                !HOLD(i,j) = H(i,j)
                !s_HOLD(i,j) = s_H(i,j)

                ! heavieside function at n+1
                H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,phiLS)
                !s_H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,s_phiLS)
              enddo
            enddo
            ! all properties are now at n+1 if they use H

          endif  ! end of trackphase

                
          ! /////////// ADVANCE THE WEDGE - WHILE CALCULATING ENERGY INPUT \\\\\\\\\\\\\\\
          
          if (wedge .EQV. .TRUE.) then      
            
            if (time < 2) then      
                n12_w_vy = -2*pi*0.2*sin(2*pi*(time + 0.5*deltaT))
            else
                n12_w_vy = 0
            endif

            ! Calculate d(n+1)                
            !n1_rb_dx = n_rb_dx + deltaT*n12_rb_vx 
            n1_w_dy = n_w_dy + deltaT*n12_w_vy
                
                
            do i=1,Nx+2
              do j=1,Ny+2
                                
                if (xwhole(i) < 0.1) then
                  w_phiLS(i,j) = -(-ywhole(j) + 0.65 + n1_w_dy)
                else
                  w_phiLS(i,j) = -(-ywhole(j) + xwhole(i) + 0.55 + n1_w_dy)
                endif
                
                ! ---- NEW HEAVIESIDE FCN ----
                ! finding the new heavieside function at each node
                w_H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,w_phiLS)                                
                                                
              enddo
            enddo        

            do i=2,Nx+1
              do j=2,Ny
                        
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1))) + &
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))


                if (0.5*(w_H(i,j)+w_H(i,j+1)) > 0) then
                                        
                  ! density at the whole node
                  !rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))                        
                  !if (s_phiLS(i,j) >= 0) then
                                        
                  u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                  u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                  v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                  v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                               
                  if (i == 2) then
                    dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
                    dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
                  elseif (i == Nx+1) then
                    dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
                    dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
                  else
                    dvdx_minus = (3*v(i,j) - 4*v(i-1,j) + v(i-2,j))/(2*deltaX)
                    dvdx_plus  = (-v(i+2,j) + 4*v(i+1,j) - 3*v(i,j))/(2*deltaX)
                  endif                                   

                  if (j == 2) then
                    dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
                    dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
                  elseif (j == Ny) then
                    dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
                    dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
                  else
                    dvdy_minus = (3*v(i,j) - 4*v(i,j-1) + v(i,j-2))/(2*deltaY)
                    dvdy_plus  = (-v(i,j+2) + 4*v(i,j+1) - 3*v(i,j))/(2*deltaY)
                  endif

                  ! dummy variable for the convective opperator
                  vstar_int = u_plus*dvdx_minus + u_minus*dvdx_plus + &
                              v_plus*dvdy_minus + v_minus*dvdy_plus

                  if ((wedge .EQV. .TRUE.) .and. (0.5*(w_H(i,j)+w_H(i,j+1)) > 0)) then                                      
                   
                     n1_w_inert = n1_w_inert + &
                                     0.5*(w_H(i,j)+w_H(i,j+1))*&
                                     rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                     deltaX*deltaY
                                                                                        
                     n1_w_conv = n1_w_conv + &
                                    0.5*(w_H(i,j)+w_H(i,j+1))*&
                                    rho_O*vstar_int*deltaX*deltaY
                                
                     n1_w_grav = n1_w_grav + 0.5*&
                                    (w_H(i,j)+w_H(i,j+1))*rho_O*deltaX*deltaY*gy

                  endif                                
                endif
              enddo
            enddo

            wedge_force = n_w_inert + n_w_conv + n_w_grav
            wedge_work = wedge_work - (wedge_force * n1_w_vy * deltaT)
            
            if (time < 2) then
                n1_w_vy = -2*pi*0.2*sin(2*pi*time)
            else
                n1_w_vy = 0      
            endif

          endif
          
          ! Fluid Energy Contribution
          if (rigidbody .EQV. .FALSE.) then
            
            total_energy = 0.0
            lin_momx = 0.0
            lin_momy = 0.0

            do i=2,Nx+1
              do j=2,Ny+1
  
                if (i < Nx+1) then
                  rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))+&
                          0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                 
                  total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*deltaX*deltaY  ! kinetic (in x)
                  lin_momx = lin_momx + rho_O*u(i,j)*deltaX*deltaY
                endif
              
                if (j < Ny+1) then
                  rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                          0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))                          

                  total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*deltaX*deltaY  ! kinetic (in y)
                  lin_momy = lin_momy + rho_O*v(i,j)*deltaX*deltaY
                endif
 
                if ((i < Nx+1) .AND. (j < Ny+1)) then 
                  rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                  total_energy = total_energy + &                        ! Pressure and gravitational potential
                                 ( -P(i,j) + rho_O*gy*(j-1)*deltaY ) * &
                                 deltaX*deltaY
                endif

              enddo
            enddo
          endif

          ! /////////// RIGID BODY DYNAMICS - UPDATED VELOCITY \\\\\\\\\\\\\\\
                   
          ! Reconstruct the level set
          if (rigidbody .EQV. .TRUE.) then

            ! Recover change in rigid body's angle

            if (n1_rbdelom(1,1) < 0) then                 ! Expand domain of asin to +/-pi

                 th_domain = sign(pi,n1_rbdelom(2,1)) - asin(n1_rbdelom(2,1))
            else    
                 th_domain = asin(n1_rbdelom(2,1))            
            endif
            
            if ((abs(th_old) > pi/2.0) .and. &            ! Track number of revolutions 
               (sign(dble(1),th_old) /= sign(dble(1),th_domain))) then
                 
                 th_counter = th_counter + sign(dble(1),th_old)            
            endif

            th_change = th_domain + 2.0*pi*th_counter
            th_old = th_domain
            
            ! Level sets

            if (rb_shape == 1) then                                
              do i=1,Nx+2
                do j=1,Ny+2
  
                  th = 0.0 + th_change 
  
                  xc = b_rbLS + n1_rb_dx
                  yc = c_rbLS + n1_rb_dy
                                                
                  ! Circle
                  s_phiLS(i,j) = a_rbLS*sqrt( (xwhole(i)-(b_rbLS + n1_rb_dx))**2.0 + &
                                 (ywhole(j)-(c_rbLS + n1_rb_dy))**2.0 ) + d_rbLS
                                                
                                                
                enddo
              enddo

            elseif (rb_shape == 2) then                                
              do i=1,Nx+2
                do j=1,Ny+2
                                        
                  xn = xo + n1_rb_dx
                  yn = yo + n1_rb_dy
                  xc = xn
                  yc = yn
                 
                  th = pi/tho + th_change
                                                        
                  s_phiLS(i,j) = -100*&
                                 (bellipse**2*(xwhole(i)*cos(th) + ywhole(j)*sin(th) - &
                                               (xn*cos(th) + yn*sin(th)))**2 + &
                                  aellipse**2*(-xwhole(i)*sin(th) + ywhole(j)*cos(th) - &
                                               (-xn*sin(th) + yn*cos(th)))**2 - &
                                  aellipse**2*bellipse**2)
                                                
                enddo
              enddo

            elseif (rb_shape == 3) then
                        
              xc = c_xo + n1_rb_dx
              yc = c_yo + n1_rb_dy
              nx1 = x1 + n1_rb_dx
              ny1 = y1 + n1_rb_dy
              nx2 = x2 + n1_rb_dx
              ny2 = y2 + n1_rb_dy
              nx3 = x3 + n1_rb_dx
              ny3 = y3 + n1_rb_dy
              nx4 = x4 + n1_rb_dx
              ny4 = y4 + n1_rb_dy
              ncgx = cgx + n1_rb_dx
              ncgy = cgy + n1_rb_dy
                                
              th = c_tho + th_change

              cost = cos(th)
              sint = sin(th)
              rot_mat(1,1) = cost
              rot_mat(1,2) = -sint
              rot_mat(2,1) = sint
              rot_mat(2,2) = cost
                                        
              ! x1-4,y-4 are ordered counterclockwise
                                
              do r=1,4
                                
                if (r == 1) then
                                        
                  ! Region 1
                  x_v(1) = xc;
                  x_v(2) = nx1;
                  x_v(3) = nx2;
                                                
                  y_v(1) = yc;
                  y_v(2) = ny1;
                  y_v(3) = ny2;
                                                
                elseif (r == 2) then
                                        
                  ! Region 2
                  x_v(1) = xc;
                  x_v(2) = nx2;
                  x_v(3) = nx3;
                  
                  y_v(1) = yc;
                  y_v(2) = ny2;
                  y_v(3) = ny3;
                                                
                elseif (r == 3) then
                                        
                  ! Region 3
                  x_v(1) = xc;
                  x_v(2) = nx3;
                  x_v(3) = nx4;
                                                
                  y_v(1) = yc;
                  y_v(2) = ny3;
                  y_v(3) = ny4;
                                                
                elseif (r == 4) then
                                        
                  ! Region 4
                  x_v(1) = xc;
                  x_v(2) = nx4;
                  x_v(3) = nx1;
                                                
                  y_v(1) = yc;
                  y_v(2) = ny4;
                  y_v(3) = ny1;
                                                
                endif
                                        
                nx_v(1) = rot_mat(1,1)*(xc-ncgx) + rot_mat(1,2)*(yc-ncgy)
                ny_v(1) = rot_mat(2,1)*(xc-ncgx) + rot_mat(2,2)*(yc-ncgy)                        
                nx_v(2) = rot_mat(1,1)*(x_v(2)-ncgx) + rot_mat(1,2)*(y_v(2)-ncgy)
                ny_v(2) = rot_mat(2,1)*(x_v(2)-ncgx) + rot_mat(2,2)*(y_v(2)-ncgy)
                nx_v(3) = rot_mat(1,1)*(x_v(3)-ncgx) + rot_mat(1,2)*(y_v(3)-ncgy)
                ny_v(3) = rot_mat(2,1)*(x_v(3)-ncgx) + rot_mat(2,2)*(y_v(3)-ncgy)

                x_v(1) = nx_v(1) + ncgx
                y_v(1) = ny_v(1) + ncgy                        
                x_v(2) = nx_v(2) + ncgx
                y_v(2) = ny_v(2) + ncgy
                x_v(3) = nx_v(3) + ncgx
                y_v(3) = ny_v(3) + ncgy
                                                
                z_v(1) = 1.0;
                z_v(2) = 0.0;
                z_v(3) = 0.0;
                        
                Det = x_v(1)*(y_v(2)*z_v(3)-y_v(3)*z_v(2)) - &
                      y_v(1)*(x_v(2)*z_v(3)-x_v(3)*z_v(2)) + &
                      z_v(1)*(x_v(2)*y_v(3)-x_v(3)*y_v(2))
                                                          
                a_det = y_v(2)*z_v(3)-y_v(3)*z_v(2) - &
                        y_v(1)*(z_v(3)-z_v(2)) + &
                        z_v(1)*(y_v(3)-y_v(2))
                                                          
                b_det = x_v(1)*(z_v(3)-z_v(2)) - &
                        x_v(2)*z_v(3)-x_v(3)*z_v(2) + &
                        z_v(1)*(x_v(2)-x_v(3))
                                                          
                c_det = x_v(1)*(y_v(2)-y_v(3)) - &
                        y_v(1)*(x_v(2)-x_v(3)) + &
                        x_v(2)*y_v(3)-x_v(3)*y_v(2)
                                                                
                a_rbLS = a_det/Det
                b_rbLS = b_det/Det
                c_rbLS = c_det/Det
                                                
                do i=1,Nx+2
                  do j=1,Ny+2
                                        
                    ! Check on the location of the grid point.
                    ! Is it in the correct region?
                                                
                    pt(1) = xwhole(i)
                    pt(2) = ywhole(j)
            
                    t(1,1) = x_v(1);
                    t(2,1) = y_v(1);
                    t(1,2) = x_v(2);
                    t(2,2) = y_v(2);
                    t(1,3) = x_v(3);
                    t(2,3) = y_v(3);
            
                    do m = 1,3,2
                                                
                      k = mod(m,3) + 1
                                                        
                      temp = ( pt(1) - t(1,m) ) * ( t(2,k) - t(2,m) ) - &
                             ( pt(2) - t(2,m) ) * ( t(1,k) - t(1,m) )
                                                                   
                      if (0.0 < temp) then
                        inside = .FALSE.
                        exit
                      endif
                        inside = .TRUE.
                    enddo
                                                
                    if (inside .EQV. .TRUE.) then
                                                
                      s_phiLS(i,j) = (-a_rbLS*xwhole(i) - b_rbLS*ywhole(j) + 1)/c_rbLS
                                                                                   
                    endif
                  enddo                
                enddo
              enddo
            endif
   
                                        
            ! ---- NEW HEAVIESIDE FCN ----
            ! finding the new heavieside function at each node
            do i=1,Nx+2
              do j=1,Ny+2
                s_H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,s_phiLS)
              enddo
            enddo
                
            ! Recalculate the forces with the new displacement
      
            ! Collect some integrals
            
            total_energy = 0.0
            lin_momx = 0.0
            lin_momy = 0.0

            do i=2,Nx
              do j=2,Ny+1
                        
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j))) + &
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))

                total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*deltaX*deltaY  ! kinetic (in x)                                
                lin_momx = lin_momx + rho_O*u(i,j)*deltaX*deltaY               
                 
                if (0.5*(s_H(i,j)+s_H(i+1,j)) > 0)  then
                        
                                        
                  u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                  u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                  v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                  v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                                                        
                  if (i == 2) then
                    dudx_minus = (u(i,j) - u(i-1,j))/deltaX
                    dudx_plus  = (u(i+1,j) - u(i,j))/deltaX
                  elseif (i == Nx) then
                    dudx_minus = (u(i,j) - u(i-1,j))/deltaX
                    dudx_plus  = (u(i+1,j) - u(i,j))/deltaX
                  else
                    dudx_minus = (3*u(i,j) - 4*u(i-1,j) + u(i-2,j))/(2*deltaX)
                    dudx_plus  = (-u(i+2,j) + 4*u(i+1,j) - 3*u(i,j))/(2*deltaX)
                  endif                                   

                  if (j == 2) then
                    dudy_minus = (u(i,j) - u(i,j-1))/deltaY
                    dudy_plus  = (u(i,j+1) - u(i,j))/deltaY
                  elseif (j == Ny+1) then
                    dudy_minus = (u(i,j) - u(i,j-1))/deltaY
                    dudy_plus  = (u(i,j+1) - u(i,j))/deltaY
                  else
                    dudy_minus = (3*u(i,j) - 4*u(i,j-1) + u(i,j-2))/(2*deltaY)
                    dudy_plus  = (-u(i,j+2) + 4*u(i,j+1) - 3*u(i,j))/(2*deltaY)
                  endif
 
                  ! dummy variable for the convective operator
                  ustar_int = u_plus*dudx_minus + u_minus*dudx_plus + &
                              v_plus*dudy_minus + v_minus*dudy_plus

                        
                                                
                  n1_xint_inert = n1_xint_inert +&
                                  0.5*(s_H(i,j)+s_H(i+1,j))*&
                                  rho_O/DeltaT*(u(i,j) - u_old(i,j))*&
                                  deltaX*deltaY
                                                
                  n1_xint_conv = n1_xint_conv +&
                                 0.5*(s_H(i,j)+s_H(i+1,j))*&
                                 rho_O*ustar_int*deltaX*deltaY
                                
                  n1_xint_grav = n1_xint_grav + &
                                 0.5*(s_H(i,j)+s_H(i+1,j))*&
                                 rho_O*deltaX*deltaY*gx
                                                
                  ! Rotations
                  ! rs is the s-distance from the location at i,j to the centroid of the object.
                  ry = ywhole(j) - netcgy
                                                
                  r1_int_inert = r1_int_inert - &
                                 ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                 rho_O/deltaT*(u(i,j) - u_old(i,j))*&
                                 deltaX*deltaY
                                                
                  r1_int_conv = r1_int_conv - &
                                ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                rho_O*ustar_int*deltaX*deltaY
                                
                  r1_int_grav = r1_int_grav - &
                                ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                rho_O*deltaX*deltaY*gx
                                
                endif
              enddo
            enddo


            do i=2,Nx+1
              do j=2,Ny
                        
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1))) + &
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                
                total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*deltaX*deltaY ! kinetic (in y)
                lin_momy = lin_momy + rho_O*v(i,j)*deltaX*deltaY

                if (0.5*(s_H(i,j)+s_H(i,j+1)) > 0)then
                                        
                  ! density at the whole node
                  !rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))                        
                  !if (s_phiLS(i,j) >= 0) then
                                        
                  u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                  u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                  v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                  v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                               
                  if (i == 2) then
                    dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
                    dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
                  elseif (i == Nx+1) then
                    dvdx_minus = (v(i,j) - v(i-1,j))/deltaX
                    dvdx_plus  = (v(i+1,j) - v(i,j))/deltaX
                  else
                    dvdx_minus = (3*v(i,j) - 4*v(i-1,j) + v(i-2,j))/(2*deltaX)
                    dvdx_plus  = (-v(i+2,j) + 4*v(i+1,j) - 3*v(i,j))/(2*deltaX)
                  endif                                   

                  if (j == 2) then
                    dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
                    dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
                  elseif (j == Ny) then
                    dvdy_minus = (v(i,j) - v(i,j-1))/deltaY
                    dvdy_plus  = (v(i,j+1) - v(i,j))/deltaY
                  else
                    dvdy_minus = (3*v(i,j) - 4*v(i,j-1) + v(i,j-2))/(2*deltaY)
                    dvdy_plus  = (-v(i,j+2) + 4*v(i,j+1) - 3*v(i,j))/(2*deltaY)
                  endif

                  ! dummy variable for the convective opperator
                  vstar_int = u_plus*dvdx_minus + u_minus*dvdx_plus + &
                              v_plus*dvdy_minus + v_minus*dvdy_plus           

                  n1_yint_inert = n1_yint_inert + &
                                  0.5*(s_H(i,j)+s_H(i,j+1))*&
                                  rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                  deltaX*deltaY
                                                                                        
                  n1_yint_conv = n1_yint_conv + &
                                 0.5*(s_H(i,j)+s_H(i,j+1))*&
                                 rho_O*vstar_int*deltaX*deltaY
                                
                  n1_yint_grav = n1_yint_grav + 0.5*&
                                 (s_H(i,j)+s_H(i,j+1))*rho_O*deltaX*deltaY*gy
                                                
                  ! Rotations
                  ! rs is the s-distance from the location at i,j to the centroid of the object.
                  rx = xwhole(i) - netcgx
                                                
                  r1_int_inert = r1_int_inert + &
                                 rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                 rho_O/deltaT*(v(i,j) - v_old(i,j))*&
                                 deltaX*deltaY
                                                                                           
                  r1_int_conv = r1_int_conv + &
                                rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                rho_O*vstar_int*deltaX*deltaY
                                
                  r1_int_grav = r1_int_grav + &
                                rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                rho_O*deltaX*deltaY*gy
                                            
                endif
              enddo
            enddo
                
            ! Pressure and gravitational contributions to total energy
            do i=2,Nx
              do j=2,Ny     
              
                ! density at the whole node
                rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                total_energy = total_energy + &
                               ( -P(i,j) + rho_O*gy*(j-1)*deltaY ) * &
                               deltaX*deltaY

              enddo
            enddo

            ! Upddate the convective angular velocity
            n1_rb_vom = n_rb_vom + 0.5*deltaT/big_J*&
                                   (r_int_inert + r_int_conv + r_int_grav) &
                                 + 0.5*deltaT/big_J*&
                                   (r1_int_inert + r1_int_conv + r1_int_grav)
                
            ! Update internal rotational elements
            if (int_type == 1) then
              
              small_theta_n1 = small_theta_n + deltaT*n12_rb2_vom 
              delta_theta = small_theta_n1 - th
              delta_vom = n12_rb2_vom - n12_rb_vom  
              int_torque = k12*delta_theta + c12*delta_vom 
              int_data(2) = ncgx                                 ! Save data for post-processing
              int_data(3) = ncgy
              int_data(4) = small_theta_n

              ! Acceleration
              n1_rb2_aom = -int_torque/small_J

              ! Update velocity
              n1_rb2_vom = n12_rb2_vom + 0.5*deltaT*n1_rb2_aom
           
              ! Energy update
              int_work = int_work + c12*(delta_vom**2)*deltaT
              total_energy = total_energy + &
                             0.5*small_J*(n1_rb2_vom**2) + &
                             0.5*k12*(delta_theta**2)

            elseif (int_type == 2) then

              ! Coordinate Transformation
              delta_vpx =  cos(th)*(n12_rb2_vx - n12_rb_vx) + sin(th)*(n12_rb2_vy - n12_rb_vy)
              delta_vpy = -sin(th)*(n12_rb2_vx - n12_rb_vx) + cos(th)*(n12_rb2_vy - n12_rb_vy)

              delta_dpx =  cos(th)*(n1_rb2_dx - n1_rb_dx) + sin(th)*(n1_rb2_dy - n1_rb_dy)
              delta_dpy = -sin(th)*(n1_rb2_dx - n1_rb_dx) + cos(th)*(n1_rb2_dy - n1_rb_dy)

              if ((delta_dpx >= 0.8*max_e) .or. (delta_dpx <= 0.8*max_w)) then
                      kx_adj = 100*kx12
              else
                      kx_adj = (.64*kx12*max_e*max_w) / ((delta_dpx - 0.8*max_e)*(delta_dpx - 0.8*max_w))
              endif

              if ((delta_dpy >= 0.8*max_n) .or. (delta_dpy <= 0.8*max_s)) then
                      ky_adj = 100*ky12
              else
                      ky_adj = (.64*ky12*max_n*max_s) / ((delta_dpy - 0.8*max_n)*(delta_dpy - 0.8*max_s))
              endif

              int_forcepx = kx_adj*delta_dpx + cx12*delta_vpx
              int_forcepy = ky_adj*delta_dpy + cy12*delta_vpy

              ! Inverse transform
              int_forcex = cos(th)*int_forcepx - sin(th)*int_forcepy
              int_forcey = sin(th)*int_forcepx + cos(th)*int_forcepy

              ! Acceleration
              n1_rb2_ax = - int_forcex/slider_mass
              n1_rb2_ay = - int_forcey/slider_mass - gy
             
              ! Update velocity
              n1_rb2_vx = n12_rb2_vx + 0.5*deltaT*n1_rb2_ax                
              n1_rb2_vy = n12_rb2_vy + 0.5*deltaT*n1_rb2_ay

              ! Energy
              int_work = int_work + deltaT * (cx12*(delta_vpx**2) + cy12*(delta_vpy**2))
              total_energy = total_energy + &
                             0.5*slider_mass*((n12_rb2_vx**2) + (n12_rb2_vy**2)) + &
                             0.5*kx12*(delta_dpx**2) + &
                             0.5*ky12*(delta_dpy**2)
              lin_momx = lin_momx + slider_mass*n12_rb2_vx
              lin_momy = lin_momy + slider_mass*n12_rb2_vy

            elseif (int_type == 3) then
             
              ! Coordinate Transformation
              ! delta_vpy = -sin(th)*(n1_rb2_vx - n1_rb_vx) + cos(th)*(n1_rb2_vy - n1_rb_vy)
              ! delta_dpy = -sin(th)*(n1_rb2_dx - n1_rb_dx) + cos(th)*(n1_rb2_dy - n1_rb_dy)

              if ((delta_dpy >= max_n) .or. (delta_dpy <= max_s)) then
                      ky_adj = 10000
              else
                      ky_adj = (ky12*max_n*max_s) / ((delta_dpy - max_n)*(delta_dpy - max_s))
              endif
 
              !!!!!!!!!!! rb_dpmove = cos(th)*(n1_rb_dx-n_rb_dx) + sin(th)*(n1_rb_dy-n_rb_dy)
              !!!!!!!!!!! rb_vpmove = cos(th)*(n1_rb_dx-n_rb_dx) + sin(th)*(n1_rb_dy-n_rb_dy)

              int_forcepy = ky12*delta_dpy + cy12*delta_vpy + cos(th)*slider_mass*gy
              
              ! External Forces
              ext_px =  cos(th)*(n1_xint_inert + n1_xint_conv + n1_xint_grav) + &
                        sin(th)*(n1_yint_inert + n1_yint_conv + n1_yint_grav + sin(th)*slider_mass*gy)
              ext_py = -sin(th)*(n1_xint_inert + n1_xint_conv + n1_xint_grav) + &
                        cos(th)*(n1_yint_inert + n1_yint_conv + n1_yint_grav + sin(th)*slider_mass*gy)

              ! Inverse Transformation
              int_forcex = - sin(th)*int_forcepy
              int_forcey =   cos(th)*int_forcepy
     
              ! Buoy Acceleration
              n1_rb_ax = cos(th)*(ext_px/(rb1_mass+slider_mass)) - sin(th)*(ext_py/rb1_mass) + &
                        int_forcex/rb1_mass - gx
              n1_rb_ay = sin(th)*(ext_px/(rb1_mass+slider_mass)) + cos(th)*(ext_py/rb1_mass) + &
                        int_forcey/rb1_mass - gy
              n1_rb_aom = 2/deltaT*(n1_rb_vom - n_rb_vom) - n_rb_aom + ((int_torque - n1_rb_vom*diff_J)/big_J)

              ! Internal Acceleration
              ! n1_rb2_ax = - int_forcex/slider_mass
              n1_rb2_ay = - int_forcepy/slider_mass

              ! Update velocity
              ! n1_rb2_vx = n12_rb2_vx + 0.5*deltaT*n1_rb2_ax                
              ! n1_rb2_vy = n12_rb2_vy + 0.5*deltaT*n1_rb2_ay
              delta_vpy = delta_vpy + 0.5*deltaT*n1_rb2_ay

              ! Energy
              int_work = int_work + deltaT * (cx12*(delta_vpx**2) + cy12*(delta_vpy**2))
              total_energy = total_energy + &
                             0.5*slider_mass*((n12_rb2_vx**2) + (n12_rb2_vy**2)) + &
                             0.5*kx12*(delta_dpx**2) + &
                             0.5*ky12*(delta_dpy**2)

            endif

            ! Compute updated acceleration ::                
            if (int_type /= 3) then
            n1_rb_ax = (n1_xint_inert + n1_xint_conv + n1_xint_grav - gx*rb_mass + int_forcex)/rb_mass
            n1_rb_ay = (n1_yint_inert + n1_yint_conv + n1_yint_grav - gy*rb_mass + int_forcey)/rb_mass
            n1_rb_aom = 2/deltaT*(n1_rb_vom - n_rb_vom) - n_rb_aom + ((int_torque - n1_rb_vom*diff_J)/big_J)
            endif

            ! Restrict motion
            if (RB_lock == 1) then
              n1_rb_ax = 0.0
              n1_rb_aom = 0.0      
            endif 

            ! Calculate final velocity ::         
            n1_rb_vx = n12_rb_vx + 0.5*deltaT*n1_rb_ax                
            n1_rb_vy = n12_rb_vy + 0.5*deltaT*n1_rb_ay

            ! Calculate Energy ::
            total_energy = total_energy + &
                         rb_mass*gy*ncgy + &
                         0.5*rb_mass*((n1_rb_vx**2)+(n1_rb_vy**2)) + &
                         0.5*big_J*(n1_rb_vom**2)
            lin_momx = lin_momx + rb_mass*n1_rb_vx
            lin_momy = lin_momy + rb_mass*n1_rb_vy

          endif
                
                
          ! PRINT diveregent errors
          !do i=2,Nx
          !        do j=2,Ny
          !                if ( abs((u(i,j)-u(i-1,j))/deltaX+(v(i,j)-v(i,j-1))/deltaY) >= 1.0e-10 ) then
          !                        print *, (u(i,j)-u(i-1,j))/deltaX+(v(i,j)-v(i,j-1))/deltaY 
          !                endif
          !        enddo
          !enddo



          ! CHECK FOR INSTABILITY due to internal body parameters (like very large k or c) by detecting
          ! large, rapid oscillations in the acceleration of the buoy.
          if (int_type == 1) then
             
             if (abs(n_rb_aom - n1_rb_aom) > aom_limit) then
                aom_counter = aom_counter + 1
             else 
                aom_counter = 0
             endif
             
             if (aom_counter > 15) then
                stability  = 0.9*stability
                aom_limit = aom_limit + 10.0
                aom_counter = 0
             endif

          endif

          !---- INCREMENT THE TIME ----
          
          ! n+1 -> n
                
          n_xint_conv  = n1_xint_conv
          n_xint_grav  = n1_xint_grav
          n_xint_inert = n1_xint_inert
        
          n_yint_conv  = n1_yint_conv
          n_yint_grav  = n1_yint_grav
          n_yint_inert = n1_yint_inert
                
          r_int_conv  = r1_int_conv
          r_int_grav  = r1_int_grav
          r_int_inert = r1_int_inert
          
          n_w_conv  = n1_w_conv
          n_w_grav  = n1_w_grav
          n_w_inert = n1_w_inert
                      
          ! Reinitialize intergral variables two different directions
                
          n1_xint_conv  = 0.0
          n1_xint_grav  = 0.0
          n1_xint_inert = 0.0
        
          n1_yint_conv  = 0.0
          n1_yint_grav  = 0.0
          n1_yint_inert = 0.0
                
          r1_int_conv  = 0.0
          r1_int_grav  = 0.0
          r1_int_inert = 0.0

          n1_w_conv  = 0.0
          n1_w_grav  = 0.0
          n1_w_inert = 0.0

          n_w_dy = n1_w_dy
                
          n_rb_dx = n1_rb_dx
          n_rb_dy = n1_rb_dy
        
          n_rb_vx = n1_rb_vx
          n_rb_vy = n1_rb_vy

          n_rb_ax = n1_rb_ax
          n_rb_ay = n1_rb_ay
                
          n_rbdelom(1,1) = n1_rbdelom(1,1)
          n_rbdelom(1,2) = n1_rbdelom(1,2)
          n_rbdelom(2,1) = n1_rbdelom(2,1)
          n_rbdelom(2,2) = n1_rbdelom(2,2)
                
          n_rb_vom = n1_rb_vom
          n_rb_aom = n1_rb_aom

          n_rb2_vom = n1_rb2_vom
          n_rb2_aom = n1_rb2_aom

          small_theta_n = small_theta_n1

          n_rb2_dx = n1_rb2_dx
          n_rb2_dy = n1_rb2_dy
        
          n_rb2_vx = n1_rb2_vx
          n_rb2_vy = n1_rb2_vy

          n_rb2_ax = n1_rb2_ax
          n_rb2_ay = n1_rb2_ay
               
          time = time + deltaT  ! increment the time to n+1
          time_count = time_count + deltaT  ! increment the print clock

          !---- SAVING DATA ----
          !if (time>=time_max .or. time_count>=time_print) then
          if (time>=time_max .or. mod(time,time_print)<deltaT) then
            time_count = 0.0
            file_count = file_count + 1
                        
            ! ---- APPLY BC'S ----
            ! apply advection BC's
            !call BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
            !                unorth,usouth,ueast,uwest,&
            !                vnorth,vsouth,veast,vwest,&
            !                deltaX,deltaY,phiLS,contC,contA,contR,contM,1)
                                                
            CALL LSREIN(Nx,Ny,deltaX,deltaY,phiLS,H,&
                        BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)
            ! apply level set BC's prior to saving
            !call BCLEVELSET(Nx,Ny,phiLS,2,2,1,1)
            !CALL BCLEVELSET(Nx,Ny,deltaX,deltay,ywhole,phiLS,&
            !                BCphisouth,BCphinorth,BCphiwest,BCphieast,d_phiLS)
            !CALL BCLEVELSET(Nx,Ny,deltaX,deltay,ywhole,s_phiLS,&
            !                BCphisouth,BCphinorth,BCphiwest,BCphieast,d_phiLS)
                                                        
            call SAVEDATA(Nx,Ny,Lx,Ly,time,x,y,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v,file_count,int_data)

          endif        ! end of data saving
                
          !---- EXIT ----
          if (time>=time_max) exit  ! exit if t is greater than tmax

        enddo  ! end of time loop


        ! FORMATS for input file
        50  FORMAT (A)
        60  FORMAT (16X, I10)
        70  FORMAT (16X, E10.4)
        80  FORMAT (16X, L)
        
        90  FORMAT (10F13.6)  !9F10.4
        95  FORMAT (2I5)

END PROGRAM MPXLIB
!********************************************************************
!*                                                                                                                                        *
!*                                                                ENODIFF                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is for advecting scalar quantities T, PhiLS, etc. Note,                *
!* this program is only used near the boundaries when advecting the        *
!* the level set field because the WENODIFF is used in the interior.*
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE ENODIFF(Nx,Ny,i,j,deltaX,deltaY,u,v,VAR,DiffX,DiffY)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN):: Nx,Ny,i,j
        
        REAL(kind=8), INTENT(IN) :: deltaX,deltaY
        REAL(kind=8), INTENT(INOUT) :: DiffX,DiffY
        REAL(kind=8), DIMENSION(Nx+1,Ny+2), INTENT(IN) :: u
        REAL(kind=8), DIMENSION(Nx+2,Ny+1), INTENT(IN) :: v
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: VAR

        REAL(kind=8) :: Dnegx,Dnegy,Dnegz,Dposx,Dposy,Dposz
        REAL(kind=8) :: DiffRight,DiffLeft,DiffTop,DiffBottom
        REAL(kind=8) :: Dnegnegx,Dposnegx,Dposposx,Dnegnegy,Dposnegy,Dposposy


        ! EXTERNAL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: MFCN

        DiffX = 0.0
        DiffY = 0.0

        ! difference operators
        Dnegx = (VAR(i,j) - VAR(i-1,j))/deltaX ! Back diff at i
        Dposx = (VAR(i+1,j) - VAR(i,j))/deltaX ! Forward diff at i
        Dnegy = (VAR(i,j) - VAR(i,j-1))/deltaY ! Back diff at j
        Dposy = (VAR(i,j+1) - VAR(i,j))/deltaY ! Forward diff at j

        IF ((i==2 .OR. i==Nx+1) .OR. (j==2 .OR. j==Ny+1)) THEN
                ! near edges I use 1st order differences
                DiffLeft = Dnegx
                DiffRight = Dposx
                DiffBottom = Dnegy
                DiffTop = Dposy
        ELSE
                Dnegnegx = (VAR(i-2,j) - 2.0*VAR(i-1,j) + VAR(i,j))/deltaX**2.0  ! Central 2diff at i-1
                Dposnegx = (VAR(i-1,j) - 2.0*VAR(i,j) + VAR(i+1,j))/deltaX**2.0  ! Central 2diff at i
                Dposposx = (VAR(i,j) - 2.0*VAR(i+1,j) + VAR(i+2,j))/deltaX**2.0  ! Central 2diff at i+1
                Dnegnegy = (VAR(i,j-2) - 2.0*VAR(i,j-1) + VAR(i,j))/deltaY**2.0  ! Central 2diff at j-1
                Dposnegy = (VAR(i,j-1) - 2.0*VAR(i,j) + VAR(i,j+1))/deltaY**2.0  ! Central 2diff at j
                Dposposy = (VAR(i,j) - 2.0*VAR(i,j+1) + VAR(i,j+2))/deltaY**2.0  ! Central 2diff at j+1        

                DiffLeft = Dnegx + deltaX/2.0*MFCN(Dposnegx, Dnegnegx)    ! Diff on left face
                DiffRight = Dposx - deltaX/2.0*MFCN(Dposnegx, Dposposx)   ! Diff on right face
                DiffBottom = Dnegy + deltaY/2.0*MFCN(Dposnegy, Dnegnegy)  ! Diff on bottom face
                DiffTop = Dposy - deltaY/2.0*MFCN(Dposnegy, Dposposy)     ! Diff on top face
        ENDIF

        IF (0.5*(u(i,j)+u(i-1,j)) > 0.0) THEN
                DiffX = DiffLeft
        ELSE
                DiffX = DiffRight
        END IF

        IF (0.5*(v(i,j)+v(i,j-1)) > 0.0) THEN
                DiffY = DiffBottom
        ELSE
                DiffY = DiffTop
        END IF


END SUBROUTINE ENODIFF
!********************************************************************
!*                                                                                                                                        *
!*                                                                MFNC                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is the switch function of the 2nd order ENO schemes                        *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION MFCN(a,b)
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: a,b

        IF (a*b > 0.0) THEN
                IF (ABS(a) <= ABS(b)) THEN
                        MFCN = a;
                ELSEIF (ABS(a) > ABS(b)) THEN
                        MFCN = b;
                END IF
        ELSE
                MFCN = 0.0;
        END IF

END FUNCTION MFCN
!********************************************************************
!*                                                                                                                                        *
!*                                                          WENODIFF                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is for advecting the smooth scalar fields away from the                *
!* boundaries.  The scheme is up to 5th        order accurate and will                * 
!* always be more the 3rd order        accurate; the accuracy depends on        * 
!* the smoothness of the field.        For details see, Osher S., Fedkiw,  *
!* R. "Level Set Methods and Dynamic Implicit Surfaces", pp.33-37,        *
!* 2003.                                                                                                                        *
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE WENODIFF(Nx,Ny,i,j,deltaX,deltaY,u,v,VAR,DiffX,DiffY)

        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j
        
        REAL(kind=8), INTENT(IN) :: deltaX,deltaY
        REAL(kind=8), INTENT(OUT) :: DiffX,DiffY
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: VAR
        REAL(kind=8), DIMENSION(Nx+1,Ny+2), INTENT(IN) :: u
        REAL(kind=8), DIMENSION(Nx+2,Ny+1), INTENT(IN) :: v

        REAL(kind=8) :: vxone,vxtwo,vxthree,vxfour,vxfive
        REAL(kind=8) :: vyone,vytwo,vythree,vyfour,vyfive
        REAL(kind=8) :: varxone,varxtwo,varxthree,varyone,varytwo,varythree 
        REAL(kind=8) :: sxone,sxtwo,sxthree,syone,sytwo,sythree
        REAL(kind=8) :: axone,axtwo,axthree,ayone,aytwo,aythree,e
        REAL(kind=8) :: wone,wtwo,wthree


        ! X-direction
        
                if (0.5*(u(i,j)+u(i-1,j))<0.0) then
                        ! forward differences, u<0
                        vxone = (VAR(i+3,j) - VAR(i+2,j))/deltaX
                        vxtwo = (VAR(i+2,j) - VAR(i+1,j))/deltaX
                        vxthree = (VAR(i+1,j) - VAR(i,j))/deltaX
                        vxfour = (VAR(i,j) - VAR(i-1,j))/deltaX
                        vxfive = (VAR(i-1,j) - VAR(i-2,j))/deltaX
                else
                        ! backward differences, u>0
                        vxone = (VAR(i-2,j) - VAR(i-3,j))/deltaX
                        vxtwo = (VAR(i-1,j) - VAR(i-2,j))/deltaX
                        vxthree = (VAR(i,j) - VAR(i-1,j))/deltaX
                        vxfour = (VAR(i+1,j) - VAR(i,j))/deltaX
                        vxfive = (VAR(i+2,j) - VAR(i+1,j))/deltaX
                endif

                varxone = vxone/3.0 - 7.0/6.0*vxtwo + 11.0/6.0*vxthree
                varxtwo = -vxtwo/6.0 + 5.0/6.0*vxthree + vxfour/3.0
                varxthree = vxthree/3.0 + 5.0/6.0*vxfour - vxfive/6.0

                ! smoothness parameters
                sxone = 13.0/12.0*(vxone - 2.0*vxtwo + vxthree)**2.0 + &
                                         0.25*(vxone - 4.0*vxtwo + 3.0*vxthree)**2.0
                sxtwo = 13.0/12.0*(vxtwo - 2.0*vxthree + vxfour)**2.0 + 0.25*(vxtwo - vxfour)**2.0
                sxthree = 13.0/12.0*(vxthree - 2.0*vxfour + vxfive)**2.0 + &
                                           0.25*(3.0*vxthree - 4.0*vxfour + 3.0*vxfive)**2.0

                ! find smoothness ratios
                e = 0.000001*MAX(vxone**2.0,vxtwo**2.0,vxthree**2.0,vxfour**2.0,vxfive**2.0) + 1.0e-15
                axone = 0.1/(sxone + e)**2.0
                axtwo = 0.6/(sxtwo + e)**2.0
                axthree = 0.3/(sxthree + e)**2.0

                ! find weights, then find derivative
                wone = axone/(axone + axtwo + axthree)
                wtwo = axtwo/(axone + axtwo + axthree)
                wthree = axthree/(axone + axtwo + axthree)
                DiffX = wone*varxone + wtwo*varxtwo + wthree*varxthree

        
        ! Y-direction
        
                if (0.5*(v(i,j)+v(i,j-1))<0.0) then
                        ! forward differences, v<0
                        vyone = (VAR(i,j+3) - VAR(i,j+2))/deltaY
                        vytwo = (VAR(i,j+2) - VAR(i,j+1))/deltaY
                        vythree = (VAR(i,j+1) - VAR(i,j))/deltaY
                        vyfour = (VAR(i,j) - VAR(i,j-1))/deltaY
                        vyfive = (VAR(i,j-1) - VAR(i,j-2))/deltaY
                else
                        ! backward differences, v>0
                        vyone = (VAR(i,j-2) - VAR(i,j-3))/deltaY
                        vytwo = (VAR(i,j-1) - VAR(i,j-2))/deltaY
                        vythree = (VAR(i,j) - VAR(i,j-1))/deltaY
                        vyfour = (VAR(i,j+1) - VAR(i,j))/deltaY
                        vyfive = (VAR(i,j+2) - VAR(i,j+1))/deltaY
                endif

                varyone = vyone/3.0 - 7.0/6.0*vytwo + 11.0/6.0*vythree
                varytwo = -vytwo/6.0 + 5.0/6.0*vythree + vyfour/3.0
                varythree = vythree/3.0 + 5.0/6.0*vyfour - vyfive/6.0

                ! smoothness parameters
                syone = 13.0/12.0*(vyone - 2.0*vytwo + vythree)**2.0 + &
                                         0.25*(vyone - 4.0*vytwo + 3.0*vythree)**2.0
                sytwo = 13.0/12.0*(vytwo - 2.0*vythree + vyfour)**2.0 + 0.25*(vytwo - vyfour)**2.0
                sythree = 13.0/12.0*(vythree - 2.0*vyfour + vyfive)**2.0 + &
                                           0.25*(3.0*vythree - 4.0*vyfour + 3.0*vyfive)**2.0

                ! find smoothness ratios
                e = 0.000001*MAX(vyone**2.0,vytwo**2.0,vythree**2.0,vyfour**2.0,vyfive**2.0) + 1.0e-15
                ayone = 0.1/(syone + e)**2.0
                aytwo = 0.6/(sytwo + e)**2.0
                aythree = 0.3/(sythree + e)**2.0

                ! find weights, then find derivative
                wone = ayone/(ayone + aytwo + aythree)
                wtwo = aytwo/(ayone + aytwo + aythree)
                wthree = aythree/(ayone + aytwo + aythree)
                DiffY = wone*varyone + wtwo*varytwo + wthree*varythree

END SUBROUTINE WENODIFF
!********************************************************************
!*                                                                                                                                        *
!*                                                            LSREIN                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                          *
!*                                                                                                                                        *
!* This is a program for reinitializing to a sign distance function        *
!* using the technique suggested by, Sussman, M., Smereka, P. and        *
!* Osher, S.J., "A Level Set Method for Computing Solutions to                *
!* Solutions to Incompressible Two-Phase Flow, J. Computational                *
!* Physics, Vol. 114, pp.146-159, 1994.                                                                * 
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE LSREIN(Nx,Ny,deltaX,deltaY,phiLS,H,&
                                  BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,num)
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,BCphisouth,BCphinorth,BCphiwest,BCphieast
        
        INTEGER :: i,j,k,step,iter

        REAL(kind=8) :: deltaTau,u

        REAL(kind=8), INTENT(IN) :: deltaX,deltaY,num
        REAL(kind=8), DIMENSION(Ny+2), INTENT(IN) :: ywhole
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(INOUT) :: H
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(INOUT) :: phiLS
        REAL(kind=8), DIMENSION(Nx+2,Ny+2) :: phiLSn,RHS
        REAL(kind=8), DIMENSION(4) :: Coef
        
        ! EXTERANL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: DISTRHS
        REAL(kind=8), EXTERNAL :: LSHEAVY


        deltaTau = deltaX/2.0
        
        Coef(1) = 1.0/4.0
        Coef(2) = 1.0/3.0
        Coef(3) = 1.0/2.0
        Coef(4) = 1.0
        
        ! I need the heavieside fcn for reinitialization
        ! finding the new heavieside function at each node
        do i=1,Nx+2
                do j=1,Ny+2
                        ! heavieside function
                        H(i,j) = LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,phiLS)
                enddo
        enddo

        ! iter is the number of artificial time steps        
        DO iter=1,4

                !////////////////// RK4 //////////////////////////
                ! PDE => d(phi)/d(Tau) = u*(1 - ||grad(phi)||)
                ! note: u corresponds to the characteristic velocity
                phiLSn = phiLS
                DO step=1,4
                                                
                        DO i=2,Nx+1
                                DO j=2,Ny+1
                                        RHS(i,j) = DISTRHS(Nx,Ny,i,j,deltaX,deltaY,phiLS,H)                
                                ENDDO
                        ENDDO
                        
                        DO i=2,Nx+1
                                DO j=2,Ny+1
                                        ! note: Coef is positive
                                        phiLS(i,j) = phiLSn(i,j) + Coef(step)*deltaTau*RHS(i,j)                                
                                ENDDO
                        ENDDO
                        
                        ! apply boundary conditions
                        !call BCLEVELSET(Nx,Ny,phiLS,2,2,1,1)
                        call BCLEVELSET(Nx,Ny,deltaX,deltaY,ywhole,phiLS,&
                                                    BCphisouth,BCphinorth,BCphiwest,BCphieast,num)
                ENDDO

        ENDDO ! end artificial time loop


END SUBROUTINE LSREIN
!********************************************************************
!*                                                                                                                                        *
!*                                                          DISTRHS                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a function used in reinitializing the level set field to        *
!* a signed distance function.                                                                                *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION DISTRHS(Nx,Ny,i,j,deltaX,deltaY,phiLS,H)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j
        REAL(kind=8), INTENT(IN) :: deltaX,deltaY
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS,H
        REAL(kind=8) :: u, NormGrad

        ! EXTERANL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: LSHEAVY, NORMDIFFPHI

        
        u = 2.0*(H(i,j) - 1.0/2.0)  ! smoothed signum function
        NormGrad = NORMDIFFPHI(Nx,Ny,i,j,deltaX,deltaY,u,phiLS)
        
        DISTRHS = u*(1.0 - NormGrad)  ! u*(1 - ||grad(phi)||)

END FUNCTION DISTRHS
!********************************************************************
!*                                                                                                                                        *
!*                                                        NORMDIFFPHI                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the norm of the gradient of        *
!* phi using the method outlined in, Sethian, J.A., "Level Set                *
!* Methods and Fast Marching Methods", p 66, 2002                                        *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION NORMDIFFPHI(Nx,Ny,i,j,deltaX,deltaY,u,phi)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: deltaX,deltaY,u
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phi
        REAL(kind=8) :: Dnegx,Dnegy,Dposx,Dposy, &
                                        DiffRight,DiffLeft,DiffTop,DiffBottom, &
                                        Dnegnegx,Dposnegx,Dposposx, &
                                        Dnegnegy,Dposnegy,Dposposy

        ! EXTERNAL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: MFCN

        ! difference operators
        Dnegx = (phi(i,j) - phi(i-1,j))/deltaX ! Back diff at i
        Dposx = (phi(i+1,j) - phi(i,j))/deltaX ! Forward diff at i
        Dnegy = (phi(i,j) - phi(i,j-1))/deltaY ! Back diff at j
        Dposy = (phi(i,j+1) - phi(i,j))/deltaY ! Forward diff at j

        IF ((i==2 .OR. i==Nx+1) .OR. (j==2 .OR. j==Ny+1)) THEN
                ! near edges I use 1st order differences
                DiffLeft = Dnegx
                DiffRight = Dposx
                DiffBottom = Dnegy
                DiffTop = Dposy

        ELSE
                Dnegnegx = (phi(i-2,j) - 2.0*phi(i-1,j) + phi(i,j))/deltaX**2.0  ! Central 2diff at i-1
                Dposnegx = (phi(i-1,j) - 2.0*phi(i,j) + phi(i+1,j))/deltaX**2.0  ! Central 2diff at i
                Dposposx = (phi(i,j) - 2.0*phi(i+1,j) + phi(i+2,j))/deltaX**2.0  ! Central 2diff at i+1
                Dnegnegy = (phi(i,j-2) - 2.0*phi(i,j-1) + phi(i,j))/deltaY**2.0  ! Central 2diff at j-1
                Dposnegy = (phi(i,j-1) - 2.0*phi(i,j) + phi(i,j+1))/deltaY**2.0  ! Central 2diff at j
                Dposposy = (phi(i,j) - 2.0*phi(i,j+1) + phi(i,j+2))/deltaY**2.0  ! Central 2diff at j+1        

                ! function M is defined in Advection
                DiffLeft = Dnegx + deltaX/2.0*MFCN(Dposnegx, Dnegnegx)    ! Diff on left face
                DiffRight = Dposx - deltaX/2.0*MFCN(Dposnegx, Dposposx)   ! Diff on right face
                
                DiffBottom = Dnegy + deltaY/2.0*MFCN(Dposnegy, Dnegnegy)  ! Diff on bottom face
                DiffTop = Dposy - deltaY/2.0*MFCN(Dposnegy, Dposposy)     ! Diff on top face
                
        ENDIF

        IF (u >= 0.0) THEN
                NORMDIFFPHI = SQRT( (MAX(DiffLeft,0.0))**2.0 + &
                                                        (MIN(DiffRight,0.0))**2.0 + &
                                                        (MAX(DiffBottom,0.0))**2.0 + &
                                                        (MIN(DiffTop,0.0))**2.0 )
        ELSE
                NORMDIFFPHI = SQRT( (MIN(DiffLeft,0.0))**2.0 + &
                                                        (MAX(DiffRight,0.0))**2.0 + &
                                                        (MIN(DiffBottom,0.0))**2.0 + &
                                                        (MAX(DiffTop,0.0))**2.0 )
        ENDIF

END FUNCTION NORMDIFFPHI
!********************************************************************
!*                                                                                                                                        *
!*                                                          LSHEAVY                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the LS heavieside function                *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION LSHEAVY(Nx,Ny,i,j,deltaX,deltaY,phiLS)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: deltaX,deltaY

        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        
        REAL(kind=8) :: pi,epsilon

        PARAMETER(pi=3.141592653589)

        
        ! smearing the heavieside function over 1.5 cells in each direction 
        epsilon = 2.0*max(deltaX,deltaY) 

        IF (phiLS(i,j) >= epsilon) THEN
                LSHEAVY = 1.0
        ELSEIF (phiLS(i,j) <= -epsilon) THEN
                LSHEAVY = 0.0
        ELSE
                LSHEAVY = 0.5 + phiLS(i,j)/(2.0*epsilon) + &
                                  sin(pi*phiLS(i,j)/epsilon)/(2.0*pi)
        ENDIF

END FUNCTION LSHEAVY
!********************************************************************
!*                                                                                                                                        *
!*                                                          LSDELTA                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the LS delta function                        *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION LSDELTA(Nx,Ny,i,j,deltaX,deltaY,phiLS)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: deltaX,deltaY
        real(kind=8) :: epsilon

        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        
        REAL(kind=8) :: pi

        PARAMETER(pi=3.141592653589)

        ! smearing the delta function over 1.5 cells in each direction 
        epsilon = 2.0*max(deltaX,deltaY)

        IF (abs(phiLS(i,j)) > epsilon) THEN
                LSDELTA = 0.0
        ELSE
                LSDELTA = 1.0/(2.0*epsilon)*( 1.0 + cos(pi*phiLS(i,j)/epsilon) )
        ENDIF

END FUNCTION LSDELTA
!********************************************************************
!*                                                                                                                                        *
!*                                                        LSNORMAL                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the normal of the level set        *
!* field.                                                                                                                        *
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE LSNORMAL(Nx,Ny,i,j,deltaX,deltaY,phiLS,normal)

        IMPLICIT NONE
        
        INTEGER, INTENT(IN):: Nx,Ny,i,j 

        REAL(kind=8) :: phi_x,phi_y,norm
        
        REAL(kind=8), INTENT(IN) :: deltaX, deltaY
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        REAL(kind=8), DIMENSION(2), INTENT(OUT) :: normal

        ! the i,j values correspond to the whole nodes

        ! first derivatives
        phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*deltaX)
        phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*deltaY)

        ! magnitude of the the gradient of phi
        norm = (phi_x**2.0 + phi_y**2.0)**0.5

        ! normal vector
        normal(1) = phi_x/norm
        normal(2) = phi_y/norm

        !normal(1) = 0.0                        !??????
        !normal(2) = 1.0                        !??????

ENDSUBROUTINE LSNORMAL
!********************************************************************
!*                                                                                                                                        *
!*                                                        LSCURVATURE                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the curvature of the level set        *
!* field.                                                                                                                        *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION LSCURVATURE(Nx,Ny,i,j,deltaX,deltaY,phiLS)

        IMPLICIT NONE
        
        INTEGER, INTENT(IN):: Nx,Ny,i,j 

        REAL(kind=8) :: phi_x,phi_y,phi_xy,phi_xx,phi_yy,norm
        !REAL(kind=8), DIMENSION(2) :: normal
        
        REAL(kind=8), INTENT(IN) :: deltaX, deltaY
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        

        ! the i,j values correspond to the whole nodes

        ! first derivatives
        phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*deltaX)
        phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*deltaY)
        
        !print *, "phi_x"
        !print *, phi_x
        !print *, "phi_y"
        !print *, phi_y
        
        ! second derivatives
        phi_xx = ( phiLS(i+1,j) - 2.0*phiLS(i,j) + phiLS(i-1,j) )/(deltaX**2.0)
        phi_yy = ( phiLS(i,j+1) - 2.0*phiLS(i,j) + phiLS(i,j-1) )/(deltaY**2.0)

        ! mixed derivative
        phi_xy = (phiLS(i+1,j+1) - phiLS(i+1,j-1))/(4.0*deltaX*deltaY) - &
                         (phiLS(i-1,j+1) - phiLS(i-1,j-1))/(4.0*deltaX*deltaY)

        ! magnitude of the the gradient of phi
        norm = (phi_x**2.0 + phi_y**2.0)**0.5

        ! normal vector
        !normal(1) = phi_x/norm
        !normal(2) = phi_y/norm

        ! curvature
        LSCURVATURE = (phi_xx*phi_y**2.0 - 2.0*phi_y*phi_x*phi_xy + phi_yy*phi_x**2.0)/norm**3.0
                        
        ! bounding the curvature to the grid scale, (-1/deltaX<=kappa<=1/deltaX)
        if ( abs(LSCURVATURE) >= 1.0/min(deltaX,deltaY) ) then
                if (LSCURVATURE < 0.0) then
                        LSCURVATURE = -1.0/min(deltaX,deltaY)
                else
                        LSCURVATURE = 1.0/min(deltaX,deltaY)
                endif
        endif

END FUNCTION LSCURVATURE
!********************************************************************
!*                                                                                                                                        *
!*                                                        BCVELOCITY                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for applying boundary conditions to velocity        *                                                                                                                        
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                                          unorth,usouth,ueast,uwest,&
                                          vnorth,vsouth,veast,vwest,&
                                          deltaX,deltaY,phiLS,contC,contA,contR,contM,par)
        implicit none

        integer :: i,j
        integer, intent(in) :: Nx,Ny,par
        integer, intent(in) :: BCnorth,BCsouth,BCeast,BCwest

        real(kind=8), intent(in) :: unorth,usouth,ueast,uwest,&
                                                                vnorth,vsouth,veast,vwest
        real(kind=8), intent(in) :: deltaX,deltaY,contC,contA,contR,contM

        real(kind=8), dimension(Nx+1,Ny+2), intent(inout) :: u
        real(kind=8), dimension(Nx+2,Ny+1), intent(inout) :: v
        real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phiLS

        real(kind=8), external :: CONTACTSPEED

        ! par is a parameter to distinguish between viscous and advection BC's
        ! par = 1 corresponds to advection
        ! par = 2 corresponds to viscous

        if (par==1) then
                ! apply advection BC's

                ! ---- NORTH BOUNDARY ----
                if (BCnorth == 1) then                ! dirclet wall
                        v(2:Nx+1,Ny+1) = 0.0                                                        ! no flow through wall
                        u(1:Nx+1,Ny+2) = 2.0*unorth - u(1:Nx+1,Ny+1)        ! ghost node horizontal velocity
                elseif (BCnorth == 2) then        ! neumman wall
                        v(2:Nx+1,Ny+1) = 0.0                                                        ! no flow through wall
                        u(1:Nx+1,Ny+2) = u(1:Nx+1,Ny+1)                                        ! horizontal slip
                elseif (BCnorth == 3) then        ! free surface
                        v(2:Nx+1,Ny+1) = v(2:Nx+1,Ny)                                        ! constant velocity
                        u(1:Nx+1,Ny+2) = u(1:Nx+1,Ny+1)                                        ! constant velocity
                elseif (BCnorth == 4) then        ! pressure
                        v(2:Nx+1,Ny+1) = v(2:Nx+1,Ny) !2.0*v(2:Nx+1,Ny)-v(2:Nx+1,Ny-1) ! extrapolation
                        u(1:Nx+1,Ny+2) = u(1:Nx+1,Ny+1)        !2.0*u(1:Nx+1,Ny+1)-u(1:Nx+1,Ny) ! extrapolation
                elseif (BCnorth == 5) then !flow
                        v(2:Nx+1,Ny+1) = vnorth                                                        ! vertical velocity
                        u(1:Nx+1,Ny+2) = 2.0*unorth - u(1:Nx+1,Ny+1)        ! ghost node horizonal velocity
                elseif (BCnorth == 6) then        ! contact line
                        v(2:Nx+1,Ny+1) = 0.0                                                        ! no flow through wall
                        j=Ny+2
                        do i=1,Nx+1
                                u(i,j)=CONTACTSPEED(Nx,Ny,i,j,deltaX,deltaY,phiLS,contC,contA,contR,contM)
                        enddo
                        
                endif


                ! ---- SOUTH BOUNDARY ----
                if (BCsouth == 1) then                ! dirclet wall
                        v(2:Nx+1,1) = 0.0                                                                ! no flow through wall
                        u(1:Nx+1,1) = 2.0*usouth - u(1:Nx+1,2)                        ! ghost node horizontal velocity
                elseif (BCsouth == 2) then        ! neumman wall
                        v(2:Nx+1,1) = 0.0                                                                ! no flow through wall
                        u(1:Nx+1,1) = u(1:Nx+1,2)                                                ! horizontal slip
                elseif (BCsouth == 3) then        ! free surface
                        v(2:Nx+1,1) = v(2:Nx+1,2)                                                 ! constant velocity                                                                
                        u(1:Nx+1,1) = u(1:Nx+1,2)                                                ! constant velocity
                elseif (BCsouth == 4) then        ! pressure
                        ! 2.0*v(2:Nx+1,2)-v(2:Nx+1,3)
                elseif (BCsouth == 5) then        ! flow boundary
                        v(2:Nx+1,1) = vsouth                                                        ! vertical velocity
                        u(1:Nx+1,1) = 2.0*usouth - u(1:Nx+1,2)                        ! ghost node horizontal velocity
                elseif (BCnorth == 6) then        ! contact line
                        v(2:Nx+1,1) = 0.0                                                                ! no flow through wall
                        j=1
                        do i=1,Nx+1
                                u(i,j)=CONTACTSPEED(Nx,Ny,i,j,deltaX,deltaY,phiLS,contC,contA,contR,contM)
                        enddo
                        
                endif
                

                ! ---- EAST BOUNDARY ----
                if (BCeast == 1) then                ! dirclet wall        
                        u(Nx+1,2:Ny+1) = 0.0                                                        ! no flow through wall
                        v(Nx+2,1:Ny+1) = 2.0*veast - v(nx+1,1:ny+1)                ! ghost node vertical velocity
                elseif (BCeast == 2) then        ! neumman wall
                        u(Nx+1,2:Ny+1) = 0.0                                                        ! no flow through wall
                        v(Nx+2,1:Ny+1) = v(nx+1,1:ny+1)                                        ! vertical slip
                elseif (BCeast == 3) then        ! free surface
                        u(Nx+1,2:Ny+1) = u(Nx,2:Ny+1)                                        ! constant velocity                                                        
                        v(Nx+2,1:Ny+1) = v(Nx+1,1:ny+1)                                        ! constant velocity
                elseif (BCeast == 4) then        ! pressure
                        ! 2.0*u(Nx,2:Ny+1)-u(Nx-1,2:Ny+1)
                elseif (BCeast == 5) then        ! flow boundary
                        u(Nx+1,2:Ny+1) = ueast                                                        ! horizontal velocity
                        v(Nx+2,1:Ny+1) = 2.0*veast - v(nx+1,1:ny+1)                ! ghost node vertical velocity
                elseif (BCeast == 6) then        ! contact line
                        u(Nx+1,2:Ny+1) = 0.0                                                        ! no flow through wall
                        i=Nx+2
                        do j=1,Ny+1
                                v(i,j)=CONTACTSPEED(Nx,Ny,i,j,deltaX,deltaY,phiLS,contC,contA,contR,contM)
                        enddo
                                
                endif
                

                ! ---- WEST BOUNDARY ----        
                if (BCwest == 1) then                ! dirclet wall
                        u(1,2:Ny+1) = 0.0                                                                ! no flow through wall
                        v(1,1:Ny+1) = 2.0*vwest - v(2,1:ny+1)                        ! ghost node vertical velocity
                elseif (BCwest == 2) then        ! neumman wall
                        u(1,2:Ny+1) = 0.0                                                                ! no flow through wall
                        v(1,1:Ny+1) = v(2,1:ny+1)                                                ! vertical slip
                elseif (BCwest == 3) then        ! free surface
                        u(1,2:Ny+1) = u(2,2:Ny+1)                                                 ! constant velocity                                                                
                        v(1,1:Ny+1) = v(2,1:ny+1)                                                ! constant velocity
                elseif (BCeast == 4) then        ! pressure
                        ! 2.0*u(2,2:Ny+1)-u(3,2:Ny+1)
                elseif (BCwest == 5) then        !flow boundary
                        u(1,2:Ny+1) = uwest                                                                ! horizontal velocity
                        v(1,1:Ny+1) = 2.0*vwest - v(2,1:ny+1)                        ! ghost node vertical velocity
                elseif (BCwest == 6) then        ! contact line
                        u(1,2:Ny+1) = 0.0                                                                ! no flow through wall
                        i=1
                        do j=2,Ny+1
                                v(i,j)=CONTACTSPEED(Nx,Ny,i,j,deltaX,deltaY,phiLS,contC,contA,contR,contM)
                        enddo
                        
                endif

        elseif (par==2) then
                ! apply viscous BCs

                ! ---- NORTH BOUNDARY ----
                if (BCnorth == 1 .or. BCnorth == 5) then        ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        u(1:nx+1,ny+2)=8.0*unorth/3.0 - 2.0*u(1:nx+1,ny+1) + u(1:nx+1,ny)/3.0
                endif

                ! ---- SOUTH BOUNDARY ----
                if (BCsouth == 1 .or. BCsouth == 5) then        ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        u(1:nx+1,1)=8.0*usouth/3.0 - 2.0*u(1:nx+1,2) + u(1:nx+1,3)/3.0
                endif
                
                ! ---- EAST BOUNDARY ----
                if (BCeast == 1 .or. BCeast == 5) then                ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        v(nx+2,1:ny+1)=8.0*veast/3.0 - 2.0*v(nx+1,1:ny+1) + v(nx,1:ny+1)/3.0
                endif
                
                ! ---- WEST BOUNDARY ----
                if (BCwest == 1 .or. BCwest == 5) then                ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        v(1,1:ny+1)=8.0*vwest/3.0 - 2.0*v(2,1:ny+1) + v(3,1:ny+1)/3.0
                endif
                
        endif

ENDSUBROUTINE BCVELOCITY
!********************************************************************
!*                                                                                                                                        *
!*                                                        BCLEVELSET                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for apply BC's to the level-set field                        *                                                                                                                        
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE BCLEVELSET(Nx,Ny,deltaX,deltaY,ywhole,phiLS,&
                                          BCphisouth,BCphinorth,BCphiwest,BCphieast,num)
        implicit none
        
        ! local integers
        integer :: i,j

        ! integers
        integer, intent(in) :: Nx,Ny,BCphisouth,BCphinorth,BCphiwest,BCphieast
        
        ! local reals
        real(kind=8) :: deltafunc

        ! reals
        real(kind=8), intent(in) :: deltaX,deltaY,num
        
        ! real arrays
        real(kind=8), dimension(Nx+2), intent(in) :: ywhole
        real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phiLS
        
        ! external functions
        real(kind=8), external :: LSDELTA

        

        ! south wall
        if (BCphisouth == 1) then
                phiLS(2:Nx+1,1) = phiLS(2:Nx+1,2)                        ! neumman 
        else
                ! extrapolation
                phiLS(2:Nx+1,1) = 2.0*phiLS(2:Nx+1,2) - phiLS(2:Nx+1,3)
        endif
        
        ! north wall
        if (BCphinorth == 1) then
                phiLS(2:Nx+1,Ny+2) = phiLS(2:Nx+1,Ny+1)                ! neumman 
        else
                ! extrapolation
                phiLS(2:Nx+1,Ny+2) = 2.0*phiLS(2:Nx+1,Ny+1) - phiLS(2:Nx+1,Ny)
        endif

        ! west wall
        if (BCphiwest == 1) then
                phiLS(1,2:Ny+1) = phiLS(2,2:Ny+1)                        ! neumman
        elseif (BCphiwest == 2) then
                ! dirclet -> making phiLS = 0 at a point and neumman everywhere else
                do j=2,Ny+1
                        deltafunc = LSDELTA(Nx,Ny,1,j,deltaX,deltaY,phiLS)
                        if (deltafunc > 0.0) then
                                ! I am fixing the interface at ywhole=num
                                phiLS(1,j) = ywhole(j) + num        
                                phiLS(2,j) = ywhole(j) + num
                        else
                                phiLS(1,j) = phiLS(2,j)
                        endif        
                enddo                                                                        
        else
                ! extrapolation
                phiLS(1,2:Ny+1) = 2.0*phiLS(2,2:Ny+1) - phiLS(3,2:Ny+1)
        endif

        ! east wall
        if (BCphieast == 1) then
                phiLS(Nx+2,2:Ny+1) = phiLS(Nx+1,2:Ny+1)                ! neumman
        elseif (BCphieast == 2) then
                ! dirclet
                do j=2,Ny+1
                        deltafunc = LSDELTA(Nx,Ny,Nx+2,j,deltaX,deltaY,phiLS)
                        if (deltafunc > 0.0) then
                                ! I am fixing the interface at ywhole=num
                                phiLS(Nx+2,j) = ywhole(j) + num
                                phiLS(Nx+1,j) = ywhole(j) + num
                        else
                                phiLS(Nx+2,j) = phiLS(Nx+1,j)
                        endif
                enddo
        else
                ! extrapolation
                phiLS(Nx+2,2:Ny+1) = 2.0*phiLS(Nx+1,2:Ny+1) - phiLS(Nx,2:Ny+1)
        endif

        ! corner BC's
        phiLS(1,1) = 0.5*(phiLS(2,1) + phiLS(1,2))
        phiLS(Nx+2,1) = 0.5*(phiLS(Nx+1,1) + phiLS(Nx+2,2))
        phiLS(1,Ny+2) = 0.5*(phiLS(1,Ny+1) + phiLS(2,Ny+2))
        phiLS(Nx+2,Ny+2) = 0.5*(phiLS(Nx+1,Ny+2) + phiLS(Nx+2,Ny+1))

ENDSUBROUTINE BCLEVELSET
!********************************************************************
!*                                                                                                                                        *
!*                                                        LSCONTACT                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for moving the contact line along a boundary        *                                                                                                        
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION CONTACTSPEED(Nx,Ny,i,j,deltaX,deltaY,phiLS,&
                                                                   contC,contA,contR,contM)

        !possible error- check the sign of the angle.

        implicit none
        
        ! integers 
        integer, intent(in):: Nx,Ny,i,j

        real(kind=8) :: phi_x,phi_y,phi_xy,phi_xx,phi_yy,norm
        real(kind=8) :: pi,angle,Ucl,deltafcn,epsilon

        real(kind=8), intent(in) :: deltaX, deltaY
        real(kind=8), intent(in) :: contC,contA,contR,contM
        real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phiLS
        
        ! local parameters
        parameter (pi = 3.141592653589)
        
        ! external functions
        real(kind=8), external :: LSDELTA
        !
        !
        !
        epsilon = 1.5*max(deltaX,deltaY)  ! smearing the deltafcn 

        if (j == Ny+2) then                        ! north wall
                phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*deltaX) ! central difference
                phi_y = ( phiLS(i,j) - phiLS(i,j-1) )/deltaY                 ! backward difference
        elseif (j == 1) then                ! south wall
                phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*deltaX) ! central difference
                phi_y = ( phiLS(i,j+1) - phiLS(i,j) )/deltaY                 ! forward difference
        elseif (i == Nx+2) then                ! east wall
                phi_x = ( phiLS(i,j) - phiLS(i-1,j) )/deltaX                  ! backward difference
                phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*deltaY)  ! central difference
        elseif (i == 1) then                ! west wall
                phi_x = ( phiLS(i+1,j) - phiLS(i,j) )/deltaX                  ! forward difference
                phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*deltaY)  ! central difference
        endif

        ! normal vector: nx = phi_x/norm, ny = phi_y/norm
        angle = pi + ATAN2(phi_x,phi_y)  ! angle of contact line
                        
        ! delta function
        deltafcn = LSDELTA(Nx,Ny,i,j,deltaX,deltaY,phiLS)
                        
        ! contact line speed
        if (angle > contA) then
                CONTACTSPEED = deltafcn*contC*(angle-contA)**contM  
        elseif (angle < contR) then
                CONTACTSPEED = -deltafcn*contC*(contR-angle)**contM  
        else
                CONTACTSPEED = 0.0
        endif
                        
ENDFUNCTION CONTACTSPEED
!********************************************************************
!*                                                                                                                                        *
!*                                                          EXPINT                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a function used in find the exponential integral function*
!* EXPINT(x) = integral from x to infinity of exp(-t)/t for x>0                *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION EXPINT(x)
        implicit none
        
        integer :: i,Num

        real(kind=8), intent(in) :: x
        real(kind=8) :: Fone, Ftwo,sum,t,delta

        ! using the trapezoidal rule to integrate the function
        Num = 2000
        delta = (10.0-x)/real(Num)
        if (delta<0.0) then
                EXPINT = 0.0
        else
                sum = 0.0
                t=x
                do i=1,Num
                        Fone = exp(-t)/t
                        Ftwo = exp(-(t+delta))/(t+delta)
                        sum = sum + 0.5*delta*(Fone + Ftwo)
                        t=t+delta
                enddo
                EXPINT = sum
        endif
        if (x==0.0) then
                EXPINT = 1.0e16
        endif
        if (x>10.0) then
                EXPINT = 0.0
        endif
ENDFUNCTION EXPINT
!********************************************************************
!*                                                                                                                                        *
!*                                                                  ERF                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a function used in find the error function                                *
!* ERF(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2)                         *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION ERF(x)
        implicit none
        
        integer :: i,Num

        real(kind=8), intent(in) :: x
        real(kind=8) :: Fone, Ftwo,sum,t,delta,pi
        
        pi = 3.1415926535897932384626433832795

        ! using the trapezoidal rule to integrate the function
        Num = 2000
        delta = x/real(Num)
        if (delta<=0.0) then
                ERF = 0.0
        else
                sum = 0.0
                t=0.0
                do i=1,Num
                        Fone = 2.0/sqrt(pi)*exp(-t**2.0)
                        Ftwo = 2.0/sqrt(pi)*exp(-(t+delta)**2.0)
                        sum = sum + 0.5*delta*(Fone + Ftwo)
                        t=t+delta
                enddo
                ERF = sum
        endif
        
ENDFUNCTION ERF
