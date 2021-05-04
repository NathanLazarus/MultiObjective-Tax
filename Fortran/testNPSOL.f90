program testNPSOL

	implicit none


	! dimensions for state variables
	INTEGER, PARAMETER :: MAX_T = 60     ! time horizon

	! number of control variables for the optimization tool
	INTEGER, PARAMETER :: g_N_vector_vars = 3
	INTEGER, PARAMETER :: MAX_NCtrlVar = g_N_vector_vars * MAX_T + 1

	! constants
	double precision, PARAMETER :: g_CloseZero = 1.0d-12
	double precision, PARAMETER :: g_pInfty = 1.0d+15



	! derivative option for objective and constraint functions
	integer :: g_derOpt
	INTEGER, PARAMETER :: DERIVATIVE_OPT = 1

	integer :: g_NVar
	double precision, DIMENSION(:), allocatable :: x0
	double precision, dimension(:), allocatable :: xsol
	double precision :: fvalue
	! bound for states and controls
	double precision, DIMENSION(MAX_NCtrlVar) :: g_LB, g_UB ! lower and upper bounds for variables for the optimization tool 

	! NPSOL specific

	double precision, DIMENSION(:, :), allocatable :: g_linA 
	double precision, DIMENSION(:), allocatable :: g_linL, g_linU, g_nlnL, g_nlnU ! constraint for NPSOL: 
										! kmin <= F(k,L)-c+eps <= kmax

	INTEGER :: g_nclin ! number of linear constraints
	INTEGER :: g_ncnln ! number of nonlinear constraints
	INTEGER :: g_ldA ! row dimension of linear constraints 
	                 ! g_linL <= g_linA*x <= g_linU
	INTEGER :: g_ldJ ! dimension of Jacobian for nonlinear constraints 
	                 ! g_nlnL <= ConFun(x) <= g_nlnU
	INTEGER :: g_ldR ! row dimension of Hessian matrix
	! length of work vectors
	INTEGER :: g_leniw
	INTEGER :: g_lenrw
	INTEGER, DIMENSION(:), allocatable :: g_iw ! work vector
	double precision, DIMENSION(:), allocatable :: g_rw ! work vector
	INTEGER :: g_nctotl
	INTEGER, DIMENSION(:), allocatable :: g_istate 
	double precision, DIMENSION(:), allocatable :: g_bl, g_bu, g_clambda, g_grad, g_Con
	double precision, DIMENSION(:,:), allocatable :: g_cJac, g_Hess

	INTEGER :: g_capT, g_T_retirement
	double precision, dimension(:), allocatable :: g_wage
	double precision :: g_R, tC, tL, tK, assets_min, g_little_r, g_Lmax
	! model parameters
	double precision :: g_beta    ! discount factor
	double precision :: g_gamma   ! used in the definition of utility function
	double precision :: g_eta     ! used in the definition of utility function



	g_capT = 60
	g_T_retirement = 45


	allocate(g_wage(g_capT))
	g_wage = wageFunction(g_capT, g_T_retirement)
	g_beta = 0.94
	g_gamma = 0.5
	g_eta = 1

	g_little_r = 0.1
	g_R = 1 + g_little_r
	assets_min = -100
	g_Lmax = 1

	tC = 0.05
	tL = 0.1
	tK = 0


	g_derOpt = 0

	g_NVar = g_N_vector_vars * g_capT + 1
	g_nclin = 0
	g_ncnln = g_capT
	allocate(xsol(g_NVar))
	allocate(x0(g_NVar))

	g_LB(1:g_capT) = g_CloseZero !lower bound on consumption
	g_LB(g_capT + 1:2 * g_capT) = 0 !lower bound on labor
	g_LB(2 * g_capT + 1:3 * g_capT) = assets_min !-g_pInfty ! lower bound on assets (borrowing constraint)
	g_LB(3 * g_capT + 1) = 0

	g_UB(1:g_capT) = g_pInfty !upper bound on consumption
	g_UB(g_capT + 1:g_capT + g_T_retirement) = g_Lmax !upper bound on labor
	g_UB(g_capT + g_T_retirement + 1:2 * g_capT) = 0 !upper bound on labor
	g_UB(2 * g_capT + 1) = 0
	g_UB(2 * g_capT + 2:3 * g_capT + 1) = g_pInfty

	x0(1:g_NVar) = 1
	x0(g_capT + g_T_retirement + 1:2 * g_capT) = 0
	x0(2 * g_capT + 1) = 0
	x0(3 * g_capT + 1) = 0
	call OptimInit()
	call OptimMethod(g_NVar, x0, xsol, fvalue)



	write(*,*) 'c'
	write(*,*) xsol(1 : g_capT)
	write(*,*) 'l'
	write(*,*) xsol(g_capT + 1 : 2 * g_capT)
	write(*,*) 'k'
	write(*,*) xsol(2 * g_capT + 1 : 3 * g_capT + 1)

contains

	function wageFunction(T, Tret)

		integer, intent(in) :: T, Tret
		double precision :: wageFunction(T)
		integer :: i

		do i = 1, Tret
			wageFunction(i) = max(1.5d0, 1/2.d0 + i * (1 - (i + 0.d0)/T)) / 16
		end do

		wageFunction(Tret + 1:T) = 0

	end function wageFunction

	subroutine OptimInit()

		g_ldA = max(1, g_nclin)  ! row dimension of linear constraints 
		                         ! g_linL <= g_linA*x <= g_linU
		g_ldJ = max(1, g_ncnln)  ! dimension of Jacobian for nonlinear constraints 
		                         ! g_nlnL <= ConFun(x) <= g_nlnU
		g_ldR = g_NVar           ! row dimension of Hessian matrix
		allocate(g_linA(g_ldA, g_NVar))
		allocate(g_linL(g_ldA), g_linU(g_ldA))
		allocate(g_nlnL(g_ldJ), g_nlnU(g_ldJ)) 

		g_nctotl = g_NVar+g_nclin+g_ncnln 
		allocate(g_istate(g_nctotl))
		allocate(g_bl(g_nctotl), g_bu(g_nctotl), g_clambda(g_nctotl)) 
		allocate(g_grad(g_NVar), g_cJac(g_ldJ,g_NVar), g_Hess(g_ldR,g_NVar), g_Con(g_ldJ))

		! length of work vectors
		g_leniw = 3*g_NVar+g_nclin+2*g_ncnln
		g_lenrw = 2*g_NVar**2+g_NVar*g_nclin+2*g_NVar*g_ncnln+20*g_NVar+11*g_nclin+21*g_ncnln
		allocate(g_iw(g_leniw), g_rw(g_lenrw)) ! work vector

		g_nlnL(1:g_ldJ) = 0 ! 0<= F(k,L)-c-Ekplus <= 0
		g_nlnU(1:g_ldJ) = 0

		! set all bounds
		g_bl(1:g_NVar) = g_LB(1:g_NVar)
		g_bu(1:g_NVar) = g_UB(1:g_NVar)
		if (g_nclin > 0) then
			g_bl((g_NVar+1):(g_NVar+g_nclin)) = g_linL(1:g_nclin)
			g_bu((g_NVar+1):(g_NVar+g_nclin)) = g_linU(1:g_nclin)
		end if
		if (g_ncnln > 0) then
			g_bl((g_NVar+g_nclin+1):g_nctotl) = g_nlnL(1:g_ncnln)
			g_bu((g_NVar+g_nclin+1):g_nctotl) = g_nlnU(1:g_ncnln)
		end if

		! g_istate = 0 ! No constrains are expected to be active, only useful for Warm Start
		call npoptn('Cold Start')
		call npoptn('Nolist')
		call npoptn('Iteration limit = 1000')
		call npoptn('Print file = 0')
		if (g_derOpt == DERIVATIVE_OPT) then
			call npoptn('Derivative level = 3')  ! 3: all gradients of objective and constraints are given
		else
			call npoptn('Derivative level = 0')  ! 0: no gradient given
		end if

	end subroutine OptimInit

	subroutine OptimMethod(nvar, x0, xsol, fvalue)

		integer :: inform, iter, nvar
		double precision, intent(IN) :: x0(nvar) 
		double precision, intent(OUT) :: xsol(nvar), fvalue
		! double precision, external :: dF_k
		! EXTERNAL NPSOLConFun, NPSOLObjFun, dU_c

		xsol = x0

		call NPSOL(nvar,g_nclin,g_ncnln,g_ldA,g_ldJ,g_ldR,g_linA,g_bl,g_bu,NPSOLConFun,NPSOLObjFun,inform,& 
			iter,g_istate,g_Con,g_cJac,g_clambda,fvalue,g_grad,g_Hess,xsol,g_iw,g_leniw,g_rw,g_lenrw)

write(*,*) 'inform'
write(*,*) inform
write(*,*) 'iter'
write(*,*) iter

		select case (inform)
			case (9)
				print *, "An input parameter was invalid"
			case (7) 		
				print *, "The function derivatives returned by funcon or fcn is incorrect"
			case (6)
	!			g_WorkerInform = WORKER_NPSOL_1st_ORDER_FAIL
				print *, "x does not satisfy the 1st-order optimality conditions"
			case (4)
				print *, "The Major iteration limit was reached"
			case (3)
				print *, "The nonlinear constraints could not be satisfied"
			case (2)
				print *, "The linear constraints and bounds could not be satisfied"
	!		case (0)
	!			print *, "Get an optimal solution x"
			case default
	!			g_WorkerInform = WORKER_NPSOL_OTHER_FAIL
	!			print *, "Terminated with other reasons"
		end select

		! NPSOL solves the minimal problem, but we are dealing with the maximal problem
		! so we change the objective function f(x) to -f(x). Now we change the optimal 
		! objective value -f(x^*) backwards f(x^*)
		fvalue = -fvalue 

		if (abs(fvalue) > g_pInfty) then
			print *, "function value is NaN"
		end if

	end subroutine OptimMethod

	! nonlinear constraint function for npsol: g_nlnL <= NPSOLConFun(x) <= g_nlnU
	subroutine NPSOLConFun(mode,ncnln,n,ldJ,needc,x,cc,cJac,nstate)

		INTEGER :: mode,ncnln,n,ldJ,nstate, i
		INTEGER :: needc(ncnln)
		double precision :: x(n), cc(ncnln), cJac(ldJ,*)

		! do i = 1, DIM
		! 	cc(i) = x(i) - 4
		! 	! cc(i) = F(g_todayK(i),x(DIM+i),g_currthetas(i)) - x(i) - x(2*DIM+i)
		! end do

		! cc(1) = x(1) ** 2 - 1
		! if (g_derOpt == DERIVATIVE_OPT) then
		! 	cJac(1,1) = 2 * x(1)
		! end if

		if (g_derOpt == DERIVATIVE_OPT) then
			cJac(1:ldJ,1:n) = 0
			do i = 1, g_capT
				cJac(i, i) = - 1 - tC
				cJac(i, g_capT + i) = (1 - tL) * g_wage(i)
				cJac(i, 2 * g_capT + i) = 1 + (1 - tK) * (1 - g_R)
				cJac(i, 2 * g_capT + i + 1) = -1
			end do
		end if
		! else
		! 	do i = 1, DIM
		! 		cc(i) = x(i) + 1
		! 	end do
		! end if
		do i = 1, g_capT
			cc(i) = g_R * x(2 * g_capT + i) + g_wage(i) * x(g_capT + i) - x(i)  - x(2 * g_capT + i + 1) - &
					tC * x(i) - tL * g_wage(i) * x(g_capT + i) - tK * (g_R - 1) * x(2 * g_capT + i)
		end do


	end subroutine NPSOLConFun

	! the objective function for NPSOL method
	subroutine NPSOLObjFun(mode,n,x,v,grad,nstate)

		INTEGER :: mode,n,nstate
		double precision :: v
		double precision :: x(n), grad(n)	
		integer :: i

		v = 0
		if (g_derOpt == DERIVATIVE_OPT) then
			grad = 0
		end if

		do i = 1, g_capT

			v = v + (g_beta ** (i - 1)) * ((1/(1.d0 - g_gamma)) * x(i) ** (1.d0 - g_gamma) - &
										   (1/(1.d0 + g_eta)) * x(g_capT + i) ** (1.d0 + g_eta))

			if (g_derOpt == DERIVATIVE_OPT) then
				grad(i) = (g_beta ** (i - 1)) * x(i) ** (- g_gamma)
				grad(g_capT + i) = (g_beta ** (i - 1)) * x(g_capT + i) ** (g_eta)
			end if

		end do
		! NPSOL solves the minimal problem, but we are dealing with the maximal problem
		! so we change the objective function f(x) to -f(x)
		v = - v
		if (g_derOpt == DERIVATIVE_OPT) then
			grad = -grad
		end if

	end subroutine NPSOLObjFun




end program testNPSOL