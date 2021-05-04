!****************************************************************************
! THIS SOFTWARE IS COPYRIGHTED BY YONGYANG CAI AND KENNETH JUDD.

! THIS SOFTWARE IS PROVIDED TO BEN SKRAINKA AND MICHAEL WILDE FOR THE PURPOSE OF PORTING IT TO 
! SUPERCOMPUTERS AT ARGONNE NATIONAL LABS. OTHERS ASSISTING BEN SKRAINKA AND MICHAEL WILDE IN
!  THIS TASK MAY ALSO HAVE ACCESS TO THIS SOFTWARE. NO OTHERS MAY USE THIS CODE FOR ANY PURPOSE
! WITHOUT THE KNOWLEDGE AND CONSENT OF KENNETH JUDD. 

! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
! WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
! TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
! ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!****************************************************************************

! global variables for NPSOL
module Optim_global_vars

	USE Common_global_vars

	! parameters used in NPSOL 
	INTEGER, PARAMETER :: g_nclin=0 ! number of linear constraints
	INTEGER :: g_ncnln ! number of nonlinear constraints
	INTEGER :: g_ldA ! row dimension of linear constraints 
	                 ! g_linL <= g_linA*x <= g_linU
	INTEGER :: g_ldJ ! dimension of Jacobian for nonlinear constraints 
	                 ! g_nlnL <= ConFun(x) <= g_nlnU
	INTEGER :: g_ldR ! row dimension of Hessian matrix
	double precision, DIMENSION(:, :), allocatable :: g_linA 
	double precision, DIMENSION(:), allocatable :: g_linL, g_linU, g_nlnL, g_nlnU ! constraint for NPSOL: 
										! kmin <= F(k,L)-c+eps <= kmax

	INTEGER :: g_nctotl
	INTEGER, DIMENSION(:), allocatable :: g_istate 
	double precision, DIMENSION(:), allocatable :: g_bl, g_bu, g_clambda, g_grad, g_Con
	double precision, DIMENSION(:,:), allocatable :: g_cJac, g_Hess
	! length of work vectors
	INTEGER :: g_leniw
	INTEGER :: g_lenrw
	INTEGER, DIMENSION(:), allocatable :: g_iw ! work vector
	double precision, DIMENSION(:), allocatable :: g_rw ! work vector

end module Optim_global_vars


! initilization for NPSOL 
subroutine OptimInit()

	use Optim_global_vars
	IMPLICIT NONE

	g_ncnln = DIM ! number of nonlinear constraints: 0<= F(k,L)-c-Ekplus <= 0
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

	g_istate = 0 ! No constrains are expected to be active, only useful for Warm Start

	call npoptn('Nolist')
	call npoptn('Iteration limit = 1000')
	call npoptn('Print file = 0')
	if (g_derOpt == DERIVATIVE_OPT) then
		call npoptn('Derivative level = 3')  ! 3: all gradients of objective and constraints are given
	else
		call npoptn('Derivative level = 0')  ! 0: no gradient given
	end if

end subroutine OptimInit


! call NPSOL
subroutine OptimMethod(nvar, x0, xsol, fvalue)

	use Optim_global_vars
	IMPLICIT NONE

	integer :: inform, iter, nvar, i
	double precision, intent(IN) :: x0(nvar) 
	double precision, intent(OUT) :: xsol(nvar), fvalue
	double precision, external :: dF_k
	EXTERNAL NPSOLConFun, NPSOLObjFun, dU_c

	xsol = x0

write(*,*) "calling NPSOL with arguments:"
write(*,*) 'nvar'
write(*,*) nvar
write(*,*) 'g_nclin'
write(*,*) g_nclin
write(*,*) 'g_ncnln'
write(*,*) g_ncnln
write(*,*) 'g_ldA'
write(*,*) g_ldA
write(*,*) 'g_ldJ'
write(*,*) g_ldJ
write(*,*) 'g_ldR'
write(*,*) g_ldR
write(*,*) 'g_linA'
write(*,*) g_linA, size(g_linA)
write(*,*) 'g_bl'
write(*,*) g_bl, size(g_bl)
write(*,*) 'g_bu'
write(*,*) g_bu, size(g_bu)
write(*,*) 'inform'
write(*,*) inform
write(*,*) 'iter'
write(*,*) iter
write(*,*) 'g_istate'
write(*,*) g_istate
write(*,*) 'g_Con'
write(*,*) g_Con
write(*,*) 'g_cJac'
write(*,*) g_cJac
write(*,*) 'g_clambda'
write(*,*) g_clambda
write(*,*) 'fvalue'
write(*,*) fvalue
write(*,*) 'g_grad'
write(*,*) g_grad
write(*,*) 'g_Hess'
write(*,*) g_Hess
write(*,*) 'xsol'
write(*,*) xsol
write(*,*) 'g_iw'
write(*,*) g_iw
write(*,*) 'g_leniw'
write(*,*) g_leniw
write(*,*) 'g_rw'
write(*,*) g_rw
write(*,*) 'g_lenrw'
write(*,*) g_lenrw

	call NPSOL(nvar,g_nclin,g_ncnln,g_ldA,g_ldJ,g_ldR,g_linA,g_bl,g_bu,NPSOLConFun,NPSOLObjFun,inform,& 
		iter,g_istate,g_Con,g_cJac,g_clambda,fvalue,g_grad,g_Hess,xsol,g_iw,g_leniw,g_rw,g_lenrw)

write(*,*) "called NPSOL and now have the following values:"
write(*,*) 'nvar'
write(*,*) nvar
write(*,*) 'g_nclin'
write(*,*) g_nclin
write(*,*) 'g_ncnln'
write(*,*) g_ncnln
write(*,*) 'g_ldA'
write(*,*) g_ldA
write(*,*) 'g_ldJ'
write(*,*) g_ldJ
write(*,*) 'g_ldR'
write(*,*) g_ldR
write(*,*) 'g_linA'
write(*,*) g_linA
write(*,*) 'g_bl'
write(*,*) g_bl
write(*,*) 'g_bu'
write(*,*) g_bu
write(*,*) 'inform'
write(*,*) inform
write(*,*) 'iter'
write(*,*) iter
write(*,*) 'g_istate'
write(*,*) g_istate
write(*,*) 'g_Con'
write(*,*) g_Con
write(*,*) 'g_cJac'
write(*,*) g_cJac
write(*,*) 'g_clambda'
write(*,*) g_clambda
write(*,*) 'fvalue'
write(*,*) fvalue
write(*,*) 'g_grad'
write(*,*) g_grad
write(*,*) 'g_Hess'
write(*,*) g_Hess
write(*,*) 'xsol'
write(*,*) xsol
write(*,*) 'g_iw'
write(*,*) g_iw
write(*,*) 'g_leniw'
write(*,*) g_leniw
write(*,*) 'g_rw'
write(*,*) g_rw
write(*,*) 'g_lenrw'
write(*,*) g_lenrw

	g_NumOptimIter = g_NumOptimIter + iter

	select case (inform)
		case (9)
			g_WorkerInform = WORKER_NPSOL_PARAMETER_FAIL 
			print *, "An input parameter was invalid"
		case (7) 		
			g_WorkerInform = WORKER_NPSOL_DERIVATIVE_FAIL 
			print *, "The function derivatives returned by funcon or fcn is incorrect"
!		case (6)
!			g_WorkerInform = WORKER_NPSOL_1st_ORDER_FAIL
!			print *, "x does not satisfy the 1st-order optimality conditions"
		case (4)
			g_WorkerInform = WORKER_NPSOL_ITER_LIMIT_FAIL
			print *, "The Major iteration limit was reached"
		case (3)
			g_WorkerInform = WORKER_NPSOL_NONLIN_CON_FAIL 
			print *, "The nonlinear constraints could not be satisfied"
		case (2)
			g_WorkerInform = WORKER_NPSOL_LIN_CON_FAIL 
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
		g_WorkerInform = WORKER_NPSOL_FUN_NAN_FAIL
	end if

end subroutine OptimMethod


! the objective function for NPSOL method
! Given today's capital stock k and estimated value function, 
! compute the value dependent on Labor and tomorrow's capital stock 
! objective function for npsol
subroutine NPSOLObjFun(mode,n,x,v,grad,nstate)

	use Optim_global_vars
	IMPLICIT NONE
	
	INTEGER :: mode,n,nstate
	double precision :: v
	double precision :: x(n), grad(n)	
	double precision :: tmp(MAX_DIM)
	integer :: i, j, d
	double precision, external :: U, ExpectValue
	external dU, ValuePlusFun, dValuePlusFun, InitValueFun, dInitValueFun

	v = U(DIM,x(1:DIM),x((DIM+1):(2*DIM)),x((2*DIM+1):(3*DIM)))
	if (g_derOpt == DERIVATIVE_OPT) then
		call dU(DIM,x(1:DIM),x((DIM+1):(2*DIM)),x((2*DIM+1):(3*DIM)), grad)
	end if

	do j = 1, g_TotalNKEps
		do i = 1, DIM
			tmp(i) = x(2*DIM+i) + g_KEps(g_KEpsIndT(j,i),i)
		end do

		if (g_derOpt == DERIVATIVE_OPT) then
			if (g_isInitStep == 0) then
				call dValuePlusFun(g_TotalThetaN, DIM, tmp, g_vPlus(1:g_TotalThetaN,j), g_dvPlus(1:g_TotalThetaN,j,1:DIM))
			else
				call dInitValueFun(DIM, tmp, g_vPlus(1,j), g_dvPlus(1,j,1:DIM))
				do i = 2, g_TotalThetaN
					g_vPlus(i,j) = g_vPlus(1,j)
					do d = 1, DIM
						g_dvPlus(i,j,d) = g_dvPlus(1,j,d)
					end do
				end do
			end if
		else
			if (g_isInitStep == 0) then
				call ValuePlusFun(g_TotalThetaN, DIM, tmp, g_vPlus(1:g_TotalThetaN,j))
			else
				call InitValueFun(DIM, tmp, g_vPlus(1,j))
				do i = 2, g_TotalThetaN
					g_vPlus(i,j) = g_vPlus(1,j)
				end do
			end if
		end if
	end do

	v = v + g_beta*ExpectValue(g_TotalThetaN,g_TotalNKEps,g_vPlus)
	if (g_derOpt == DERIVATIVE_OPT) then
		do i = 1, DIM
			grad(2*DIM+i) = grad(2*DIM+i) + g_beta*ExpectValue(g_TotalThetaN,g_TotalNKEps,g_dvPlus(:,:,i))
		end do	
	end if

	! NPSOL solves the minimal problem, but we are dealing with the maximal problem
	! so we change the objective function f(x) to -f(x)
	v = - v
	if (g_derOpt == DERIVATIVE_OPT) then
		grad = -grad
	end if

end subroutine NPSOLObjFun


! nonlinear constraint function for npsol: g_nlnL <= NPSOLConFun(x) <= g_nlnU
subroutine NPSOLConFun(mode,ncnln,n,ldJ,needc,x,cc,cJac,nstate)

	use Optim_global_vars
	IMPLICIT NONE
	
	INTEGER :: mode,ncnln,n,ldJ,nstate, i
	INTEGER :: needc(ncnln)
	double precision :: x(n), cc(ncnln), cJac(ldJ,*)
	double precision, external :: F, dF_L

	if (g_derOpt == DERIVATIVE_OPT) then
		cJac(1:ldJ,1:n) = 0
		do i = 1, DIM
			cc(i) = F(g_todayK(i),x(DIM+i),g_currthetas(i)) - x(i) - x(2*DIM+i)
			cJac(i,i) = -1
			cJac(i,DIM+i) = dF_L(g_todayK(i),x(DIM+i),g_currthetas(i))
			cJac(i,2*DIM+i) = -1
		end do
	else
		do i = 1, DIM
			cc(i) = F(g_todayK(i),x(DIM+i),g_currthetas(i)) - x(i) - x(2*DIM+i)
		end do
	end if

end subroutine NPSOLConFun


subroutine EndOptim()

	use Optim_global_vars

	deallocate(g_linA)
	deallocate(g_linL, g_linU)
	deallocate(g_nlnL, g_nlnU)
	deallocate(g_istate)
	deallocate(g_bl, g_bu, g_clambda) 
	deallocate(g_grad, g_cJac, g_Hess, g_Con)
	deallocate(g_iw, g_rw) 

end subroutine EndOptim


function ExpectValue(n,m,v)

	use Optim_global_vars
	IMPLICIT NONE

	integer, intent(IN) :: n, m
	double precision, intent(IN) :: v(n,m)
	double precision :: ExpectValue
	integer :: k, j

	ExpectValue = 0
	do k = 1, n
		do j = 1, m
			ExpectValue = ExpectValue + v(k,j)*g_TotalProb(k,j)
		end do
	end do

end function ExpectValue



subroutine ValuePlusFun(n,ddim,x,vplus)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE

	integer, intent(IN) :: n, ddim		
	double precision, INTENT(IN) :: x(ddim)
	double precision, INTENT(OUT) :: vplus(n)
	integer :: i

	call ChebyBasisFunc(g_NumBasis,ddim,x(1:ddim),g_kmins(1:ddim),g_kmaxs(1:ddim),g_deg,g_phi)
	do i = 1, n
		vplus(i) = dot_product(g_ChebyCoefs(1:g_NumBasis,i),g_phi)
	end do

end subroutine ValuePlusFun


subroutine dValuePlusFun(n,ddim,x,vplus,dvplus)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE

	integer, intent(IN) :: n, ddim		
	double precision, INTENT(IN) :: x(ddim)
	double precision, INTENT(OUT) :: vplus(n), dvplus(n,ddim)
	integer :: i, d

	call ChebyBasisFuncD(g_NumBasis, ddim,x(1:ddim),g_kmins(1:ddim),g_kmaxs(1:ddim),g_deg,g_phi,g_dphi)
	do i = 1, n
		vplus(i) = dot_product(g_ChebyCoefs(1:g_NumBasis,i),g_phi)
		do d = 1, ddim
			dvplus(i,d) = dot_product(g_ChebyCoefs(1:g_NumBasis,i), g_dphi(1:g_NumBasis,d))
		end do
	end do

end subroutine dValuePlusFun


! initial guess of value function
subroutine InitValueFun(ddim,x,vplus)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE

	integer, intent(IN) :: ddim	
	double precision, INTENT(IN) :: x(ddim)
	double precision, INTENT(OUT) :: vplus
	double precision :: Lss(ddim), c(ddim), theta0
	integer :: i

	double precision, external :: U, F

	if (g_initValueOpt == ZERO_INIT_VALUEFUN) then
		vplus = 0
	else
		if (g_initValueOpt == STEADY_STATE_INIT_VALUEFUN) then	
		        theta0 = 1.0
			do i = 1, ddim
				Lss(i) = 1.0 ! when theta=1, the steady state is kss=Lss=1 by choosing g_BB = (1-g_alpha) * g_A^(1-g_gamma)
				c(i) = F(x(i),Lss(i),theta0) - x(i)
			end do

			vplus = U(ddim,c,Lss,x)/(1-g_beta)
		elseif (g_initValueOpt == QUADRATIC_INIT_VALUEFUN) then
			vplus = 0
			do i = 1, ddim
				vplus = vplus - 100*(x(i)-1.0)**2
			end do
		else if (g_initValueOpt == UTILITY_INIT_VALUEFUN) then
			vplus = 0
			do i = 1, ddim
				vplus = vplus + myPower(x(i),1-g_gamma)/(1-g_gamma)
			end do
		else
			print *, "The option for an initial guess of value function is not known"			
			g_WorkerInform = WORKER_NOT_DONE_FEATURE
		end if
	end if

end subroutine InitValueFun


! initial guess of value function and its derivative
subroutine dInitValueFun(ddim,x,vplus,dvplus)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE
	
	integer, intent(IN) :: ddim
	double precision, INTENT(IN) :: x(ddim)
	double precision, INTENT(OUT) :: vplus, dvplus(ddim)
	double precision :: Lss(ddim), c(ddim), theta0, dc(ddim), dUv(3*ddim)
	integer :: i

	double precision, external :: U, F, dF_k
	external dU

	if (g_initValueOpt == ZERO_INIT_VALUEFUN) then
		vplus = 0
		dvplus = 0
	else
		if (g_initValueOpt == STEADY_STATE_INIT_VALUEFUN) then	
		        theta0 = 1.0
			do i = 1, ddim
				Lss(i) = 1.0 ! when theta=1, the steady state is kss=Lss=1 by choosing g_BB = (1-g_alpha) * g_A^(1-g_gamma)
				c(i) = F(x(i),Lss(i),theta0) - x(i)
				dc(i) = dF_k(x(i),Lss(i),theta0) - 1
			end do

			vplus = U(ddim,c,Lss,x)/(1-g_beta)
			call dU(ddim,c,Lss,x,dUv)
			do i = 1, ddim
				dvplus(i) = (dUv(i)*dc(i)+dUv(2*ddim+i)) /(1-g_beta)
			end do
		elseif (g_initValueOpt == QUADRATIC_INIT_VALUEFUN) then
			vplus = 0
			do i = 1, ddim
				vplus = vplus - 100*(x(i)-1.0)**2
				dvplus(i) = -200*(x(i)-1.0)
			end do
		else if (g_initValueOpt == UTILITY_INIT_VALUEFUN) then
			vplus = 0
			do i = 1, ddim
				vplus = vplus + myPower(x(i),1-g_gamma)/(1-g_gamma)
				dvplus(i) = myPower(x(i),-g_gamma)
			end do
		else
			print *, "The option for an initial guess of value function is not known"			
			g_WorkerInform = WORKER_NOT_DONE_FEATURE
		end if
	end if

end subroutine dInitValueFun


! Utility function
function U(n,c,L,kPlus)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE
	
	integer, intent(IN) :: n
	double precision, intent(IN) :: c(n), L(n), kPlus(n)
	double precision :: U
	integer :: i, j

	U = 0
	do i = 1, n-1
		do j = i+1, n
			U = U + (kPlus(i)-kPlus(j))**2
		end do
	end do
	U = g_mu*U 

	do i = 1, n
		U = U + myPower(c(i),1-g_gamma)/(1-g_gamma) - g_BB*myPower(L(i),1+g_eta)/(1+g_eta)
	end do

end function U



! gradient of Utility function
subroutine dU(n,c,L,kPlus,dv)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE
	
	integer, intent(IN) :: n
	double precision, intent(IN) :: c(n), L(n), kPlus(n)
	double precision :: dv(3*n)
	integer :: i, j

	dv = 0
	do i = 1, n-1
		do j = i+1, n
			dv(2*n+i) = dv(2*n+i) + 2*(kPlus(i)-kPlus(j))
			dv(2*n+j) = dv(2*n+j) - 2*(kPlus(i)-kPlus(j))
		end do
	end do
	dv = g_mu*dv

	do i = 1, n
!		dv(i) = dv(i) + myPower(c(i),-g_gamma)
!		dv(n+i) = dv(n+i) - g_BB*myPower(L(i),g_eta)
		dv(i) = dv(i) + dmyPower(c(i),1-g_gamma)/(1-g_gamma)
		dv(n+i) = dv(n+i) - g_BB*dmyPower(L(i),1+g_eta)/(1+g_eta)
	end do

end subroutine dU


! gradient of Utility function
subroutine dU_c(n,c,L,kPlus,dv)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE
	
	integer, intent(IN) :: n
	double precision, intent(IN) :: c(n), L(n), kPlus(n)
	double precision :: dv(n)
	integer :: i

	do i = 1, n
!		dv(i) = myPower(c(i),-g_gamma)
		dv(i) = dmyPower(c(i),1-g_gamma)/(1-g_gamma)
	end do

end subroutine dU_c


! production function 
function F(k,L,theta)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE

	double precision, INTENT(IN) :: k, L, theta
	double precision :: F

	F = k + theta * g_A * myPower(k,g_alpha)*myPower(L,1-g_alpha)

end function F

! partial derivative of production function 
function dF_k(k,L,theta)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE

	double precision, INTENT(IN) :: k, L, theta
	double precision :: dF_k

!	dF_k = 1 + theta * g_A * g_alpha * myPower(k,g_alpha-1)*myPower(L,1-g_alpha)
	dF_k = 1 + theta * g_A * dmyPower(k,g_alpha)*myPower(L,1-g_alpha)

end function dF_k


! partial derivative of production function 
function dF_L(k,L,theta)

	use Optim_global_vars
	use Common_Funs
	IMPLICIT NONE

	double precision, INTENT(IN) :: k, L, theta
	double precision :: dF_L

!	dF_L = theta * g_A * myPower(k,g_alpha)* (1-g_alpha)*myPower(L,-g_alpha)
	dF_L = theta * g_A * myPower(k,g_alpha)* dmyPower(L,1-g_alpha)

end function dF_L

