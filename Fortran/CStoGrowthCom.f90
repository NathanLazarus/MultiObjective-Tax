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

! constants and common global variables which are independent of optimization packages
module Common_global_vars

	! dimensions for state variables
	INTEGER :: DIM
	INTEGER, PARAMETER :: MAX_DIM = 20     ! dimension for capitals

	! number of control variables for the optimization tool
	INTEGER :: g_NVar 
	INTEGER, PARAMETER :: MAX_NCtrlVar = 3*MAX_DIM  

	! comparison table nodes
	integer :: g_tablenum 
	INTEGER, PARAMETER :: MAX_TABLENUM = 200
	double precision :: g_kTable(MAX_TABLENUM) ! comparison table nodes
	integer :: g_TotalTableN
	integer, dimension(:,:), allocatable :: g_TableIndT
	double precision, dimension(:), allocatable :: g_fvnew, g_fvold
		
	! interpolation/approximation methods
	integer :: g_approxOpt      ! approximation option
	INTEGER, PARAMETER :: TENSOR_CHEBY_METHOD = 1
	
	! initial guess of value function
	integer :: g_initValueOpt
	INTEGER, PARAMETER :: ZERO_INIT_VALUEFUN = 0 ! zero function 
	INTEGER, PARAMETER :: STEADY_STATE_INIT_VALUEFUN = 1  ! u(c,l)/(1-beta), 
						! where l=lss, and c is chosen such that k^+=k, 
						! by assuming theta=1
	INTEGER, PARAMETER :: QUADRATIC_INIT_VALUEFUN = 2 ! -100*(k-kss)^2
	INTEGER, PARAMETER :: UTILITY_INIT_VALUEFUN = 3 ! k^(1-gamma)/(1-gamma)

	! constants
	double precision, PARAMETER :: PI = 3.14159265358979
	double precision, PARAMETER :: g_CloseZero = 1.0d-12
	double precision, PARAMETER :: g_pInfty = 1.0d+15

	! model parameters
	double precision :: g_beta    ! discount factor
	double precision :: g_gamma   ! used in the definition of utility function
	double precision :: g_eta     ! used in the definition of utility function
	double precision :: g_mu      ! penalty paramters for cross term (k1-k2)^2
	double precision :: g_alpha   ! used in the definition of production function
	double precision :: g_A, g_BB ! used in the definition of production function

	! bound for states and controls
	double precision :: g_LBrate, g_UBrate
	double precision :: g_kmins(MAX_DIM), g_kmaxs(MAX_DIM) ! bounds for capitals
	double precision, DIMENSION(MAX_NCtrlVar) :: g_LB, g_UB ! lower and upper bounds for variables for the optimization tool 

	! nodes of continuous state variables for value function approximation
	integer :: g_N ! number of nodes
	INTEGER :: g_TotalN
	integer :: g_ActiveN, g_startN, g_endN
	double precision, dimension(:), ALLOCATABLE :: g_kNodes ! nodes
	double precision, dimension(:,:), ALLOCATABLE :: g_kNodess
	integer, dimension(:,:), allocatable :: g_NodeIndT
	double precision :: g_todayK(MAX_DIM) ! today's capitals

	! discrete states
	integer :: g_sparsePMopt
	integer :: g_thetanum(MAX_DIM)  ! numbers of possible values of thetas
	integer :: g_thetaind(MAX_DIM) ! indices of today's thetas
	integer :: g_TotalThetaN
	double precision, dimension(:,:), ALLOCATABLE :: g_thetas ! the lists of possible values of thetas
	double precision, dimension(:,:,:), ALLOCATABLE :: g_PM ! probability matrices for thetas
	double precision, dimension(:), allocatable :: g_thetaPrVec
	double precision :: g_currthetas(MAX_DIM)

	! the lists of possible values of capital epsilon and its probability
	INTEGER, PARAMETER :: MAX_NUM_EPSILON = 20
	integer :: g_NumKEps(MAX_DIM)	! number of capital's epsilon
	integer :: g_TotalNKEps
	double precision, dimension(MAX_NUM_EPSILON, MAX_DIM) :: g_KEps, g_KEpsProb 
	integer, dimension(:,:), allocatable :: g_KEpsIndT

	! derivative option for objective and constraint functions
	integer :: g_derOpt
	INTEGER, PARAMETER :: DERIVATIVE_OPT = 1
	double precision, dimension(:,:,:), allocatable :: g_dvPlus

	integer :: g_deg     ! degree of global/local polynomials
	integer :: g_selfind
	integer :: g_NumOptimIter ! number of iterations for optimization 
	INTEGER :: g_NumBasis
	integer :: g_isInitStep
	double precision :: g_err
	double precision, dimension(:), allocatable :: g_x0, g_xsol ! the optimal solution for the first node
	double precision, dimension(:), allocatable :: g_fvalue
	double precision, dimension(:,:), allocatable :: g_vPlus
	double precision, dimension(:,:), allocatable :: g_TotalProb

	! global variables for Chebyshev approximation	
	integer, dimension(:,:), allocatable :: g_CoefIndT
	double precision :: g_extkmins(MAX_DIM), g_extkmaxs(MAX_DIM) ! extended-kmin and extended-kmax for extended 
								! Chebyshev polynomial such that the first and 
								! last Chebyshev nodes are kmin and kmax respectively
	double precision, dimension(:), allocatable :: g_phi
	double precision, dimension(:,:), ALLOCATABLE :: g_ChebyCoefs ! coefficients
	double precision, dimension(:,:), ALLOCATABLE :: g_dphi, g_ChebyT, g_ChebyDT, g_ChebyDDT ! Chebyshev matrix
	double precision, dimension(:,:), ALLOCATABLE :: g_ProdChebyT, g_ProdChebyDT
	double precision, dimension(:), ALLOCATABLE :: g_ChebyCoefs1, g_oldChebyCoefs1 ! coefficients for one theta

	! communication code between Worker and Master
	integer :: g_WorkerInform
	integer, parameter :: WORKER_SUCCESS = 0
	integer, parameter :: WORKER_NPSOL_SUCCESS = 0
	integer, parameter :: WORKER_NPSOL_FUN_NAN_FAIL = 10
	integer, parameter :: WORKER_NPSOL_PARAMETER_FAIL = 19
	integer, parameter :: WORKER_NPSOL_DERIVATIVE_FAIL = 17
	integer, parameter :: WORKER_NPSOL_1st_ORDER_FAIL = 16
	integer, parameter :: WORKER_NPSOL_ITER_LIMIT_FAIL = 14
	integer, parameter :: WORKER_NPSOL_NONLIN_CON_FAIL = 13
	integer, parameter :: WORKER_NPSOL_LIN_CON_FAIL = 12
	integer, parameter :: WORKER_NPSOL_OTHER_FAIL = 20
	integer, parameter :: WORKER_OPTION_FAIL = 2
	integer, parameter :: WORKER_INIT_FILE_OPEN_FAIL = 3
    	integer, parameter :: WORKER_APPROX_FILE_OPEN_FAIL = 4
	integer, parameter :: WORKER_NODES_FILE_OPEN_FAIL = 5
	integer, parameter :: WORKER_INIT_PARAMETER_FAIL = 6
	integer, parameter :: WORKER_NOT_DONE_FEATURE = 7

end module Common_global_vars


! common functions which are independent of model expression and optimization packages
module Common_Funs

	USE Common_global_vars
	implicit none

contains

	function myPower(x,a)

	!****************************************************************************
	! redefine y=x.^a function such that when x is negative 
	! the function still produces a double precision value y
	!****************************************************************************

		double precision, INTENT(IN) :: x, a
		double precision :: myPower
		double precision :: eps, tmp

		if (a >= -1) then
			eps = 1e-6
		else if (a >= -2) then
			eps = 2e-5
		else if (a >= -4) then
			eps = 3e-4
		else if (a >= -4) then
			eps = 1.6e-3
		else if (a >= -5) then
			eps = 4.6e-3
		else if (a >= -6) then
			eps = 1e-2
		else if (a >= -7) then
			eps = 0.018
		else if (a >= -8) then
			eps = 0.028
		else if (a >= -9) then
			eps = 0.04
		else
			eps = 0.1
		end if

		if (x > eps) then
			myPower = x**a
		else
			tmp = (x-eps)/eps
			myPower = (eps**a)*(1 + a*tmp + 0.5*a*(a-1)*tmp*tmp)
		end if
	
	end function myPower


	function dmyPower(x,a)

	!****************************************************************************
	! derivative of myPower: x^a
	!****************************************************************************

		double precision, INTENT(IN) :: x, a
		double precision :: dmyPower
		double precision :: eps, tmp

		if (a >= -1) then
			eps = 1e-6
		else if (a >= -2) then
			eps = 2e-5
		else if (a >= -4) then
			eps = 3e-4
		else if (a >= -4) then
			eps = 1.6e-3
		else if (a >= -5) then
			eps = 4.6e-3
		else if (a >= -6) then
			eps = 1e-2
		else if (a >= -7) then
			eps = 0.018
		else if (a >= -8) then
			eps = 0.028
		else if (a >= -9) then
			eps = 0.04
		else
			eps = 0.1
		end if

		if (x > eps) then
			dmyPower = a*x**(a-1)
		else
			tmp = (x-eps)/eps
			dmyPower = a*(eps**(a-1))*(1 + (a-1)*tmp)
		end if
	
	end function dmyPower


	function Factorial(N)

		integer :: Factorial, N, i

		Factorial = 1
		do i = 2, N
			Factorial = Factorial*i
		end do

	end function Factorial


	function NChooseK(N,K)

		integer :: NChooseK, N, K

		if ((K>N) .OR. (N<=0) .OR. (K<0)) then
			NChooseK = 0
			print *, N, K
			stop "Error: K>N or N<=0 or K<0 in NChooseK"
		elseif ((K==N) .OR. (K==0)) then
			NChooseK = 1
		elseif ((K==N-1) .OR. (K==1)) then
			NChooseK = N
		elseif ((K==N-2) .OR. (K==2)) then
			NChooseK = N*(N-1)/2
		elseif ((K==N-3) .OR. (K==3)) then
			NChooseK = N*(N-1)*(N-2)/6
		elseif ((K==N-4) .OR. (K==4)) then
			NChooseK = N*(N-1)*(N-2)*(N-3)/24
		else
			NChooseK = Factorial(N) / Factorial(N-K) / Factorial(K)
		end if

	end function NChooseK


	function LinfNorm(N,V1,V2)

		INTEGER, INTENT(IN) :: N
		double precision, DIMENSION(N), INTENT(IN) :: V1, V2
		double precision :: LinfNorm
		INTEGER :: k

		LinfNorm = 0
		do k = 1, N
			LinfNorm = max(LinfNorm, abs((V1(k)-V2(k))/(1+abs(V1(k)))))
		end do

	end function LinfNorm


!**************************************************************************
! functions and subroutines for multidimensional handling
!**************************************************************************

	subroutine SubscriptUnif(n,d,num,inds)

	!****************************************************************************
	! Given a number n, return a d-dimensional positive integer 
	! vector inds such that 
	! n = sum((inds(i)-1)*num**(d-i), {i,1,d-1}) + inds(d)
	!****************************************************************************

		integer, intent(IN) :: n, d, num
		integer :: inds(d), prod, i, m

		if (d == 1) then
			inds = n
		elseif (d == 2) then
			inds(1) = ceiling(1.0*n/num) 
			inds(2) = n - (inds(1)-1)*num
		else
			m = n
			do i = 1, d-2
				prod = num**(d-i)
				inds(i) = ceiling(1.0*m/prod) 
				m = m - (inds(i)-1)*prod
			end do
			inds(d-1) = ceiling(1.0*m/num) 
			inds(d) = m - (inds(d-1)-1)*num
		end if

	end subroutine SubscriptUnif


	subroutine Subscript(n,d,nums,inds)

	!****************************************************************************
	! Given a number n, and d-dimensional vectors nums, return a d-dimensional 
	! vector inds such that 
	! n = sum((inds(i)-1)*product(nums((i+1):d)), {i,1,d-1}) + inds(d)
	! where n>=1, 1<=inds(i)<=nums(i) for 1<=i<=d
	!****************************************************************************

		integer, intent(IN) :: n, d
		integer, intent(IN) :: nums(d)
		integer :: inds(d), prod, i , m

		if (d == 1) then
			inds = n
		elseif (d == 2) then
			inds(1) = ceiling(1.0*n/nums(2)) 
			inds(2) = n - (inds(1)-1)*nums(2)			
		else 
			m = n
			do i = 1, d-2
				prod = product(nums((i+1):d))
				inds(i) = ceiling(1.0*m/prod) 
				m = m - (inds(i)-1)*prod
			end do
			inds(d-1) = ceiling(1.0*m/nums(d)) 
			inds(d) = m - (inds(d-1)-1)*nums(d)			
		end if

	end subroutine Subscript


	subroutine CrossProd(n,d,nd,inds,x,v)

	!****************************************************************************
	! Given a n*d matrix x, return a nd-dimensional vector v 
	! such that v(i)=product(x(inds(i,j),j), {j,1,d}) 
	! for i corresponding to inds(i,1:d), i=1,nd, inds(i,j)=1,n
	!****************************************************************************

		integer, intent(IN) :: n, d, nd
		integer, intent(IN) :: inds(nd,d)
		double precision, intent(IN) :: x(n,d)
		double precision, intent(OUT) :: v(nd)
		integer :: i, j

		do i = 1, nd
			v(i) = x(inds(i,1),1)
			do j = 2, d
				v(i) = v(i)*x(inds(i,j),j)
			end do
		end do

	end subroutine CrossProd


!**************************************************************************
! functions and subroutines for Chebyshev polynomial approximation
!**************************************************************************

	! recursive Chebyshev evaluation at x 
	subroutine ChebyEvalRec(z, deg, T)

		INTEGER, INTENT(IN) :: deg
		double precision, INTENT(IN) :: z
		double precision, DIMENSION(deg+1), INTENT(OUT) :: T
		integer :: i

		T(1) = 1
		T(2) = z
		do i = 2, deg
			T(i+1) = 2*z*T(i)-T(i-1)
		end do

	end subroutine ChebyEvalRec


	! recursive Chebyshev evaluation at x 
	subroutine dChebyEvalRec(z, dz, deg, T, dT)

		INTEGER, INTENT(IN) :: deg
		double precision, INTENT(IN) :: z, dz
		double precision, DIMENSION(deg+1), INTENT(OUT) :: T, dT
		integer :: i

		T(1) = 1
		dT(1) = 0
		T(2) = z
		dT(2) = dz
		do i = 2, deg
			T(i+1) = 2*z*T(i)-T(i-1)
			dT(i+1) = 2*dz*T(i)+2*z*dT(i)-dT(i-1)
		end do

	end subroutine dChebyEvalRec


	! recursive Chebyshev evaluation at x 
	subroutine ddChebyEvalRec(z, dz, deg, T, dT, ddT)

		INTEGER, INTENT(IN) :: deg
		double precision, INTENT(IN) :: z, dz
		double precision, DIMENSION(deg+1), INTENT(OUT) :: T, dT, ddT
		integer :: i

		T(1) = 1
		dT(1) = 0
		ddT(1) = 0
		T(2) = z
		dT(2) = dz
		ddT(2) = 0
		do i = 2, deg
			T(i+1) = 2*z*T(i)-T(i-1)
			dT(i+1) = 2*dz*T(i)+2*z*dT(i)-dT(i-1)
			ddT(i+1) = 4*dz*dT(i)+2*z*ddT(i)-ddT(i-1)
		end do

	end subroutine ddChebyEvalRec



	subroutine ChebyCoefs(d, N, nd, inds, W, coefs)

	!****************************************************************************
	!ChebyCoefs finds coefficients for tensor or complete Chebyshev approximation.
	!   ChebyCoefs takes as input data vector W, where W are the function values; 
	!   The input 'd' is dimension, 'N' is number of grids in one dimension, 
	!   'nd' is the number of all coefficients of 
	!   Chebyshev polynomial approximation. The function finds Chebyshev coefficients 
	!   using Algorithm 6.4 on p. 238 of Judd (1998), and returns a  
	!   vector of length-(g_NumBasis) containing the coefficients.
	!****************************************************************************

		INTEGER, INTENT(IN) :: d, N, nd, inds(N**d,d)
		double precision, INTENT(IN) :: W(N**d)  
		double precision, INTENT(OUT) :: coefs(nd) 
		integer :: i, j, k
		double precision :: Numer
		double precision, dimension(:), allocatable :: TT
		double precision, dimension(:,:), allocatable :: ChebyT1

		allocate(TT(N**d), ChebyT1(N,d))

		! Computes Chebyshev coefficients
		do i = 1, nd
			do j = 1, N
				do k = 1, d
					ChebyT1(j,k) = g_ChebyT(j,g_CoefIndT(i,k))
				end do
			end do
			call CrossProd(N,d,N**d,inds,ChebyT1,TT) 

			Numer = dot_product(W,TT)
			k = 0
			do j = 1, d
				if (g_CoefIndT(i,j) == 1) then
					k = k+1
				end if
			end do
			coefs(i) = (2**(d-k))*Numer/(N**d)
		end do

		deallocate(TT, ChebyT1)

	end subroutine ChebyCoefs


	subroutine GetChebyNodes(N,deg,a,b)

		INTEGER, intent(IN) :: N, deg
		double precision, intent(IN) :: a, b
		INTEGER :: i, j
		double precision :: z, eps

		z = -cos(PI/(2*N))
		eps = -(z+1)*(b-a)/2/z
		g_extkmins(1:DIM) = a - eps
		g_extkmaxs(1:DIM) = b + eps

		do i = 1, N
			! compute the N Chebyshev interpolation nodes on [-1,1]
			z = -cos((2*i-1)*PI/(2*N))

			! compute Chebyshev matrix T
			g_ChebyT(i,1) = 1
			g_ChebyT(i,2) = z
			do j = 2, deg
				g_ChebyT(i,j+1) = 2*z*g_ChebyT(i,j)-g_ChebyT(i,j-1)
			end do

			! Adjust the nodes to the [a,b] interval:
			g_kNodes(i) = (z+1)*(g_extkmaxs(1)-g_extkmins(1))/2 + g_extkmins(1)
		end do

	end subroutine GetChebyNodes


	subroutine GetChebyNodesD(N,deg,a,b)

		INTEGER, intent(IN) :: N, deg
		double precision, intent(IN) :: a, b
		INTEGER :: i, j, k, j1, j2
		double precision :: z, eps, dz
		double precision, dimension(:,:), allocatable :: ChebyT1

		z = -cos(PI/(2*N))
		eps = -(z+1)*(b-a)/2/z
		g_extkmins(1:DIM) = a - eps
		g_extkmaxs(1:DIM) = b + eps
		dz = 2/(g_extkmaxs(1)-g_extkmins(1))

		do i = 1, N
			! compute the N Chebyshev interpolation nodes on [-1,1]
			z = -cos((2*i-1)*PI/(2*N))

			! compute Chebyshev matrix T
			g_ChebyT(i,1) = 1
			g_ChebyDT(i,1) = 0
			g_ChebyDDT(i,1) = 0
			g_ChebyT(i,2) = z
			g_ChebyDT(i,2) = dz
			g_ChebyDDT(i,2) = 0
			do j = 2, deg
				g_ChebyT(i,j+1) = 2*z*g_ChebyT(i,j)-g_ChebyT(i,j-1)
				g_ChebyDT(i,j+1) = 2*dz*g_ChebyT(i,j) + 2*z*g_ChebyDT(i,j) - g_ChebyDT(i,j-1)
				g_ChebyDDT(i,j+1) = 4*dz*g_ChebyDT(i,j) + 2*z*g_ChebyDDT(i,j) - g_ChebyDDT(i,j-1) 
			end do

			! Adjust the nodes to the [a,b] interval:
			g_kNodes(i) = (z+1)/dz + g_extkmins(1)
		end do

		allocate(ChebyT1(g_deg+1,DIM))

		do k = 1, g_TotalN
			do i = 1, g_deg+1
				do j = 1, DIM
					ChebyT1(i,j) = g_ChebyT(g_NodeIndT(k,j),i)
				end do
			end do
			call CrossProd(g_deg+1,DIM,g_NumBasis,g_CoefIndT,ChebyT1, g_ProdChebyT(k,:))
		end do

		do k = 1, g_TotalN
			do i = 1, g_NumBasis
				do j1 = 1, DIM	
					g_ProdChebyDT((j1-1)*g_TotalN+k,i) = g_ChebyDT(g_NodeIndT(k,j1),g_CoefIndT(i,j1))
					do j2 = 1, DIM
						if (j2 /= j1) then
							g_ProdChebyDT((j1-1)*g_TotalN+k,i) = g_ProdChebyDT((j1-1)*g_TotalN+k,i)*g_ChebyT(g_NodeIndT(k,j2),g_CoefIndT(i,j2))
						end if
					end do
				end do
			end do
		end do

		deallocate(ChebyT1)

	end subroutine GetChebyNodesD


	subroutine ChebyBasisFunc(n,d,x,a,b,deg,phi)

	!****************************************************************************
	!ChebyBasisFunc evaluates the n-length vector phi of values of all basis functions 
	!at x in [a,b] in d-dimensional space
	!****************************************************************************

		integer, intent(IN) :: n, d, deg
		double precision, intent(IN) :: x(d), a(d), b(d) 
		double precision, intent(OUT) :: phi(n)
		double precision :: x1, zx, Tx(deg+1,d)
		integer :: i

		! x must be in domain [a, b].
		do i = 1, d
			if (x(i) < a(i)) then
				x1 = a(i)
			else if (x(i) > b(i)) then
				x1 = b(i)
			else
				x1 = x(i)
			end if

			! normalizes x
			zx = 2*(x1-g_extkmins(i))/(g_extkmaxs(i)-g_extkmins(i)) - 1

			call ChebyEvalRec(zx, deg, Tx(1:(deg+1),i))
		end do

		call CrossProd(deg+1,d,n,g_CoefIndT,Tx,phi)

	end subroutine ChebyBasisFunc


	subroutine ChebyBasisFuncD(n,d,x,a,b,deg,phi,dphi)

	!****************************************************************************
	!dChebyBasisFunc evaluates the nxd matrix dphi of gradients of all basis functions 
	!at x in [a,b] in d-dimensional space
	!****************************************************************************

		integer, intent(IN) :: n, d, deg
		double precision, intent(IN) :: x(d), a(d), b(d) 
		double precision, intent(OUT) :: phi(n), dphi(n,d)
		double precision :: x1, zx, Tx(deg+1,d), dTx(deg+1,d), TT(deg+1,d)
		integer :: i

		! x must be in domain [a, b].
		do i = 1, d
			if (x(i) < a(i)) then
				x1 = a(i)
			else if (x(i) > b(i)) then
				x1 = b(i)
			else
				x1 = x(i)
			end if

			! normalizes x
			zx = 2*(x1-g_extkmins(i))/(g_extkmaxs(i)-g_extkmins(i)) - 1

			call dChebyEvalRec(zx, 2/(g_extkmaxs(i)-g_extkmins(i)), deg, Tx(1:(deg+1),i), dTx(1:(deg+1),i))
		end do

		call CrossProd(deg+1,d,n,g_CoefIndT,Tx,phi)
		do i = 1, d
			TT = Tx
			TT(1:(deg+1),i) = dTx(1:(deg+1),i)
			call CrossProd(deg+1,d,n,g_CoefIndT,TT,dphi(1:n,i))
		end do

	end subroutine ChebyBasisFuncD

end module Common_Funs