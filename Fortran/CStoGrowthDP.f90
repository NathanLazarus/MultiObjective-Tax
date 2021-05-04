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

! This program is trying to solve the stochastic optimal growth problem
! with given utility function u(w) using dynamic programming:
! 
! V(k,theta) = maximize EV : u(c,L,k^+)+beta*E[V(k^+ + epsilon,theta^+) | theta] ;
! subject to  Balance {i in 1..d}: k^+[i] = k[i] + f(k[i],L[i],theta) - c[i];
! where beta is the discount factor, d is the dimension,
! theta^+ is a random variable having a serial correlation with theta, 
! epsilon is a random variable with a given distribution, 
! u(c,L,k^+) = sum{i in 1..d} {c[i]^{1-gamma}/(1-gamma)+B*L[i]^(1+eta)/(1+eta)} 
!        + mu * sum{i in 1..(d-1), j in (i+1)..d} {(k^+[i] - k^+[j])^2},  
! and f(k,L,theta)=theta*A*k^alpha*L^{1-alpha}.

PROGRAM WStoDP

	USE Optim_global_vars
	use Common_Funs
	IMPLICIT NONE 

	INTEGER :: i, j
	double precision, dimension(:), allocatable :: x0, xsol
	EXTERNAL OptimInit, OptimMethod

	g_WorkerInform = WORKER_SUCCESS

	! the worker needs to read the initial setup data and the coefficients data file. 
	call ReadInputData()

	allocate(x0(g_NVar), xsol(g_NVar))
	call AllocateNodeTable()

	! worker use the initial data to set optimization package parameters.
	call OptimInit()

	g_NumOptimIter = 0
	x0 = g_x0
	do i = 1, g_ActiveN
		do j = 1, DIM
			g_todayK(j) = g_kNodes(g_NodeIndT(i,j))
		end do

		if (i == 1) then
			call npoptn('Cold Start')
		else
			call npoptn('Warm Start')
		end if

		! workers use optimization packages to get the optimal solution and objective value
		call OptimMethod(g_NVar, x0, xsol, g_fvalue(i)) 

		if (i==1) then
			g_xsol = xsol
		end if
		x0 = xsol
	end do

	call EndOptim()

	if (g_ActiveN == g_TotalN) then
		! Compute the coefficents of the functions
		call ChebyCoefs(DIM,g_N,g_NumBasis,g_NodeIndT,g_fvalue,g_ChebyCoefs1)

		! evaluate the relative change of function values at test nodes
		call EvalRelativeChange()
	end if

	call WriteSolVal()

	! Deallocate spaces 
	call DeallocateSpaces()
	deallocate(x0, xsol)

!*****************
!** Subprograms **
!*****************

contains

	subroutine initialize()

		integer :: i

		g_TotalNKEps = product(g_NumKEps(1:DIM))
		g_TotalN = g_N**DIM
		allocate(g_kNodes(g_N), g_x0(g_NVar), g_xsol(g_NVar), g_thetaPrVec(g_TotalThetaN)) 
		allocate(g_vPlus(g_TotalThetaN,g_TotalNKEps),g_dvPlus(g_TotalThetaN,g_TotalNKEps,DIM))
		allocate(g_KEpsIndT(g_TotalNKEps,DIM))
		allocate(g_TotalProb(g_TotalThetaN, g_TotalNKEps))

		do i = 1, g_TotalNKEps
			call Subscript(i,DIM,g_NumKEps(1:DIM),g_KEpsIndT(i,1:DIM))
		end do

		! Allocate spaces dependent on approximation methods
		call AllocateApproxSpacesInit()

	end subroutine initialize


	subroutine AllocateNodeTable()

		integer :: i

		allocate(g_fvalue(g_ActiveN))
		! the comparison table nodes 
		if (g_ActiveN == g_TotalN) then
			do i = 1, g_tablenum
				g_kTable(i) = g_kmins(1) + (g_kmaxs(1)-g_kmins(1)) * (i-1) / (g_tablenum-1)
			end do
			g_TotalTableN = g_tablenum**DIM
			allocate(g_fvnew(g_TotalTableN), g_fvold(g_TotalTableN))
			allocate(g_TableIndT(g_TotalTableN,DIM))
		end if
		allocate(g_NodeIndT(g_ActiveN,DIM))

		do i = 1, g_ActiveN
			call SubscriptUnif(g_startN+i,DIM,g_N,g_NodeIndT(i,1:DIM))
		end do

		if (g_ActiveN == g_TotalN) then
			do i = 1, g_TotalTableN
				call SubscriptUnif(i,DIM,g_tablenum,g_TableIndT(i,1:DIM))
			end do
		end if

		call GenTotalProb()

	end subroutine AllocateNodeTable


	subroutine GenTotalProb()

		integer :: i, k, inds(MAX_DIM), n
		double precision, dimension(:), allocatable :: KEpsPrVec

		allocate(KEpsPrVec(g_TotalNKEps)) 

		n = maxval(g_NumKEps(1:DIM))
		call CrossProd(n,DIM,g_TotalNKEps,g_KEpsIndT,g_KEpsProb(1:n,1:DIM),KEpsPrVec)

		if (g_sparsePMopt == 0) then
			do k = 1, g_TotalThetaN
				call Subscript(k,DIM,g_thetanum(1:DIM),inds(1:DIM))
				g_thetaPrVec(k) = g_PM(inds(1),g_thetaind(1),1)
				do i = 2, DIM
					g_thetaPrVec(k) = g_thetaPrVec(k)*g_PM(inds(i),g_thetaind(i),i)
				end do			
			end do
		end if

		do k = 1, g_TotalThetaN
			do i = 1, g_TotalNKEps
				g_TotalProb(k,i) = g_thetaPrVec(k)*KEpsPrVec(i)
			end do
		end do
		
		deallocate(KEpsPrVec)

	end subroutine GenTotalProb


	! Allocate spaces dependent on approximation methods and initialize
	subroutine AllocateApproxSpacesInit()

		integer :: i

		g_NumBasis = (g_deg+1)**DIM

		allocate(g_ChebyCoefs1(g_NumBasis),g_dphi(g_NumBasis,DIM))
		allocate(g_ChebyT(g_N, g_deg+1), g_ChebyCoefs(g_NumBasis,g_TotalThetaN))
		allocate(g_phi(g_NumBasis), g_CoefIndT(g_NumBasis,DIM))

		do i = 1, g_NumBasis
			call SubscriptUnif(i,DIM,g_deg+1,g_CoefIndT(i,1:DIM))
		end do
		call GetChebyNodes(g_N,g_deg,g_kmins(1),g_kmaxs(1))

	end subroutine AllocateApproxSpacesInit


	! Deallocate spaces 
	subroutine DeallocateSpaces()

		if (g_ActiveN == g_TotalN) then
			deallocate(g_fvnew, g_fvold, g_TableIndT)
		end if
		deallocate(g_vPlus, g_dvPlus, g_TotalProb, g_KEpsIndT)
		deallocate(g_NodeIndT, g_thetaPrVec)
		deallocate(g_kNodes, g_x0, g_xsol, g_fvalue)

		if (g_sparsePMopt == 0) then
			deallocate(g_thetas, g_PM)
		end if

		deallocate(g_ChebyT, g_ChebyCoefs)
		deallocate(g_ChebyCoefs1, g_phi, g_dphi, g_CoefIndT)

		if ((g_selfind == 0) .AND. (g_ActiveN == g_TotalN)) then
			deallocate(g_oldChebyCoefs1)
		end if

	end subroutine DeallocateSpaces


	! evaluate the relative change of function values at test nodes
	subroutine EvalRelativeChange()

		integer :: i, j
		double precision :: tmp(MAX_DIM)
		external InitValueFun

		do j = 1, g_TotalTableN
			do i = 1, DIM
				tmp(i) = g_kTable(g_TableIndT(j,i))
			end do

			if (g_isInitStep == 1) then
				call InitValueFun(DIM,tmp,g_fvold(j))
			end if

			call ChebyBasisFunc(g_NumBasis,DIM,tmp,g_kmins(1:DIM),g_kmaxs(1:DIM),g_deg,g_phi)
			if (g_isInitStep == 0) then
				if (g_selfind /= 0) then
					g_fvold(j) = dot_product(g_ChebyCoefs(1:g_NumBasis, g_selfind),g_phi)
				else
					g_fvold(j) = dot_product(g_oldChebyCoefs1(1:g_NumBasis),g_phi)
				end if

				g_fvnew(j) = dot_product(g_ChebyCoefs1(1:g_NumBasis),g_phi)
			end if
		end do

		g_err = LinfNorm(g_TotalTableN,g_fvold,g_fvnew)

	end subroutine EvalRelativeChange


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! Communication Subroutines !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine ReadInputData()

		integer :: j, d, maxn
		integer :: OpenStatus, i, k, i2, i3, i4
		character*10 :: cmd_line_arg
		integer :: task_id

		! open the file "ApproxData.txt"
		open (unit = 30, file = "ApproxData.txt", status = "old", &
			  action = "read", position = "rewind", iostat = OpenStatus)
		if (OpenStatus > 0) then
			g_WorkerInform = WORKER_APPROX_FILE_OPEN_FAIL 
			call WriteSolVal()
			print *, " *** Cannot open file ApproxData ***"
			stop 	
		end if
		
		read (30, *) DIM, g_NVar
		if (DIM > MAX_DIM) then
			g_WorkerInform = WORKER_INIT_PARAMETER_FAIL
			call WriteSolVal()
			print *, " *** the dimension is too large ***"
			stop
		end if

		read (30, *) g_approxOpt, g_derOpt 
		if (g_approxOpt /= TENSOR_CHEBY_METHOD) then
			g_WorkerInform = WORKER_NOT_DONE_FEATURE
			call WriteSolVal()
			print *, "This version only supports tensor Chebyshev polynomials"
			stop
		end if

		read (30,*) g_N, g_deg		
		read (30,*) g_beta, g_gamma, g_eta, g_mu, g_alpha

		read (30,*) g_sparsePMopt
		if (g_sparsePMopt == 0) then
			do d = 1, DIM
				read (30, *) g_thetanum(d)
			end do

			maxn = g_thetanum(1)
			do d = 2, DIM
				maxn = max(maxn, g_thetanum(d))
			end do
			allocate(g_thetas(maxn,DIM), g_PM(maxn,maxn,DIM))

			do d = 1, DIM
				do i = 1, g_thetanum(d)
					read (30, *) g_thetas(i,d)
				end do

				do i = 1, g_thetanum(d)
					do j = 1, g_thetanum(d)
						read (30, *) g_PM(i,j,d)
					end do
				end do
			end do
			g_TotalThetaN = product(g_thetanum(1:DIM))
		else
			read (30,*) g_TotalThetaN  
		end if

		read (30, *) g_A, g_BB


		do i = 1, DIM
			read (30, *) g_kmins(i)
		end do
		do i = 1, DIM
			read (30, *) g_kmaxs(i)
		end do

		do i = 1, g_NVar
			read (30, *) g_LB(i)
		end do
		do i = 1, g_NVar
			read (30, *) g_UB(i)
		end do

		do j = 1, DIM
			read (30,*) g_NumKEps(j)
			if (g_NumKEps(j) > MAX_NUM_EPSILON) then
				g_WorkerInform = WORKER_INIT_PARAMETER_FAIL
				call WriteSolVal()
				print *, " *** the number of epsilons is too large ***"
				stop
			end if
		end do

		do j = 1, DIM
			do i = 1, g_NumKEps(j)
				read (30,*) g_KEps(i,j)
			end do

			do i = 1, g_NumKEps(j)
				read (30,*) g_KEpsProb(i,j)
			end do
		end do

		read (30,*) g_tablenum
		if (g_tablenum > MAX_TABLENUM) then
			g_WorkerInform = WORKER_INIT_PARAMETER_FAIL
			call WriteSolVal()
			print *, " *** the number of test grids is too large ***"
			stop
		end if

		read (30,*) g_initValueOpt


		! initialize: set global variables, allocate spaces, etc.
		call initialize()

		! the index of current theta (which may be not in the non-zero PM list when P(i,i)=0), for dense PM, it it task_id
		read (30, *) g_selfind 
		read (30, *) g_isInitStep
		read (30, *) g_startN
		read (30, *) g_endN
		g_ActiveN = g_endN - g_startN + 1

		do i = 1, g_NVar
			read (30, *) g_x0(i)
		end do

		if (g_sparsePMopt == 0) then
			call getarg(1, cmd_line_arg)
			read (cmd_line_arg, *) task_id
			g_selfind = task_id
			call Subscript(task_id, DIM, g_thetanum, g_thetaind(1:DIM))
			do i = 1, DIM
			! 	# read (30, *) g_thetaind(i)
			! 	call getarg(i, cmd_line_arg)
			! 	read (cmd_line_arg, *) g_thetaind(i)
				g_currthetas(i) = g_thetas(g_thetaind(i),i)
			end do
		else
			do i = 1, DIM
				read (30, *) g_currthetas(i)
			end do

			read (30, *) g_TotalThetaN
			do i = 1, g_TotalThetaN
				read (30, *) g_thetaPrVec(i)
			end do
		end if

		if (g_isInitStep == 0) then
			do k = 1, g_TotalThetaN
				do i = 1, g_NumBasis
					read (30, *) g_ChebyCoefs(i,k)
				end do
			end do
		end if

		if ((g_selfind == 0) .AND. (g_ActiveN == g_TotalN)) then
			allocate(g_oldChebyCoefs1(g_NumBasis))
			do k = 1, g_NumBasis
				read (30, *) g_oldChebyCoefs1(k)
			end do
		end if

		! close the file
		close(30)

	end subroutine ReadInputData

	subroutine WriteSolVal()

		integer :: OpenStatus, i, i2, i3, i4
		character*224 :: SolFileName
		character*200 :: char_theta

		write(char_theta, "(100i4.4)") g_thetaind(1:DIM)
		SolFileName = "SolValTheta" // trim(char_theta) // ".txt"
		! open the file SolFileName
		open (unit = 50, file = SolFileName, status = "old", &
			  action = "write", position = "rewind", iostat = OpenStatus)
		if (OpenStatus > 0) then 
			open (unit = 50, file = SolFileName, status = "new", &
				  action = "write", position = "rewind", iostat = OpenStatus)
			if (OpenStatus > 0) then 
				print *," *** Worker cannot open SolVal file ***"	
				STOP
			end if
		end if

		write (50, *) g_WorkerInform
		if (g_WorkerInform == WORKER_SUCCESS) then
			write (50, *) g_NumOptimIter

			if (g_ActiveN == g_TotalN) then
				do i = 1, g_NumBasis
					write (50, *) g_ChebyCoefs1(i)
				end do

				write (50,*) g_err
			else 
				do i = 1, g_ActiveN
					write (50, *) g_fvalue(i)	
				end do
			end if

			do i = 1, g_NVar
				write (50, *) g_xsol(i)
			end do
		end if

		! close the file
		close(50)

	end subroutine WriteSolVal


!********************
!** End of Program **
!********************
end program WStoDP
