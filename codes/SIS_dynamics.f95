
!-------------------------------------------------------------------
! MODULE: SIS dynamics
! AUTHOR: Ivan Casanovas Rodríguez
! DATE: July, 2024
! DESCRIPTION: This module reads a undirected and unweighted network
!			   given in edge list format and simulates the Susceptible
!			   Infected-Susceptible model with the Gillepsie algorithm.
!			   An initial number of infected nodes and an infection
!			   rate are required as input parameters.
!-------------------------------------------------------------------


MODULE SIS_DYNAMICS

USE RANDOM_NUMBER_GENERATOR

IMPLICIT NONE

CONTAINS  ! define all the functions & subroutines


!-------------------------------------------------------------------
!                             SIS MODEL                                                  
!-------------------------------------------------------------------
! SUBROUTINE: SIS_model
! PURPOSE: THe function applies the SIS model with the Gillepsie al-
!		   gorithm supposing that the infection/recovery process are
!		   Poisson processes. Given an initial number of infected
!		   nodes, it returns the empirical prevalence as a function 
!		   of time.
! INPUTS: 
!	- E: number of edges
!   - N: number of nodes
!	- Ni0: number of initial infected nodes
!   - D: list of degrees for each node
!   - V: list of neighbors for each node
!	- ptr: list with the positions, in V, of the first and last
!			   neighbors for each node
!	- lambda: infection rate
!	- delta: recovery rate
!	- tmax: total time of the simulation
! OUTPUT:
!	- times: times at which events occur
!   - rhos: empirical prevalence as a function of time
!-------------------------------------------------------------------


SUBROUTINE SIS_model(E,N,Ni0,D,V,ptr,lambda,delta,tmax,times_print,rhos_print)

	integer*8, intent(in) 					:: E,N,Ni0            			! int. variables
	integer*8, dimension(:), intent(in) 	:: D,V,ptr					! int. arrays (1D)
	real*8, intent(in) 						:: lambda,delta   			! real variables
	real*8, intent(in) 						:: tmax						! real variables

	integer*8, dimension(:), allocatable 	:: infected_nodes,stat		! int. arrays (1D)
	integer*8, dimension(:,:), allocatable 	:: active_links				! int. arrays (2D)
	real*8, dimension(:), allocatable 		:: times,rhos				! real arrays (1D)
	real*8, dimension(:), allocatable 		:: times_print,rhos_print	! real arrays (1D)

	integer*8 								:: i,j,k,imcs,zero_index	! auxiliar int.
	integer*8 								:: node,node1,node2			! auxiliar int.
	integer*8 								:: Ni,Ea,Nmcs				! auxiliar int.
	real*8 									:: rnd_num,p,pi				! auxiliar reals 
	real*8 									:: t,tau					! auxiliar reals

	! 				RANDOM NUMBER INITIALIZATION (SEED) 				
	!-------------------------------------------------------------------
	integer, dimension(:), allocatable :: state
	integer :: state_size, iseed
	integer :: cput	
	real*8 :: rnd												
	CALL random_seed(size=state_size)
	allocate(state(state_size))
	CALL system_clock(cput)
	state = cput
	CALL random_seed(put=state)
	CALL random_number(rnd)
	iseed = nint(rnd*200.0d0 + 0.5d0)
	CALL setr1279(iseed)

	!                            PARAMETERS                                                       
	!-------------------------------------------------------------------
	Nmcs = 1000000					! number of 'Monte Carlo' steps

	allocate(times(Nmcs))			! time recording for the change of prevalence
	allocate(rhos(Nmcs))			! empirical prevalence evolution
	allocate(infected_nodes(N))		! list with infected nodes
	allocate(stat(N))				! list with the status of nodes (infected-1- or susceptible-0)
	allocate(active_links(E,2))		! list with active links

	Ni = Ni0						! initial number of infected nodes
	Ea = 0 							! initial number of active links 
	infected_nodes = 0				! list with infected nodes
	stat = 0						! list with the status of nodes (infected-1- or susceptible-0)
	active_links = 0				! list with active links
	times = 0 						! event times
	rhos = 0						! empirical prevalence


	!                        INITIAL CONDITIONS                                               
	!-------------------------------------------------------------------

	! setting the initial list with the infected nodes and the status
	i = 0
	DO WHILE (i < Ni0)
		! choose a random node and infect it if possible
		j = ceiling(r1279()*N)
		IF (stat(j) == 0) THEN
			i = i + 1
			infected_nodes(i) = j
			stat(j) = 1
		END IF
	END DO 

	! setting the initial list with the active links
	DO i = 1, Ni0
		j = infected_nodes(i) ! infected node
		DO k = ptr(j), ptr(j+N) ! for all neighbors V(k) of node j
			if (stat(V(k)) == 0) then
				Ea = Ea + 1
				if (Ea > E) then
					print*, 'ERROR: Active links exceeds number of edges (initialization)'
					stop
				end if
				active_links(Ea,1) = j ! infected node
				active_links(Ea,2) = V(k) ! susceptible neighbor
			end if
		END DO		
	END DO

	rhos(1) = dble(Ni0) / N
	times(1) = 0.d0

	!                     		    DYNAMICS
	!-------------------------------------------------------------------
	t = 0.d0
	imcs = 1

	DO WHILE (t <= tmax)

		! probability of an event (infection or recovery) 
		if (Ni == 0) then
			p = Ea*lambda
		else if (Ni == N) then
			p = Ni*delta
		else
			p = Ni*delta + Ea*lambda
		end if
		if (p == 0.d0) exit ! end of the dynamics when no events are possible

		! generate a time for next event (Poisson process)
		tau = -log(r1279()) / p 
		t = t + tau
		imcs = imcs + 1
		
		! probability of choosing infection as a type event
		if (Ni == N) then
			pi = 0.d0
		else
			pi =  Ea * lambda / p
		end if

		rnd_num = r1279() ! random number


		!                     		   INFECTIONS
		!-------------------------------------------------------------------
		! if the type event is an INFECTION
		IF (rnd_num <= pi) THEN

			! choose randomly an active link and infect the susceptible node
			node = active_links(ceiling(r1279()*Ea),2)
			stat(node) = 1

			! update the list of infected nodes
			Ni = Ni + 1
			infected_nodes(Ni) = node

			! update the list of active links, dropping out links with both nodes infected and setting the last ones to these positions
			k = 1
			DO
				if (active_links(k,2) == node) then
					node1 = active_links(Ea,1)
					node2 = active_links(Ea,2)
					active_links(k,1) = node1
					active_links(k,2) = node2
					active_links(Ea,1) = 0
					active_links(Ea,2) = 0
					Ea = Ea - 1
				else
					k = k + 1
				end if
				if (k > Ea) exit
			END DO

			! update the list of active links, summing the susceptible neighbors of the new infected node
			DO j = ptr(node), ptr(node+N) 
				if (stat(V(j)) == 0) then
					Ea = Ea + 1
					if (Ea > E) then
						print*, 'ERROR: Active links exceeds number of edges (infection)'
						stop
					end if
					active_links(Ea,1) = node
					active_links(Ea,2) = V(j)
				end if
			END DO
		
		!                     		   RECOVERIES
		!-------------------------------------------------------------------
		! if the type event is a RECOVERY
		ELSE IF (rnd_num > pi) THEN

			! choose randomly one of the infected nodes to recover
			node = infected_nodes(ceiling(r1279()*Ni))
			stat(node) = 0
			Ni = Ni - 1
			
			! update the list of infected nodes, dropping the recovered one and setting the last one to that position
			j = infected_nodes(Ni + 1)
			infected_nodes(Ni + 1) = 0
			DO k = 1, Ni
				if (infected_nodes(k) == node) then
					infected_nodes(k) = j
					exit
				end if
			END DO

			! update the list of active links, dropping links with both nodes susceptible and setting the last ones to these positions
			k = 1
			DO
				if (active_links(k,1) == node) then
					node1 = active_links(Ea,1)
					node2 = active_links(Ea,2)
					active_links(k,1) = node1
					active_links(k,2) = node2
					active_links(Ea,1) = 0
					active_links(Ea,2) = 0
					Ea = Ea - 1
				else 
					k = k + 1
				end if
				if (k > Ea) exit
			END DO
			
			! update the list of active links, summing the infected neighbors of the new susceptible node
			DO j = ptr(node), ptr(node+N) 
				if (stat(V(j)) == 1) then
					Ea = Ea + 1
					if (Ea > E) then
						print*, 'ERROR: Active links exceeds number of edges (recovery)'
						stop
					endif
					active_links(Ea,1) = V(j) 
					active_links(Ea,2) = node
				end if
			END DO

		END IF

		times(imcs) = t
		rhos(imcs) = dble(Ni) / N
		
	END DO

	!                         CUTTING THE ARRAYS
	!-------------------------------------------------------------------

	! if finding a 0 in the prevalence
    zero_index = -1
    DO imcs = 1, Nmcs - 1
        if ((times(imcs) == 0.d0) .and. (times(imcs+1) == 0.d0)) then
            zero_index = imcs
            exit
        end if
    END DO

	! if not finding a 0 in the prevalence
    if (zero_index == -1) then
        zero_index = Nmcs + 1
    end if

    ! create new arrays with no zeros in prevalence
    allocate(times_print(zero_index-1))
	allocate(rhos_print(zero_index-1))
    times_print = times(1:zero_index-1)
    rhos_print = rhos(1:zero_index-1)

END SUBROUTINE SIS_model


END MODULE SIS_DYNAMICS
