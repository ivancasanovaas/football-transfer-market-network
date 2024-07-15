
!-------------------------------------------------------------------
! MODULE: Equilibrium Models
! AUTHOR: Ivan Casanovas Rodríguez
! DATE: July, 2024
! DESCRIPTION: This module contains two functions to maximally 
!			   randomize the original undirected andunweighted
!			   network with two methods, mantaining the degree
!			   distribution. The random rewiring process (RW) &
!			   the configuration model (CM) are both presented.
!-------------------------------------------------------------------


MODULE EQUILIBRIUM_MODELS

USE RANDOM_NUMBER_GENERATOR

IMPLICIT NONE

CONTAINS  ! define all the functions & subroutines


!-------------------------------------------------------------------
!                   RANDOM REWIRING PROCESS (RW)                     
!-------------------------------------------------------------------
! FUNCTION: random_rewiring_process
! PURPOSE: Selects 2 different links randomly and exchange their 
!		   connections if possible (avoid self-loops & repeated links).
! INPUTS: 
!   - E_list: integer array variable 2D containing the edge list
! OUTPUT:
!   - E_list_rw: integer array variable 2D containing the rewired 
!				 edge list
!-------------------------------------------------------------------

FUNCTION random_rewiring_process(E_list) result(E_list_rw)

	integer*8 									:: E,i,j,k,l,m,iter ! auxiliar int.
	integer*8 									:: link1,link2,link	! int. variables
	integer*8, dimension(:,:), intent(in) 		:: E_list			! int. arrays (2D)
	integer*8, dimension(:,:), allocatable		:: E_list_rw		! int. arrays (2D)
	logical										:: valid_edge		! boolean variables
	
	! RANDOM NUMBER INITIALIZATION (SEED)
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

	E = size(E_list, 1)

	! creating a copy of the edge list
	allocate(E_list_rw(E,2))
	DO i = 1, E
		E_list_rw(i,1) = E_list(i,1)
		E_list_rw(i,2) = E_list(i,2)
	END DO	

	iter = 0
	
	DO 
		valid_edge = .true.

		! 2 random links selected
		link1 = ceiling(r1279()*E)
		link2 = ceiling(r1279()*E)

		! nodes that are in the links
		i = E_list_rw(link1,1)
		j = E_list_rw(link2,1)
		k = E_list_rw(link1,2)
		l = E_list_rw(link2,2)

		! avoid self-connections
		if ((i == l) .or. (j == k)) then
			valid_edge = .false.
		end if

		! avoid multiple connections (the ones that already exist)
		DO link = 1, E
			if ((E_list_rw(link,1) == i) .and. (E_list_rw(link,2) == l)) then
				valid_edge = .false.
			end if
			if ((E_list_rw(link,1) == j) .and. (E_list_rw(link,2) == k)) then
				valid_edge = .false.
			end if
		END DO

		! set the proposed rewiring if possible
		if (valid_edge .eqv. .true.) then
			iter = iter + 1
			E_list_rw(link1,2) = l
			E_list_rw(link2,2) = k
		end if

		! stop when a 2E rewirings have been done
		if (iter == 2*E) then
			exit
		end if
	END DO

END FUNCTION random_rewiring_process


!-------------------------------------------------------------------
!                        CONFIGURATION MODEL (CM)                     
!-------------------------------------------------------------------
! FUNCTION: configuration_model
! PURPOSE: Produces maximally random networks with a preassigned
!		   degree sequence for the nodes.
! INPUTS: 
!   - E_list: integer array variable 2D containing the edge list
!   - D: list of degrees for each node
! OUTPUT:
!   - E_list_cm: integer array variable 2D containing the edge list 
!				 generated from a preassigned degree sequence
!   - D_cm: list of reassigned degrees for each node
!-------------------------------------------------------------------

FUNCTION configuration_model(E_list,D) result(E_list_cm)

	integer*8 									:: i,j,k		    ! auxiliar int.
	integer*8 									:: E,N,iter		 	! int. variables
	integer*8 									:: link,stub1,stub2 ! int. variables
	integer*8, dimension(:), intent(in) 		:: D				! int. arrays (1D)
	integer*8, dimension(:), allocatable 		:: stubs			! int. arrays (1D)
	integer*8, dimension(:,:), intent(in) 		:: E_list			! int. arrays (2D)
	integer*8, dimension(:,:), allocatable		:: E_list_cm		! int. arrays (2D)
	logical										:: valid_edge		! boolean variables
	
	! RANDOM NUMBER INITIALIZATION (SEED)
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

	E = size(E_list, 1)
	N = size(D)

	! creating a copy of the edge list
	allocate(E_list_cm(E,2))
	DO i = 1, E
		E_list_cm(i,1) = E_list(i,1)
		E_list_cm(i,2) = E_list(i,2)
	END DO
	E_list_cm = 0
	
	! creating the node stubs (half-links) list
	allocate(stubs(2*E))
	k = 1
    DO i = 1, N
        DO j = 1, D(i)
            stubs(k) = i
            k = k + 1
        END DO
    END DO

	! creating the edge list from the stubs list
	iter = 0
	DO
		valid_edge = .true.

		! choosing two stubs uniformly at random and propose to connect them
		stub1 = stubs(ceiling(r1279()*2*E))
		stub2 = stubs(ceiling(r1279()*2*E))

		! avoid self-connections
		if (stub1 == stub2) then
			valid_edge = .false.
		end if

		! avoid multiple connections (the ones that already exist)
		DO link = 1, E
			if ((E_list_cm(link,1) == stub1) .and. (E_list_cm(link,2) == stub2)) then
				valid_edge = .false.
			end if
			if ((E_list_cm(link,1) == stub2) .and. (E_list_cm(link,2) == stub1)) then
				valid_edge = .false.
			end if
		END DO

		! set the proposed connection if possible
		if (valid_edge .eqv. .true.) then
			iter = iter + 1
			E_list_cm(iter,2) = stub1
			E_list_cm(iter,2) = stub2
		end if

		! stop when all the E links are set
		if (iter == E) then
			exit
		end if
    END DO

	deallocate(stubs)

END FUNCTION configuration_model



END MODULE EQUILIBRIUM_MODELS
