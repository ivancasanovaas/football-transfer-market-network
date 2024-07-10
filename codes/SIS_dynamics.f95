
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

IMPLICIT NONE

CONTAINS  ! define all the functions & subroutines

!-------------------------------------------------------------------
!                             SIS MODEL                                                  
!-------------------------------------------------------------------
! FUNCTION: SIS_model
! PURPOSE: THe function applies the SIS model with the Gillepsie al-
!		   gorithm supposing that the infection/recovery process are
!		   Poisson processes. Given an initial number of infected
!		   nodes, it returns the empirical prevalence as a function 
!		   of time.
! INPUTS: 
!	- E: number of edges
!   - N: number of nodes
!   - D: list of degrees for each node
!   - V: list of neighbors for each node
!	- ptr: list with the positions, in V, of the first and last
!			   neighbors for each node
!	- lambda: infection rate
!	- delta: recovery rate
!	- Ni0: initial number of infected nodes
!	- T: total number of events (an event is an infection or a recovery)
! OUTPUT:
!   - rho: empirical prevalence as a function of time of nodes
!-------------------------------------------------------------------

FUNCTION SIS_model(E,N,D,V,ptr,lambda,delta,Ni0,T) result(rho)

	integer*8								:: i,j,k,time				! auxiliar int.
	integer*8, intent(in)					:: E,N						! int. variables
	integer*8								:: T,Ni0,Ni,Ea				! int. variables
	real*8									:: rnd_num					! auxiliar reals
	real*8 									:: lambda,delta,p	 		! real variables
	integer*8, dimension(:), intent(in) 	:: D,V,ptr					! int. arrays (1D)
	integer*8, dimension(:), allocatable 	:: infected_nodes,stat   	! int. arrays (1D)
	integer*8, dimension(:,:), allocatable 	:: E_list,active_links		! int. arrays (2D)
	real*8, dimension(:), allocatable 		:: rho						! real arrays (1D)


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


	!                           PARAMETERS
	!-------------------------------------------------------------------

	allocate(rho(T))				! empirical prevalence
	allocate(infected_nodes(N))		! list with infected nodes
	allocate(stat(N))				! list with the status of nodes (infected-1- or susceptible-0)
	allocate(active_links(E,2))		! list with active links
	
	Ni = Ni0						! initial number of infected nodes
	Ea = 0 							! initial number of active links 
	infected_nodes = 0				! list with infected nodes
	stat = 0						! list with the status of nodes (infected-1- or susceptible-0)
	active_links = 0				! list with active links
	rho = 0							! empirical prevalence


	!                     INITIAL INFECTED NODES
	!-------------------------------------------------------------------

	DO i = 1, Ni0

		! set the intial list with infected nodes and their status
		j = ceiling(r1279()*N)
		infected_nodes(i) = j
		stat(j) = 1

		! set the initial list with active links
		DO k = ptr(j), ptr(j+N) ! for all neighbors V(k) of node j
			if (stat(V(k)) == 0) then
				Ea = Ea + 1
				active_links(Ea,1) = j ! infected node
				active_links(Ea,2) = V(k) ! susceptible neighbor
			end if
		END DO
	END DO


	!                     		    DYNAMICS
	!-------------------------------------------------------------------

	! for T time steps
	DO time = 1, T

		p = Ea*lambda / (Ni*delta + Ea*lambda) ! probability of choosing infection as a type event
		rnd_num = r1279() ! random number

	!                     		   INFECTIONS
	!-------------------------------------------------------------------

		! if the type event is an INFECTION
		IF (rnd_num <= p) THEN

			! choose randomly an active link and infect the susceptible node
			i = active_links(ceiling(r1279()*Ea),2)
			stat(i) = 1

			! update the list of infected nodes
			Ni = Ni + 1
			infected_nodes(Ni) = i

			! update the list of active links, summing the susceptible neighbors of the new infected node
			DO j = ptr(i), ptr(i+N) 
				if (stat(V(j)) == 0) then
					Ea = Ea + 1
					active_links(Ea,1) = i ! infected node
					active_links(Ea,2) = V(j) ! susceptible neighbor
				end if
			END DO

			! update the list of active links, dropping out links with both nodes infected and setting the last ones to these positions
			DO k = 1, Ea
				if (stat(active_links(k,2)) == 1) then
					active_links(k,1) = active_links(Ea,1)
					active_links(k,2) = active_links(Ea,2)
					active_links(Ea,1) = 0
					active_links(Ea,2) = 0
					Ea = Ea - 1
				end if
			END DO
		END IF

	!                     		   RECOVERIES
	!-------------------------------------------------------------------

		! if the type event is a RECOVERY
		IF (rnd_num > p) THEN
			
			! choose randomly one of the infected nodes to recover
			i = infected_nodes(ceiling(r1279()*Ni))
			stat(i) = 0

			! update the list of infected nodes, dropping the recovered one and setting the last one to that position
			j = infected_nodes(Ni)
			stat(j) = 0
			DO k = 1, Ni
				if (infected_nodes(k) == i) then
					infected_nodes(k) = j
				end if
			END DO
			infected_nodes(Ni) = 0
			Ni = Ni - 1

			! update the list of active links, dropping out links with both nodes susceptible and setting the last ones to these positions
			DO k = 1, Ea
				if (stat(active_links(k,1)) == 0) then
					active_links(k,1) = active_links(Ea,1)
					active_links(k,2) = active_links(Ea,2)
					active_links(Ea,1) = 0
					active_links(Ea,2) = 0
					Ea = Ea - 1
				end if
			END DO

			! update the list of active links, summing the infected neighbors of the new susceptible node
			DO j = ptr(i), ptr(i+N) 
				if (stat(V(j)) == 1) then
					Ea = Ea + 1
					active_links(Ea,1) = V(j) ! infected neighbor
					active_links(Ea,2) = i ! susceptible node
				end if
			END DO
		END IF
		rho(time) = dble(Ni) / N
	END DO

END FUNCTION SIS_model



!-------------------------------------------------------------------
!					   RANDOM NUMBER GENERATOR
!-------------------------------------------------------------------
! FUNCTIONS & SUBROUTINES: r1279(),setr1279(iseed), ran2(idum)
! PURPOSE: Generates a number maximally random
! INPUTS: 
!   - iseed: seed selected
! OUTPUT:
!   - random number
!-------------------------------------------------------------------

FUNCTION r1279()

    IMPLICIT NONE
    INCLUDE "r1279block.h"
    REAL    r1279, inv_max
    REAL    INV_MAXINT 
    PARAMETER (INV_MAXINT = 1.0/2147483647.)

    ioffset = iand(ioffset + 1, 2047)
    irand(ioffset) = (irand(index1(ioffset))*irand(index2(ioffset)))
    r1279 = ishft(irand(ioffset), -1) * INV_MAXINT

END FUNCTION

SUBROUTINE setr1279(iseed)

    IMPLICIT	NONE
    INCLUDE "r1279block.h"
    INTEGER	ibit, ispoke, one_bit, iseed, localseed, NBITM1
    PARAMETER (NBITM1 = 31)
	!
	!	Initialize ioffset. This will be increased by (1 mod 2048) for
	!	each random number which is called. 
	!
    ioffset = 0
	!
	!	Set up the two arrays which give locations of the two random
	!	numbers which will be multiplied to get the new random number
	!
    do ispoke = 0, 2047
	index1(ispoke) = iand(ispoke - 1279, 2047)
	index2(ispoke) = iand(ispoke - 418, 2047)
    end do
	!
	!	set up the initial array of 2048 integer random numbers
	!	Each bit is separately initialized using ran2 from numerical recipes
	!
    localseed = -abs(iseed)


	! Matteo:  I have substituted lshift with ishft


    do ispoke = 0, 2047

	irand(ispoke) = 0
	do ibit = 0, NBITM1
	    one_bit = 0
	    if (ran2(localseed) > 0.5) one_bit = 1
	    irand(ispoke) = ior(irand(ispoke), ishft(one_bit, ibit))
	end do
	irand(ispoke) = 2 * irand(ispoke) + 1

    end do

END SUBROUTINE

FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
 11     continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
END FUNCTION


END MODULE SIS_DYNAMICS
