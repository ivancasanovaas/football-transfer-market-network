
!-------------------------------------------------------------------
! MODULE: Network Structural Properties
! AUTHOR: Ivan Casanovas Rodríguez
! DATE: July, 2024
! DESCRIPTION: This module contains functions & subroutines to read 
!			   an undirected and unweighted network (given in edge 
!			   list format), and calculate its structural properties
!			   (number of edges and nodes, list of degrees, list of 
!			   neighbors and pointer, degree distributions, average
!			   nearest neighbors degree and clustering coefficient).
!-------------------------------------------------------------------


MODULE STRUCTURAL_PROPERTIES

IMPLICIT NONE

CONTAINS ! define all the functions & subroutines


!-------------------------------------------------------------------
!                        NUMBER OF EDGES (E)                        
!-------------------------------------------------------------------
! FUNCTION: num_edges
! PURPOSE: Computes the number of edges given the edge list file
! INPUT: 
!   - edgelist_filepath: file name containing the edge list
! OUTPUTS:
!   - N: number of nodes
!-------------------------------------------------------------------

FUNCTION num_edges(edgelist_filepath) result(E)

	character(len=*), intent(in) 			:: edgelist_filepath	! character (*)
	character(len=256) 						:: line					! character (256)
	integer*8								:: ios					! auxiliar int.
	integer*8								:: E					! int. variables

	E = 0

	OPEN(0, file=trim(adjustl(edgelist_filepath)), status='old', action='read', iostat=ios)
		DO
			read(0, '(A)', iostat=ios) line
			if (ios /= 0) then
				if (ios == -1) then
					exit  ! end of the file
				else
					print *, 'ERROR reading the file.'
					stop
				end if
			end if
			E = E + 1
		END DO
	CLOSE(0)

END FUNCTION num_edges


!-------------------------------------------------------------------
!                        NUMBER OF NODES (N)                        
!-------------------------------------------------------------------
! FUNCTION: num_nodes
! PURPOSE: Computes the number of nodes given the edge list file
! INPUT: 
!   - edgelist_filepath: file name containing the edge list
! OUTPUTS:
!   - N: number of nodes
!-------------------------------------------------------------------

FUNCTION num_nodes(edgelist_filepath) result(N)

	character(len=*), intent(in)			:: edgelist_filepath	! character (*)
	integer*8 								:: i,j,ios 				! auxiliar int.
	integer*8			 					:: N					! int. variables

	N = 0

	OPEN(0, file=trim(adjustl(edgelist_filepath)), status='old', action='read', iostat=ios)
		DO 
			read(0, *, iostat=ios) i,j
			if (ios /= 0) then
				if (ios == -1) then
					exit  ! end of the file
				else
					print *, 'ERROR reading the file.'
					stop
				end if
			end if
			if (N < i) then
				N = i
			end if
			if (N < j) then
				N = j
			end if
		END DO
	CLOSE(0)

END FUNCTION num_nodes


!-------------------------------------------------------------------
!                        EDGE LIST (E_list)
!-------------------------------------------------------------------
! FUNCTION: edges_list
! PURPOSE: Reads and stores the edge list of the network, given the
!		   edge list file
! INPUTS: 
!   - edgelist_filepath: file name containing the edge list
!	- E: number of edges
! OUTPUT:
!   - E_list: integer array variable 2D containing the edge list
!-------------------------------------------------------------------

FUNCTION edges_list(edgelist_filepath,E) result(E_list)

	character(len=*), intent(in)			:: edgelist_filepath	! character (*)
	integer*8 								:: i,j,ij,ios			! auxiliar int.
	integer*8 								:: E					! int. variables
	integer*8, dimension(:,:), allocatable	:: E_list(:,:)			! int. arrays (2D)

	allocate(E_list(E,2))
	E_list = 0

    OPEN(0, file=trim(adjustl(edgelist_filepath)), status='old', action='read', iostat=ios)
		DO	ij = 1, E
			read(0, *, iostat=ios) i,j
			E_list(ij,1) = i
			E_list(ij,2) = j
		END DO
    CLOSE(0)

END FUNCTION edges_list


!-------------------------------------------------------------------
!                        LIST OF DEGREES (D)
!-------------------------------------------------------------------
! FUNCTION: node_degree
! PURPOSE: Computes the list of degrees of the network
! INPUT: 
!   - E_list: integer array variable 2D containing the edge list
! OUTPUT:
!   - D: list of degrees for each node
!-------------------------------------------------------------------

FUNCTION node_degree(E_list) result(D)

	integer*8 								:: i,j,ij				! auxiliar int.
	integer*8								:: E,N					! int. variables
	integer*8, dimension(:,:), intent(in)	:: E_list				! int. arrays (2D)
	integer*8, dimension(:), allocatable 	:: D					! int. arrays (1D)

	E = size(E_list, 1)
	N = maxval(E_list)
	allocate(D(N))
	D = 0

	DO ij = 1, E
		i = E_list(ij,1)
		j = E_list(ij,2)
		if (i /= j) then
			D(i) = D(i) + 1
			D(j) = D(j) + 1
		end if
	END DO

END FUNCTION node_degree


!-------------------------------------------------------------------
!                       LIST OF NEIGHBORS (V)
!-------------------------------------------------------------------
! SUBROUTINE: node_neighbors
! PURPOSE: Computes the list of degrees of the network
! INPUT: 
!   - E_list: integer array variable 2D containing the edge list
! OUTPUTS:
!   - V: list of neighbors for each node
!	- ptr: list with the positions, in V, of the first and last
!			   neighbors for each node
!-------------------------------------------------------------------

SUBROUTINE node_neighbors(E_list,D,V,ptr)

	integer*8 								:: i,j,ij				! auxiliar int.
	integer*8								:: E,N					! int. variables
	integer*8, dimension(:,:), intent(in)	:: E_list				! int. arrays (2D)
	integer*8, dimension(:), intent(in) 	:: D 					! int. arrays (1D)
	integer*8, dimension(:), allocatable 	:: V,ptr			! int. arrays (1D)

	E = size(E_list, 1)
	N = maxval(E_list)
	allocate(V(2*E))
	allocate(ptr(2*N))
	V = 0
	ptr = 0

	! initialiting the first pointer of each node (frozen numbers)
	ptr(1) = 1
	ptr(N+1) = 1
	DO i = 1, (N-1)
		ptr(i+1) = ptr(i) + D(i)
		ptr(i+1+N) = ptr(i+1)
	ENDDO

	! filling the neighbors in the position of V according to the second pointer of each node
	DO ij = 1, E
		i = E_list(ij,1)
		j = E_list(ij,2)
		if (i /= j) then
			V(ptr(i+N)) = j
			V(ptr(j+N)) = i
			ptr(i+N) =  ptr(i+N) + 1
			ptr(j+N) =  ptr(j+N) + 1
		endif
	ENDDO

END SUBROUTINE node_neighbors


!-------------------------------------------------------------------
!                DEGREE DISTRIBUTIONS (ddd,cdd,ccdd)
!-------------------------------------------------------------------
! SUBROUTINE: degree_distributions
! PURPOSE: Calculates the degree distributions of the network
!		   (direct, cumulative and complementary cumulative)
! INPUT: 
!   - D: list of degrees for each node
! OUTPUTS:
!   - ddd: direct degree distribution
!	- cdd: cumulative degree distribution
!	- ccdd: complementary cumulative degree distribution 
!-------------------------------------------------------------------

SUBROUTINE degree_distributions(D,ddd,cdd,ccdd)

	integer*8 								:: i,j,k				! auxiliar int.
	integer*8			  					:: N					! int. variables
	integer*8, dimension(:), intent(in)		:: D					! int. arrays (1D)
	real*8, dimension(:), allocatable 		:: ddd,cdd,ccdd 		! real arrays (1D)

	N = size(D)
	allocate(ddd(N))
	allocate(cdd(N))
	allocate(ccdd(N))
	ddd = 0.d0
	cdd = 0.d0
	ccdd = 0.d0

	! computing the probabilities of each distribution
	DO i = 1, N
		DO j = 1, N
			if (D(i) == j) then
				ddd(j) = ddd(j) + 1.d0 
			end if
			if (D(i) <= j) then
				cdd(j) = cdd(j) + 1.d0 
			end if
			if (D(i) >= j) then
				ccdd(j) = ccdd(j) + 1.d0 
			end if
		END DO
	END DO

	! normalizing to the highest possible degree (number of nodes)
	DO k = 1, N
		ddd(k) = ddd(k) / N
		cdd(k) = cdd(k) / N
		ccdd(k) = ccdd(k) / N
	END DO

END SUBROUTINE degree_distributions


!-------------------------------------------------------------------
!              AVERAGE NEAREST NEIGHBORS DEGREE (knn)
!-------------------------------------------------------------------
! FUNCTION: average_nearest_neighbors_degree
! PURPOSE: Computes the average nearest neighbors degree
! INPUTS: 
!   - D: list of degrees for each node
!   - V: list of neighbors for each node
!	- ptr: list with the positions, in V, of the first and last
!			   neighbors for each node
! OUTPUT:
!   - knn: list of average nearest neighbors degrees as a function
!			of the degree
!-------------------------------------------------------------------

FUNCTION average_nearest_neighbors_degree(D,V,ptr) result(knn)

	integer*8 								:: i,j,k				! auxiliar int.
	integer*8			  					:: N,N_k,k_sum,k2_sum	! integer variables
	real*8 									:: k_exp,k2_exp,kappa 	! real variables			
	integer*8, dimension(:), intent(in)		:: D,V,ptr 				! int. arrays (1D)
	real*8, dimension(:), allocatable       :: knn 					! real arrays (1D)

	N = size(D)
	if (.not.allocated(knn)) then
        allocate(knn(N))
    end if
	knn = 0.d0

	! computing the normalization factor (kappa)
	k_sum = 0
	k2_sum = 0
	DO i = 1, N
		k_sum = k_sum + D(i)
		k2_sum = k2_sum + D(i)**2
	END DO
	k_exp = dble(k_sum)/N
	k2_exp = dble(k2_sum)/N
	kappa = k2_exp/k_exp

	DO k = 1, N ! for all possible degrees
		N_k = 0 ! nodes with degree k
		DO i = 1, N ! for all nodes
			if (D(i) == k) then ! if node i has degree k
				N_k = N_k + 1
				DO j = ptr(i), ptr(i+N) ! for all neighbors of i
					knn(k) = knn(k) + dble(D(V(j))) ! sum the degree of the neighbor
				END DO
			end if
		END DO
		if (N_k > 0) then
			knn(k) = knn(k) / k / N_k / kappa 
		end if
	END DO

END FUNCTION average_nearest_neighbors_degree


!-------------------------------------------------------------------
!                    CLUSTERING COEFFICIENT (c)
!-------------------------------------------------------------------
! FUNCTION: clustering_coefficient
! PURPOSE: Computes the average clustering coefficient for each k
! INPUTS: 
!   - D: list of degrees for each node
!   - V: list of neighbors for each node
!	- ptr: list with the positions, in V, of the first and last
!			   neighbors for each node
! OUTPUT:
!   - c: list of clustering coefficients as a function of the degree
!-------------------------------------------------------------------

FUNCTION clustering_coefficient(D,V,ptr) result(c)

	integer*8 								:: i,j,k,l,m			! auxiliar int.
	integer*8			  					:: N,N_k				! integer variables
	integer*8, dimension(:), intent(in)		:: D,V,ptr 				! int. arrays (1D)
	real*8, dimension(:), allocatable 		:: c	 				! real arrays (1D)

	N = size(D)
	if (.not.allocated(c)) then
        allocate(c(N))
    end if	
	c = 0.d0

	DO k = 2, N ! for all possible degrees (nodes with k>1)
		N_k = 0 ! nodes with degree k
		DO i = 1, N ! for all nodes
			if (D(i) == k) then ! if node i has degree k
				N_k = N_k + 1
				DO j = ptr(i), ptr(i+N) ! for all neighbors V(j) of node i
					DO l = ptr(V(j)), ptr(V(j)+N) ! for all neighbors V(l) of node V(j)
						DO m = ptr(V(l)), ptr(V(l)+N) ! for all neighbors V(m) of node V(l)
							if (V(m) == i) then ! if node i is a neighbor of node V(m)
								c(k) = c(k) + 1.d0 ! count a 'triangle'
							end if
						END DO
					END DO 
				END DO
			end if
		END DO
		if (N_k > 0) then
			c(k) = c(k) * 2.d0 / (k*(k-1)) / N_k
		end if	
	END DO

END FUNCTION clustering_coefficient


END MODULE STRUCTURAL_PROPERTIES
