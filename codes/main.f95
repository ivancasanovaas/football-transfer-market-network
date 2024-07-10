
!-------------------------------------------------------------------
! PROGRAM: Main Program
! AUTHOR: Ivan Casanovas Rodríguez
! DATE: July, 2024
! DESCRIPTION: This program reads a undirected and unweighted network
!			   given in edge list format and executes the 3 modules:
!			   structural properties, equilibrium models & SIS dynamics.
!			   The different results are written in .txt files located 
!			   in the given output files path.
!-------------------------------------------------------------------


PROGRAM MAIN

USE STRUCTURAL_PROPERTIES
USE EQUILIBRIUM_MODELS
USE SIS_DYNAMICS

IMPLICIT NONE

character(len=256) 						:: edgelist_filepath		! characters (256)
character(len=256) 						:: output_filespath			! characters (256)
character(len=256)						:: dir,fname_SIS,a,b		! characters (256)
logical 								:: istat					! boolean variables
integer*8 								:: i,j,k,Niter			 	! auxiliar int.
integer*8 								:: E,N,k_sum,k2_sum			! int. variables (properties)
integer*8 								:: Ni0,T,M					! int. variables (SIS)
real*8 									:: k_exp,k2_exp				! real variables (properties)
real*8									:: lambda,delta 			! real variables (SIS)
integer*8, dimension(:), allocatable 	:: D,V,ptr					! int. 1D arrays (properties)
integer*8, dimension(:,:), allocatable 	:: E_list,E_list_em			! int. 2D arrays (properties)
real*8, dimension(:), allocatable 		:: ddd,cdd,ccdd 			! real 1D arrays (properties)
real*8, dimension(:), allocatable 		:: knn,c,knn_em,c_em		! real 1D arrays (properties)
real*8, dimension(:), allocatable 		:: rho						! real 1D arrays (SIS)


!-------------------------------------------------------------------
!           READING THE INPUT PATHS & CREATE DIRECTORIES                                                                    
!-------------------------------------------------------------------

CALL read_paths(edgelist_filepath,output_filespath)

CALL create_directory(trim(adjustl(output_filespath)) // "/properties")
CALL create_directory(trim(adjustl(output_filespath)) // "/SIS")
CALL create_directory(trim(adjustl(output_filespath)) // "/SIS/SIS-life_time")
CALL create_directory(trim(adjustl(output_filespath)) // "/SIS/SIS-rho_st")
CALL create_directory(trim(adjustl(output_filespath)) // "/SIS/SIS-rho_t")



!-------------------------------------------------------------------
!                      STRUCTURAL PROPERTIES
!-------------------------------------------------------------------

!                  NUMBER OF EDGES (E) & NODES (N)
!-------------------------------------------------------------------

E = num_edges(edgelist_filepath)
N = num_nodes(edgelist_filepath)
print *, "Number of edges:", E
print *, "Number of nodes:", N


!                        EDGE LIST (E_list)
!-------------------------------------------------------------------

E_list = edges_list(edgelist_filepath,E)

! writing the edge list in a file
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/edge_list.txt")
	DO i = 1, E
		write(0,fmt='(i6,5X,i6)') E_list(i,1), E_list(i,2)
	END DO
CLOSE(0)


!                        LIST OF DEGREES (D)
!-------------------------------------------------------------------

D = node_degree(E_list)

! writing the node degree in a file
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/degrees.txt")
write(0,fmt='(A6,5X,A6)') '# node', 'degree'
	DO i = 1, N
		write(0,fmt='(i6,5X,i6)') i, D(i)
	END DO
CLOSE(0)

! computing the average node degree
k_sum = 0
k2_sum = 0
DO i = 1, N
	k_sum = k_sum + D(i)
	k2_sum = k2_sum + D(i)**2
END DO
k_exp = dble(k_sum)/N
k2_exp = dble(k2_sum)/N
print*, '<k> = ', k_exp
print*, '<k2>/<k> = ', k2_exp/k_exp



!                       LIST OF NEIGHBORS (V)
!-------------------------------------------------------------------

CALL node_neighbors(E_list,D,V,ptr)

! writing the node neighbors in a file
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/neighbors.txt")
	DO i = 1, 2*E
		write(0,fmt='(i6)') V(i)
	ENDDO
CLOSE(0)

! writing the node pointers in a file
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/pointer.txt")
write(0,fmt='(A6,4X,A7)') '# node', 'pointer'
	DO i = 1, 2*N
		if (i <= N) then
			write(0,fmt='(i6,5X,i6)') i, ptr(i)
		else 
			write(0,fmt='(i6,5X,i6)') i-N, ptr(i)
		end if
	ENDDO
CLOSE(0)


!                DEGREE DISTRIBUTIONS (ddd,cdd,ccdd)
!-------------------------------------------------------------------

CALL degree_distributions(D,ddd,cdd,ccdd)

! opening files to store the distributions
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/ddd.txt")
OPEN(1,file=trim(adjustl(output_filespath))//"/properties/cdd.txt")
OPEN(2,file=trim(adjustl(output_filespath))//"/properties/ccdd.txt")


write(0,fmt='(A8,16X,A3)') '# degree','ddd'
write(1,fmt='(A8,16X,A3)') '# degree','cdd'
write(2,fmt='(A8,15X,A4)') '# degree','ccdd'

	! writing the distributions in a file
	DO k = 1, N
		write(0,fmt='(i8,5X,f14.8)') k, ddd(k)
		write(1,fmt='(i8,5X,f14.8)') k, cdd(k)
		if (ccdd(k) /= ccdd(k+1)) then
			write(2,fmt='(i8,5X,f14.8)') k, ccdd(k)
		end if
	END DO

CLOSE(0)
CLOSE(1)
CLOSE(2)


!              AVERAGE NEAREST NEIGHBORS DEGREE (knn)
!-------------------------------------------------------------------

knn = average_nearest_neighbors_degree(D,V,ptr)

! writing the average nearest neighbors degree in a file
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/knn.txt")
write(0,fmt='(A8,16X,A3)') '# degree','knn'
	DO k = 1, N
		write(0,fmt='(i8,5X,f14.8)') k, knn(k)
	ENDDO
CLOSE(0)


!                    CLUSTERING COEFFICIENT (c)
!-------------------------------------------------------------------

 c = clustering_coefficient(D,V,ptr)

! writing the average clustering coefficiet in a file
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/c.txt")
write(0,fmt='(A8,18X,A1)') '# degree','c'
	DO k = 2, N
		write(0,fmt='(i8,5X,f14.8)') k, c(k)
	ENDDO
CLOSE(0)


!-------------------------------------------------------------------
!                    	  EQUILIBRIUM MODELS
!-------------------------------------------------------------------

!                   RANDOM REWIRING PROCESS (RW)                     
!-------------------------------------------------------------------

allocate(knn_em(N))
allocate(c_em(N))

! run 'Niter' RW networks and calculate its properties (knn,c)
Niter = 100
DO i = 1, Niter

	deallocate(V,ptr,knn,c)

	E_list_em = random_rewiring_process(E_list)
	CALL node_neighbors(E_list_em,D,V,ptr)
	knn = average_nearest_neighbors_degree(D,V,ptr)
 	c = clustering_coefficient(D,V,ptr)
	
	DO k = 1, N
		knn_em(k) = knn_em(k) + knn(k)  
		c_em(k) = c_em(k) + c(k)  
	END DO
	
	deallocate(E_list_em)

END DO

! average the properties over the 'Niter' RW models
DO k = 1, N
	knn_em(k) = knn_em(k) / Niter
	c_em(k) = c_em(k) / Niter
END DO

! writing the average nearest neighbors degree & clustering coefficient in a file of 1 and 100 RW models
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/knn_rw1.txt")
OPEN(1,file=trim(adjustl(output_filespath))//"/properties/knn_rw100.txt")
OPEN(2,file=trim(adjustl(output_filespath))//"/properties/c_rw1.txt")
OPEN(3,file=trim(adjustl(output_filespath))//"/properties/c_rw100.txt")

	write(0,fmt='(A8,16X,A3)') '# degree','knn'
	write(1,fmt='(A8,18X,A1)') '# degree','knn'
	write(2,fmt='(A8,18X,A1)') '# degree','c'
	write(3,fmt='(A8,18X,A1)') '# degree','c'
	DO k = 1, N
		write(0,fmt='(i8,5X,f14.8)') k, knn(k)
		write(1,fmt='(i8,5X,f14.8)') k, knn_em(k)
		if (k > 1) then
			write(2,fmt='(i8,5X,f14.8)') k, c(k)
			write(3,fmt='(i8,5X,f14.8)') k, c_em(k)
		end if
	ENDDO
CLOSE(0)
CLOSE(1)
CLOSE(2)
CLOSE(3)


!                        CONFIGURATION MODEL (CM)                     
!-------------------------------------------------------------------


! run 'Niter' RW networks and calculate its properties (knn,c)
Niter = 100
DO i = 1, Niter

	deallocate(V,ptr,knn,c)

	E_list_em = random_rewiring_process(E_list)
	CALL node_neighbors(E_list_em,D,V,ptr)
	knn = average_nearest_neighbors_degree(D,V,ptr)
 	c = clustering_coefficient(D,V,ptr)
	
	DO k = 1, N
		knn_em(k) = knn_em(k) + knn(k)  
		c_em(k) = c_em(k) + c(k)  
	END DO
	
	deallocate(E_list_em)

END DO

! average the properties over the 'Niter' RW models
DO k = 1, N
	knn_em(k) = knn_em(k) / Niter
	c_em(k) = c_em(k) / Niter
END DO

! writing the average nearest neighbors degree & clustering coefficient in a file of 1 and 100 RW models
OPEN(0,file=trim(adjustl(output_filespath))//"/properties/knn_cm1.txt")
OPEN(1,file=trim(adjustl(output_filespath))//"/properties/knn_cm100.txt")
OPEN(2,file=trim(adjustl(output_filespath))//"/properties/c_cm1.txt")
OPEN(3,file=trim(adjustl(output_filespath))//"/properties/c_cm100.txt")

	write(0,fmt='(A8,16X,A3)') '# degree','knn'
	write(1,fmt='(A8,18X,A1)') '# degree','knn'
	write(2,fmt='(A8,18X,A1)') '# degree','c'
	write(3,fmt='(A8,18X,A1)') '# degree','c'
	DO k = 1, N
		write(0,fmt='(i8,5X,f14.8)') k, knn(k)
		write(1,fmt='(i8,5X,f14.8)') k, knn_em(k)
		if (k > 1) then
			write(2,fmt='(i8,5X,f14.8)') k, c(k)
			write(3,fmt='(i8,5X,f14.8)') k, c_em(k)
		end if
	ENDDO
CLOSE(0)
CLOSE(1)
CLOSE(2)
CLOSE(3)


!-------------------------------------------------------------------
!                           SYS DYNAMICS                     
!-------------------------------------------------------------------

!                      NETWORK CHARACTERIZATION
!-------------------------------------------------------------------

E = num_edges(edgelist_filepath)
N = num_nodes(edgelist_filepath)
E_list = edges_list(edgelist_filepath,E)
D = node_degree(E_list)
CALL node_neighbors(E_list,D,V,ptr)


!                      NETWORK CHARACTERIZATION
!-------------------------------------------------------------------

E = num_edges(edgelist_filepath)
N = num_nodes(edgelist_filepath)
E_list = edges_list(edgelist_filepath,E)
D = node_degree(E_list)
deallocate(V,ptr)
CALL node_neighbors(E_list,D,V,ptr)


!                            PARAMETERS
!-------------------------------------------------------------------

T = N ! number of events simulated
M = 100 ! number of realizations
delta = 1 ! recovery rate


!        SINGLE REALIZATIONS BELOW/ABOVE THE CRITICAL POINT        
!-------------------------------------------------------------------
Ni0 = int(0.25*N) 												

DO i = 1, 6
	lambda = 0.0001*10**i 
	WRITE(a, '(f14.8)') nint(lambda*1000.0) / 1000.0
	fname_SIS = trim(adjustl(output_filespath))//"/SIS/SIS-rho_t/SIS-rho_t_lambda"// &
	trim(adjustl(a))// ".txt"
	rho = SIS_model(E,N,D,V,ptr,lambda,delta,Ni0,T)
	OPEN(0,file = trim(fname_SIS))
	DO j = 1, T
		write(0, '(f12.6)') (rho(j))
	END DO
END DO


!       		  	    LIFE-TIME DISTRIBUTION       		  	           			 
!-------------------------------------------------------------------
Ni0 = int(0.25*N) 

DO i = 1, 100 - 1
	lambda = EXP(log(0.001) + i * (log(1000.) - log(0.001)) / (100 - 1)) ! infection rate
	WRITE(a, '(f14.8)') lambda
	fname_SIS = trim(adjustl(output_filespath))//"/SIS/SIS-life_time/SIS-life_time_lambda"// trim(adjustl(a))// ".txt"
	OPEN(0,file = trim(fname_SIS))
	DO j = 1, M 
		rho = SIS_model(E,N,D,V,ptr,lambda,delta,Ni0,T)
		DO k = 1, T
			if (rho(k) == 0) then
				write(0, '(i0)') k
				exit
			end if
		END DO
	END DO
END DO


!            REALIZATIONS OVER DIFFERENT INFECTION RATES             
!          AND DIFFERENT INITIAL NUMBER OF INFECTED NODES           
!-------------------------------------------------------------------

DO i = 0, 2
	Ni0 = int(0.25*i*N) ! initial number of infected nodes
	if (i == 0) then
		Ni0 = 0.1*N
	end if
	DO j = 0, 100 - 1
		lambda = EXP(log(0.001) + j * (log(1000.) - log(0.001)) / (100 - 1)) ! infection rate
		write(a, '(f14.8)') lambda
		write(b, '(f14.8)') 0.25
		if (i == 0) then
			write(b, '(f14.8)') 0.1
		endif
		fname_SIS = trim(adjustl(output_filespath))//"/SIS/SIS-rho_st/SIS-rho_st_"// &
		trim(adjustl(b))//"N_lambda"// trim(adjustl(a))// ".txt"
		OPEN(0,file = trim(fname_SIS))
		DO k = 1, M
			rho = SIS_model(E,N,D,V,ptr,lambda,delta,Ni0,T) ! empirical prevalence as a function of time
			write(0, '(f12.6)') (rho(T))
		END DO 
		CLOSE(0)
	END DO
END DO



!-------------------------------------------------------------------
!                      FUNCTIONS & SUBROUTINES                      
!-------------------------------------------------------------------

CONTAINS ! define all the functions & subroutines


!-------------------------------------------------------------------
! SUBROUTINE: read_paths
! PURPOSE: Returns the paths of the edge list and the output files,
!		   once both are given in the command line
! INPUT: 
!   - command line for the execution of the programm
! OUTPUTS:
!   - edgelist_filepath: file name containing the edge list
!   - output_filespath: path for the output files
!-------------------------------------------------------------------

SUBROUTINE read_paths(edgelist_filepath,output_filespath)

	character(len=256) 							:: edgelist_filepath	! characters (256)
	character(len=256) 							:: output_filespath		! characters (256)
	integer*8 									:: ios					! auxiliar int.	
	logical 									:: file_exists			! boolean variables		

	! reading command argument lines
	CALL get_command_argument(1, edgelist_filepath)
	CALL get_command_argument(2, output_filespath)

	! checking that both paths are provided
	if (len_trim(edgelist_filepath) == 0 .or. len_trim(output_filespath) == 0) then
		print *, 'ERROR: Provide the paths of the edge list file and the output files.'
		stop
	end if

	! getting the file paths
	edgelist_filepath = trim(adjustl(edgelist_filepath))
	output_filespath = trim(adjustl(output_filespath))

	! checking if the edge list file exists
	OPEN(0, file=edgelist_filepath, status='old', action='read', iostat=ios)
		if (ios /= 0) then
			print *, 'ERROR opening the file.'
			stop
		end if
	CLOSE(0)

	! checking if the output files path exists
	INQUIRE(file=output_filespath, exist=file_exists)
	if (file_exists) then
		continue
	else
		print*,'ERROR: Output files path does not exist.'
	end if

END SUBROUTINE read_paths


!-------------------------------------------------------------------
! SUBROUTINE: create_directory
! PURPOSE: Creates a directory folder with the corresponding name in 
!		   the specified directory path
! INPUT: 
!   - dir_name: name of the directory
! OUTPUTS:
!   - creation of the directory
!-------------------------------------------------------------------

SUBROUTINE create_directory(dir_name)

    character(len=*), intent(in) :: dir_name
	character(len=256) :: command
	logical :: dir_exists
    integer			   :: istat

	inquire(file=dir_name, exist=dir_exists)
	if (.not. dir_exists) then	
		command = 'mkdir -p ' // trim(dir_name)
		call system(command)
	end if

END SUBROUTINE create_directory


END PROGRAM MAIN
