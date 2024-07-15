
!-------------------------------------------------------------------
! MODULE: Random Number Generator
! AUTHOR: Ivan Casanovas Rodríguez
! DATE: July, 2024
! DESCRIPTION: This module generates a random number by initializing
!			   the system wwith a seed. 
! WARNING: It is necessary to include the file "r1279block.h" in the
!		   same folder of this module.
!-------------------------------------------------------------------


MODULE RANDOM_NUMBER_GENERATOR

IMPLICIT NONE

CONTAINS  ! define all the functions & subroutines


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


END MODULE RANDOM_NUMBER_GENERATOR
