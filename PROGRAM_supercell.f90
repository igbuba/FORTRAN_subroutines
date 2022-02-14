!
!  This module contains subroutines and functions for compute_force.90 code
!	10th Feb., 2021
!
!  *Updated 14th Feb, 2022
! 
!  This program compute the supercell for MD or PIOUD-LD using 
!  primitive unit cell and q-space phonon grid from mat2R
!
!-------------------------------------------------------------------------
  !
  !
PROGRAM compute_supercell			! 
  USE kinds,    ONLY 	   : DP		! use dble precision		
  USE input_fc, ONLY 	   : read_fc2, forceconst2_grid, ph_system_info, & ! from program module "input_fc" 
                       	     read_system, aux_system		           ! read the fllwn info
  USE random_numbers, ONLY : randy 
  USE asr2_module,    ONLY : impose_asr2
  USE constants,      ONLY : BOHR_RADIUS_ANGS

  IMPLICIT NONE
  REAL(kind=DP),ALLOCATABLE :: aa(:,:), bb(:,:), r(:,:)
  REAL(DP),PARAMETER        :: ANGS_TO_BOHR = 1/BOHR_RADIUS_ANGS
  INTEGER       	    :: i, j, k, l, uni, ios, unit_sc
  !
  TYPE(ph_system_info)      :: S, SC   ! Define this TYPE() as S 
  TYPE(forceconst2_grid)    :: fc2	
  !
  CALL read_fc2("mat2R", S, fc2) ! Read file "mat2R" and put info into S and fc2
  CALL aux_system(S)
  CALL impose_asr2( "simple", S%nat, fc2) ! Apply asr to FCs
  !
  ALLOCATE(aa(3,3))
  aa(:,1) = S%at(:,1)*fc2%nq(1)
  aa(:,2) = S%at(:,2)*fc2%nq(2)
  aa(:,3) = S%at(:,3)*fc2%nq(3)
  ALLOCATE(bb(3,3))
  bb(:,1) = S%bg(:,1)/fc2%nq(1) ! vectors in reciprocal lattice of the sc 
  bb(:,2) = S%bg(:,2)/fc2%nq(2)
  bb(:,3) = S%bg(:,3)/fc2%nq(3)
  !
  ALLOCATE(r(3,S%nat*fc2%n_R))
  OPEN(newunit=unit_sc,file="supercell.dat",status="unknown") 
  WRITE(unit_sc,*) "&SYSTEM "
  WRITE(unit_sc,*) "ibrav =", 0
  WRITE(unit_sc,*) "ntyp  =", LEN(S%atm)
  WRITE(unit_sc,*) "nat   =", fc2%n_R*S%nat
  WRITE(unit_sc,*) "/ "
  WRITE(unit_sc,*) "CELL_PARAMETERS bohr"
  WRITE(unit_sc,*) S%at(:,1)*fc2%nq(1)*S%alat
  WRITE(unit_sc,*) S%at(:,2)*fc2%nq(2)*S%alat
  WRITE(unit_sc,*) S%at(:,3)*fc2%nq(3)*S%alat
  !
  WRITE(unit_sc,*) "ATOMIC_POSITIONS crystal"
  PRINT*, "------------------------------------------------------------"
  PRINT*, "======= Supercell from mat2R grid ========"
  PRINT*, 
   l = 0
   DO j = 1, fc2%n_R
   DO i = 1, S%nat
   !DO k = 1, 3
	l = l + 1
        r(:,l) =  S%tau(:,i) + fc2%xR(:,j) 	   ! CRYSTAL    
        !r(:,l) =  S%tau(:,i) + fc2%xR(:,j)*S%alat 	! BOHR   
        !r(:,l) =  S%tau(:,i) + fc2%xR(:,j)*S%alat*BOHR_RADIUS_ANGS ! ANGSTROM
        CALL cryst_to_cart(1, r(:,l), bb, -1) ! only when crystal unit above
  WRITE(unit_sc,'(a3,3f18.10)') S%atm(S%ityp(i)), r(:,l) ! S%tau(:,i)tau_sc(i_at:k) + n_R(:,i) 
   ENDDO
   ENDDO
   CLOSE(unit_sc)
  !
  DEALLOCATE(r) 
  DEALLOCATE(bb) 
  !
END PROGRAM compute_supercell
