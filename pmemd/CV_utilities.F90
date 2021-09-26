! <compile=optimized>
!===============================================================================
!				 The CV utilities Module
!===============================================================================
! 
!	Contains utilities for CV calculation.
! 
!   Subroutines:
!
!   prepare_CVs
!	update_CVs
!
!===============================================================================
module CV_utilities_mod

        use CV_module_mod, only : CV_bond, CV_angle, CV_dihedral, CV_pointplane
        use string_utilities_mod, only : next_unit, write_error
        use multiCV_module_mod

		use prmtop_dat_mod, only : natom

	public
	private :: CV
	
	type CV
		integer :: CV_type
		integer, dimension(4) :: atoms
		real*8, dimension(2) :: box
		logical :: periodic
		real*8  :: period
		integer :: mindex
	end type CV
	
	integer :: nCV   !number of collective variables
	integer :: msize !size of the metric tensor array (nCV*(nCV+1)/2)

    integer :: n_used_atoms !number of atoms involved in any CV
    integer, dimension(:), allocatable :: used_atoms !atoms involved in any CV
	
	real*8, dimension(:), allocatable :: CVs, M
	real*8, dimension(:,:), allocatable :: Jacobian
	type(CV), dimension(:), allocatable :: CV_list

contains
	

	!================================
	subroutine prepare_CVs( filename )
	
		character(len=*), intent(in) :: filename
		
		integer :: i, u
		integer, dimension(4) :: atoms
		real*8, dimension(2) :: box
		character*20 :: COLVAR_type
		character*200 :: atoms_file

		namelist /COLVAR/ COLVAR_type, atoms, box, atoms_file
	
		u = next_unit()
		open(unit=u, file=filename, status="old")
		read(u,*) nCV
		msize = nCV*(nCV+1)/2
		allocate(CVs(nCV), CV_list(nCV), Jacobian(natom*3,nCV), M(msize))
		
		!---read the definitions of CVs---
		do i = 1, nCV
			box = (/ -huge(0._8), huge(0._8) /)
			read(u,NML=COLVAR)
			CV_list(i)%periodic = .false.
			select case ( trim( adjustl( COLVAR_type ) ) )
				case ("BOND")
					CV_list(i)%CV_type = 1
				case ("ANGLE")
					CV_list(i)%CV_type = 2
				case ("DIHEDRAL")
					CV_list(i)%CV_type = 3
					CV_list(i)%periodic = .true.
					CV_list(i)%period = 360._8
				case ("PPLANE")
					CV_list(i)%CV_type = 4
				case ("MBOND")
					CV_list(i)%CV_type = MBOND_TYPE
					call prepare_multiCV(MBOND_TYPE, atoms_file, CV_list(i)%mindex)
				case ("MANGLE")
					CV_list(i)%CV_type = MANGLE_TYPE
					call prepare_multiCV(MANGLE_TYPE, atoms_file, CV_list(i)%mindex)
				case ("MDIHEDRAL")
					CV_list(i)%CV_type = MDIHEDRAL_TYPE
					CV_list(i)%periodic = .true.
					CV_list(i)%period = 360._8
					call prepare_multiCV(MDIHEDRAL_TYPE, atoms_file, CV_list(i)%mindex)
				case ("MPPLANE")
					CV_list(i)%CV_type = MPPLANE_TYPE
					call prepare_multiCV(MPPLANE_TYPE, atoms_file, CV_list(i)%mindex)

			end select
			CV_list(i)%atoms = atoms
			CV_list(i)%box = box
		end do
		!---------------------------------
		close(u)
        call populate_used_atoms
	
	end subroutine prepare_CVs
	!================================



    !================================
    subroutine populate_used_atoms

        integer, dimension(4), parameter :: n_CV_atoms = (/ 2, 3, 4, 4 /)
        logical, dimension(natom) :: used_atoms_mask

        integer :: i, j

        used_atoms_mask = .false.
        do i = 1, nCV
            do j = 1, n_CV_atoms(CV_list(i)%CV_type)
                used_atoms_mask(CV_list(i)%atoms(j)) = .true.
            end do
        end do

        n_used_atoms = count(used_atoms_mask)
        allocate(used_atoms(n_used_atoms))

        j = 1
        do i = 1, natom
            if (.not. used_atoms_mask(i)) cycle
            used_atoms(j) = i
            j = j + 1
        end do

    end subroutine populate_used_atoms
    !================================

	
	
	!================================
	subroutine update_CVs(x)
	
		real*8, dimension(natom*3), intent(in) :: x
		
		integer :: i
		
		Jacobian = 0._8
		do i = 1, nCV
			select case ( CV_list(i)%CV_type )
				case (1)
					call CV_bond( x, CV_list(i)%atoms, CVs(i), Jacobian(:,i) )
				case (2)
					call CV_angle( x, CV_list(i)%atoms, CVs(i), Jacobian(:,i) )
				case (3)
					call CV_dihedral( x, CV_list(i)%atoms, CVs(i), Jacobian(:,i) )
				case (4)
					call CV_pointplane( x, CV_list(i)%atoms, CVs(i), Jacobian(:,i) )
				case (MBOND_TYPE, MANGLE_TYPE, MDIHEDRAL_TYPE, MPPLANE_TYPE)
					call calculate_multiCV(CV_list(i)%mindex, x, CVs(i), Jacobian(:,i) )
				end select
		end do
		
		call get_M(Jacobian, M)
   
	end subroutine update_CVs
	!================================

	
	
	!================================
	subroutine to_box(A)
	
		real*8, dimension(nCV), intent(inout) :: A
		
		integer :: i
		
		do i = 1, nCV
			A(i) = min(CV_list(i)%box(2), &
					   max(CV_list(i)%box(1), &
						   A(i)))
		end do
	
	end subroutine to_box
	!================================
	
	
	
	!================================
	!Projects a path in the CV space so that there are no jumps in CV values due to periodic coordinates
	subroutine to_continuous(A)
	
		real*8, dimension(:,:), intent(inout) :: A
		
		integer :: i, n
		
		if (size(A, dim=1) /= nCV) call write_error("to_continuous: size mismatch")
		n = size(A, dim=2)
		
		do i = 2, n
			A(:,i) = A(:,i-1) + map_periodic(A(:,i) - A(:,i-1))
		end do
		
	
	end subroutine to_continuous
	!================================
	
	
	
	!================================
	!Projects a path in the CV space so that all the values of CVs are in [-period/2:period/2] range
	subroutine to_period(A)
	
		real*8, dimension(:,:), intent(inout) :: A
		
		integer :: i, n
		
		if (size(A, dim=1) /= nCV) call write_error("to_continuous: size mismatch")
		n = size(A, dim=2)

		do i = 1, n
			A(:,i) = map_periodic(A(:,i))
		end do
		
	end subroutine to_period
	!================================
	
	
	
	!================================
	subroutine get_M(Jacobian, M)

		use prmtop_dat_mod, only : atm_mass

		real*8, dimension(natom*3,nCV), intent(in) :: Jacobian
		real*8, dimension(msize), intent(out) :: M
		
		integer :: i, j, k, l, idx

		idx = 0
		M = 0._8
		do i = 1, nCV
			do j = 1, i
				idx = idx + 1
				do l = 1, n_used_atoms
                    k = used_atoms(l)
                    if (atm_mass(k) == 0) cycle
					M(idx) = M(idx) + &
							 dot_product( Jacobian(k*3-2:k*3,i), &
										  Jacobian(k*3-2:k*3,j) ) / &
							 atm_mass(k)
				end do
			end do
		end do  

	end subroutine get_M
	!================================
	
	
	
	!================================
	real*8 function CV_distance(A, B, M)
	
		real*8, dimension(nCV), intent(in) :: A, B
		real*8, dimension(:), intent(in) :: M
		
		CV_distance = sqrt(CV_distance2(A, B, M))
	
	end function CV_distance
	!================================
	
	
	
	!================================
	real*8 function CV_distance2(A, B, M)
	
		real*8, dimension(nCV), intent(in) :: A, B
		real*8, dimension(:), intent(in) :: M
		
		CV_distance2 = len2_M(map_periodic(A-B), M)
	
	end function CV_distance2
	!================================
	
	
	
	!================================
	subroutine normalize_M(A, M)
	
		real*8, dimension(:), intent(inout) :: A
		real*8, dimension(:), intent(in) :: M
		
		A = A/len_M(A, M)
	
	end subroutine normalize_M
	!================================

	

	!================================
	real*8 function dist_M(A, B, M)
	
		real*8, dimension(:), intent(in) :: A, B
		real*8, dimension(:), intent(in) :: M
		
		dist_M = len_M(A-B, M)
	
	end function dist_M
	!================================
	
	

	!================================
	real*8 function len_M(A, M)
	
		real*8, dimension(:), intent(in) :: A
		real*8, dimension(:), intent(in) :: M
		
		len_M = sqrt(len2_M(A, M))
	
	end function len_M
	!================================
	
	
	
	!================================
	real*8 function len2_M(A, M)
	
		real*8, dimension(:), intent(in) :: A
		real*8, dimension(:), intent(in) :: M
		
		len2_M = dot_product_M(A, A, M)
	
	end function len2_M
	!================================	
	


	!================================
	real*8 function dot_product_M(A, B, M)
	
		real*8, dimension(:), intent(in) :: A, B
		real*8, dimension(:), intent(in) :: M
		
		dot_product_M = dot_product(A, matmulp(M, B))
	
	end function dot_product_M
	!================================
	
	
	
	!================================
	function matmulp(M, A)
	
		real*8, dimension(:), intent(in) :: M, A
		real*8, dimension(size(A)) :: matmulp
		
		integer :: i, j, n, idx
		
		n = size(A)
		if (size(M) /= n*(n+1)/2) call write_error("matmulp: array size mismatch")
		
		matmulp = 0._8
		idx = 0
		do i = 1, n
			do j = 1, i-1
				idx = idx + 1
				matmulp(i) = matmulp(i) + M(idx)*A(j)
				matmulp(j) = matmulp(j) + M(idx)*A(i)
			end do
			idx = idx + 1
			matmulp(i) = matmulp(i) + M(idx)*A(i)
		end do
		
	end function matmulp
	!================================
	
	
	
	!================================
	function map_periodic(A)
	
		real*8, dimension(nCV), intent(in) :: A
		real*8, dimension(nCV) :: map_periodic
		
		integer :: i
		
		map_periodic = A
		do i = 1, nCV
			if (CV_list(i)%periodic) &
				map_periodic(i) = A(i) - sign( CV_list(i)%period, A(i) ) * &
											 int( abs(A(i)/CV_list(i)%period) + 0.5 )
		end do
	
	end function map_periodic
	!================================
	
	
end module CV_utilities_mod
