module multiCV_module_mod

    use CV_module_mod
    use string_utilities_mod, only : next_unit

    public :: MBOND_TYPE, MANGLE_TYPE, MDIHEDRAL_TYPE, MPPLANE_TYPE, &
              prepare_multiCV, set_multiCV_used_atoms_mask, calculate_multiCV
    private

    integer, parameter :: MBOND_TYPE = 5, &
                          MANGLE_TYPE = 6, &
                          MDIHEDRAL_TYPE = 7, &
                          MPPLANE_TYPE = 8

    integer, dimension(MBOND_TYPE:MPPLANE_TYPE) :: TYPE_NCENTERS = (/ 2, 3, 4, 4 /)
    real*8, dimension(12) :: centers, centers_gradient


    type atom_list
        integer :: size
        integer, dimension(:), allocatable :: atoms
    end type atom_list

    type multiCV
		integer :: CV_type, ncenters
		type(atom_list), dimension(4) :: centers
    end type multiCV


    integer :: nmultiCV
    type(multiCV), dimension(:), target, allocatable :: multiCV_list

contains


    subroutine extend_multiCV_list

        type(multiCV), dimension(:), allocatable :: new_list
        integer :: i, current_size

        if (allocated(multiCV_list)) then
            current_size = size(multiCV_list)
            allocate(new_list(current_size + 1))
            do i = 1, current_size
                new_list(i) = multiCV_list(i)
            end do
            deallocate(multiCV_list)
            call move_alloc(new_list, multiCV_list)
        else
            allocate(multiCV_list(1))
        end if

    end subroutine extend_multiCV_list


    subroutine prepare_multiCV(CV_type, filename, index)

        integer,intent(in) :: CV_type
        character(len=*), intent(in) :: filename
        type(multiCV), pointer :: CV

        integer, intent(out) :: index

        integer :: i, natoms, u

        call extend_multiCV_list
        index = size(multiCV_list)
        CV => multiCV_list(index)

        u = next_unit()
        open(unit=u, file=filename, status="old")

        CV%CV_type = CV_type
        CV%ncenters = TYPE_NCENTERS(CV_type)
        do i = 1, CV%ncenters
            read(u,*) natoms
            CV%centers(i)%size = natoms
            allocate(CV%centers(i)%atoms(natoms))
            read(u,*) CV%centers(i)%atoms
        end do

        close(u)

    end subroutine prepare_multiCV


    subroutine set_multiCV_used_atoms_mask(index, used_atoms_mask)

        integer, intent(in) :: index
        logical, dimension(:), intent(inout) :: used_atoms_mask

        integer :: i, j
        type(multiCV), pointer :: CV

        CV => multiCV_list(index)

        do i = 1, CV%ncenters
            do j = 1, CV%centers(i)%size
                used_atoms_mask(CV%centers(i)%atoms(j)) = .true.
            end do
        end do

    end subroutine set_multiCV_used_atoms_mask


    subroutine calculate_multiCV(index, x, value, gradient)

        integer, intent(in) :: index
        real*8, dimension(:), intent(in) :: x
        real*8, intent(out) :: value
        real*8, dimension(:), intent(inout) :: gradient

        integer :: i
        type(multiCV), pointer :: CV

        CV => multiCV_list(index)

        do i = 1, CV%ncenters
            call set_xyz(i, centers, get_center(CV%centers(i)%atoms, x))
        end do

        centers_gradient = 0.
        select case (CV%CV_type)
            case (MBOND_TYPE)
                call CV_bond(centers, (/ 1, 2 /), value, centers_gradient)
            case (MANGLE_TYPE)
                call CV_angle(centers, (/ 1, 2, 3 /), value, centers_gradient)
            case (MDIHEDRAL_TYPE)
                call CV_dihedral(centers, (/ 1, 2, 3, 4 /), value, centers_gradient)
            case (MPPLANE_TYPE)
                call CV_pointplane(centers, (/ 1, 2, 3, 4 /), value, centers_gradient)
            end select

        do i = 1, CV%ncenters
            gradient = gradient + unpack_gradient(size(gradient) / 3, &
                                                  CV%centers(i)%atoms, &
                                                  get_xyz(i, centers_gradient))
        end do

    end subroutine calculate_multiCV


    function get_center(indexes, x) result(center)

        integer, dimension(:), intent(in) :: indexes
        real*8, dimension(:), intent(in) :: x
        real*8, dimension(3) :: center, reference, xyz

        integer :: i

        center = 0.

        !All xyz will be projected to the nearest image to the reference,
        !to ensure that the center is calculated correctly in periodic systems
        reference = get_xyz(indexes(1), x)

        do i = 1, size(indexes)
            xyz = get_xyz(indexes(i), x) - reference
            call map_periodic(xyz)
            center = center + reference + xyz
        end do

        center = center / real(size(indexes), kind=8)

    end function get_center


    function unpack_gradient(natoms, indexes, center_gradient) result(gradient)

        integer, intent(in) :: natoms
        integer, dimension(:), intent(in) :: indexes
        real*8, dimension(:), intent(in) :: center_gradient

        integer :: i
        real*8, dimension(natoms*3) :: gradient

        gradient = 0.
        do i = 1, size(indexes)
            call set_xyz(indexes(i), gradient, center_gradient)
        end do
        gradient = gradient / real(size(indexes), kind=8)

    end function unpack_gradient


    function get_xyz(index, x) result(xyz)

        integer, intent(in) :: index
        real*8, dimension(:), intent(in) :: x
        real*8, dimension(3) :: xyz

        xyz = x(index*3-2:index*3)

    end function get_xyz


    subroutine set_xyz(index, x, xyz)
        integer, intent(in) :: index
        real*8, dimension(:), intent(inout) :: x
        real*8, dimension(3), intent(in) :: xyz

        x(index*3-2:index*3) = xyz

    end subroutine set_xyz

end module multiCV_module_mod
