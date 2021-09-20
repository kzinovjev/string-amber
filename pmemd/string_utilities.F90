module string_utilities_mod

	public
	
contains
	
	!================================
	subroutine read_string(file, string)
	
		character(len=*), intent(in) :: file
		real*8, dimension(:,:), intent(out) :: string
		
		integer :: u
		
		u = next_unit()
		open(unit=u, file=trim(file), status="old")
		read(u,*) string
		close(u)
	
	end subroutine read_string
	!================================
	
	
	
	!===========================
	integer function next_unit()

		integer :: i
		logical :: opened

		next_unit = -1

		do i = 20, 100
			inquire(i, opened=opened)
			if (.not. opened) then
				next_unit = i
				exit
			end if
		end do
	
	end function next_unit
	!===========================
	
	
	
	!================================
	subroutine matinv( M, Minv )
		real*8, dimension(:), intent(in) :: M
		real*8, dimension(:), intent(out) :: Minv
		
		real*8, dimension(:,:), allocatable :: Mtmp

		integer :: i, j, idx, n, inf

		if (size(M) /= size( Minv) ) call write_error("matinv: size mismatch")
		
		n = int(sqrt(real(size(M)*2))) !size(M)=n*(n+1)/2 => n**2 < size(M)*2 < (n+1)**2
		allocate(Mtmp(n,n))
		
		idx = 0
		do i = 1, n
			do j = 1, i
				idx = idx + 1
				Mtmp(i,j) = M(idx)
			end do
		end do

		call dpotrf( 'L', n, Mtmp, n, inf )
		if( inf /= 0 ) then
			write(*,*) "dpotrf in matinv failed: ", inf
			stop
		end if

		call dpotri( 'L', n, Mtmp, n, inf )
		if( inf /= 0 ) then
			write(*,*) "dpotri in matinv failed: ", inf
			stop
		end if
		
		Minv = 0._8
		idx = 0
		do i = 1, n
		    do j = 1, i
				idx = idx+1
				Minv(idx) = Mtmp(i,j)
		    end do
		end do

		
		
		deallocate(Mtmp)

	end subroutine matinv
	!================================
	
	
	!================================
	subroutine write_error(text)
	
		character(len=*) :: text
		
		write(*,*) text
		stop
	
	end subroutine write_error
	!================================



	!================================
	integer function lines( filename )

		character(len=*), intent(in) :: filename
		integer :: stat, u

		u = next_unit()
		open(unit=u, file=filename, status="old")
		lines = 0
		do
			read(u,*,end=10)
			lines = lines + 1
		end do
10		close(u)

	end function lines
	!================================



end module string_utilities_mod
