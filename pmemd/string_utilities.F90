module string_utilities

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
	
	
	
	!================================
	subroutine coordinates_write( x, file )

		use nblist, only: a, b, c, alpha, beta, gamma
#include "../include/memory.h"			
	
		real*8, dimension(:), intent(in) :: x
		character( len = * ), intent(in) :: file
		integer :: unit
		
		unit = next_unit()
		open( unit = unit, file = file, status = "replace" )
		write(unit,"(A)") "default"
		write(unit,"(I6)") natom
		write(unit,"(6F12.7)") x(1:natom*3)
		write(unit,"(6F12.7)") a, b, c, alpha, beta, gamma
		close(unit)
	
	end subroutine coordinates_write
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
#include "../include/dprec.fh"
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

		call D_OR_S()potrf( 'L', n, Mtmp, n, inf )
		if( inf /= 0 ) then
			write(*,*) "dpotrf in matinv failed: ", inf
			stop
		end if

		call D_OR_S()potri( 'L', n, Mtmp, n, inf )
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
	
	
	
	!================================
	integer function lines2( filename )
#include "parallel.h"
		character(len=*), intent(in) :: filename
		integer :: stat, u
		character*200 :: tmpfile
		
		write(tmpfile,*) masterrank
		tmpfile = "string_stilities_lines_tmp_file"//trim(adjustl(tmpfile))
		call system("wc -l "//trim(adjustl(filename))//" > "//trim(tmpfile))
		u = next_unit()
		open(unit=u, file=trim(tmpfile), status="old")
		read(u,*) lines2
		close(u)
		call system("rm -rf "//trim(tmpfile))
		
		write(*,*) lines2
		flush(6)

	end function lines2
	!================================



end module string_utilities
