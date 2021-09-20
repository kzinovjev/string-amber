module splines_utilities_mod

public
private :: solve
	
contains

	
	real*8 function spline_value( x, coef )
	
		real*8, intent(in) :: x
		real*8, dimension(:,:),intent(in) :: coef
		
		integer :: idx
		real*8 :: x_max, y_max, k
		
		idx = size( coef, dim = 1 )
	  
	x_max = spline_max(coef)
		if( x < coef(1,1) ) then
			spline_value = coef(1,4)*(x-coef(1,1)) + &
						   coef(1,5)
			return
		elseif( x > x_max ) then
			k = spline_der( x_max, coef )
			y_max = coef(idx,2)*(x_max-coef(idx,1))**3 + &
					coef(idx,3)*(x_max-coef(idx,1))**2 + &
					coef(idx,4)*(x_max-coef(idx,1)) + &
					coef(idx,5)
			spline_value = k*(x-x_max) + y_max
			return
		end if
		
		do while( x < coef(idx,1) )
			idx = idx - 1
		end do
		
		spline_value = coef(idx,2)*(x-coef(idx,1))**3 + &
					   coef(idx,3)*(x-coef(idx,1))**2 + &
					   coef(idx,4)*(x-coef(idx,1)) + &
					   coef(idx,5)
	
	end function spline_value



	real*8 function spline_min(coef)
	
		real*8, dimension(:,:),intent(in) :: coef
		
		spline_min = coef(1,1)
		
	end function spline_min



	real*8 function spline_max(coef)

	real*8, dimension(:,:),intent(in) :: coef

	integer :: n

	n = size( coef, dim = 1 )

		if( coef(n,2) /= 0._8 ) then
		spline_max = coef(n,1)-coef(n,3)/(coef(n,2)*3) !Since natural Cubic Spines are used
		else
			spline_max = coef(n,1)*2-coef(n-1,1)
		end if

	end function spline_max
	
	
	
	real*8 function spline_der( x, coef ) !Derivative dy(x)/dx
	
		real*8, intent(in) :: x
		real*8, dimension(:,:),intent(in) :: coef
		
		integer :: idx
		
		idx = size( coef, dim = 1 )
		do while( x < coef(idx,1) )
			idx = idx - 1
		end do
		
		spline_der = coef(idx,2)*3*(x-coef(idx,1))**2 + &
					 coef(idx,3)*2*(x-coef(idx,1)) + &
					 coef(idx,4)
	
	end function spline_der
  
	
	
	real*8 function spline_integrate( x_min, x_max, coef )
	
		real*8, intent(in) :: x_min, x_max
		real*8, dimension(:,:),intent(in) :: coef
		
		real*8 :: tmp
		
		integer :: idx_min, idx_max, i
		
		idx_min = size( coef, dim = 1 )
		do while( x_min < coef(idx_min,1) )
			idx_min = idx_min - 1
		end do
		
		idx_max = size( coef, dim = 1 )
		do while( x_max < coef(idx_max,1) )
			idx_max = idx_max - 1
		end do
		
		tmp = 0._8
		
		if( idx_min == idx_max ) then
		
			tmp = tmp + coef(idx_min,2)*((x_max-coef(idx_max,1))**4-(x_min-coef(idx_min,1))**4)/4._8 + &
						coef(idx_min,3)*((x_max-coef(idx_max,1))**3-(x_min-coef(idx_min,1))**3)/3._8 + &
						coef(idx_min,4)*((x_max-coef(idx_max,1))**2-(x_min-coef(idx_min,1))**2)/2._8 + &
						coef(idx_min,5)*((x_max-coef(idx_max,1))-(x_min-coef(idx_min,1)))
					  
		else
		
			tmp = tmp + coef(idx_max,2)*((x_max-coef(idx_max,1))**4)/4._8 + &
						coef(idx_max,3)*((x_max-coef(idx_max,1))**3)/3._8 + &
						coef(idx_max,4)*((x_max-coef(idx_max,1))**2)/2._8 + &
						coef(idx_max,5)*((x_max-coef(idx_max,1)))
			tmp = tmp + coef(idx_min,2)*((coef(idx_min+1,1)-coef(idx_min,1))**4-(x_min-coef(idx_min,1))**4)/4._8 + &
						coef(idx_min,3)*((coef(idx_min+1,1)-coef(idx_min,1))**3-(x_min-coef(idx_min,1))**3)/3._8 + &
						coef(idx_min,4)*((coef(idx_min+1,1)-coef(idx_min,1))**2-(x_min-coef(idx_min,1))**2)/2._8 + &
						coef(idx_min,5)*((coef(idx_min+1,1)-coef(idx_min,1))-(x_min-coef(idx_min,1)))
			
			do i = idx_min+1,idx_max-1
				tmp = tmp + coef(i,2)*(coef(i+1,1)-coef(i,1))**4/4._8 + &
							coef(i,3)*(coef(i+1,1)-coef(i,1))**3/3._8 + &
							coef(i,4)*(coef(i+1,1)-coef(i,1))**2/2._8 + &
							coef(i,5)*(coef(i+1,1)-coef(i,1))			 
			end do
		
		end if
		
		spline_integrate = tmp
	
	end function spline_integrate



	subroutine CubicSplines( x, y, coef )

		real*8, dimension(:), intent(in) :: x, y
		real*8, dimension(size(x)-1,5), intent(out) :: coef

		integer :: i, n
	
		real*8, dimension(size(x)-1) :: h
		real*8, dimension(size(x)-2) :: d, b
		real*8, dimension(size(x)) :: s

		n = size(x)
		
		do i = 1, n-1
			h(i) = x(i+1)-x(i)
		end do
		
		do i = 1, n-2
			d(i) = (h(i)+h(i+1))*2
			b(i) = (y(i+2)-y(i+1))/h(i+1) - (y(i+1)-y(i))/h(i)
		end do
		b = b * 6
		
		s(1) = 0._8
		s(n) = 0._8
		s(2:n-1) = Thomas( n-2, d, h(1:n-2), h(2:n-1), b )
		
		coef(:,1) = x(1:n-1)
		coef(:,2) = (s(2:n)-s(1:n-1))/(h(1:n-1)*6)
		coef(:,3) = s(1:n-1)/2._8
		coef(:,4) = (y(2:n)-y(1:n-1))/h(1:n-1) - ((s(2:n)+2*s(1:n-1))/6._8)*h(1:n-1)
		coef(:,5) = y(1:n-1)

	end subroutine CubicSplines



	subroutine CubicSplinesND( x, y, coef )

		real*8, dimension(:), intent(in) :: x
		real*8, dimension(:,:), intent(in) :: y
		real*8, dimension(size(x)-1,5,size(y, dim=1)), intent(out) :: coef

		integer :: i, d

		d = size(y, dim=1)
		do i = 1, d
			call CubicSplines(x, y(i,:), coef(:,:,i))
		end do

	end subroutine CubicSplinesND



	function spline_value_ND(x, coef)

		real*8, intent(in) :: x
		real*8, dimension(:,:,:), intent(in) :: coef
		real*8, dimension(size(coef, dim=3)) :: spline_value_ND

		integer :: i, d
		
		d = size(coef, dim=3)
		do i = 1, d
			spline_value_ND(i) = spline_value(x, coef(:,:,i))
		end do

	end function spline_value_ND
	
	
	
	real*8 function spline_grad_mod_ND(x, coef)

		real*8, intent(in) :: x
		real*8, dimension(:,:,:), intent(in) :: coef

		spline_grad_mod_ND = sqrt(sum(spline_der_ND(x, coef)**2))
	
	end function spline_grad_mod_ND
	

	
	function spline_der_ND(x, coef)

		real*8, intent(in) :: x
		real*8, dimension(:,:,:), intent(in) :: coef
		real*8, dimension(size(coef, dim=3)) :: spline_der_ND

		integer :: i, d
		
		d = size(coef, dim=3)
		do i = 1, d
			spline_der_ND(i) = spline_der(x, coef(:,:,i))
		end do

	end function spline_der_ND
	
	
	
	subroutine CubicSplinesFitND( x, y, coef )

		real*8, dimension(:), intent(in) :: x
		real*8, dimension(:,:), intent(in) :: y
		real*8, dimension(:,:,:) :: coef
		
		integer :: d
		d = size(y, dim=1)
		
		do i = 1, d
			call CubicSplinesFit(x, y(i,:), coef(:,:,i))
		end do
	
	end subroutine CubicSplinesFitND
	
	
	
	subroutine CubicSplinesFit( a, b, coef )
	
		real*8, dimension(:), intent(in) :: a, b	  !Vectors with coordinates of data points
		real*8, dimension(:,:), intent(out) :: coef !coef(:,1) - x_i values, coef(:,2:5) - spline coefficients
		
		integer :: n, np	!number of data points
		integer :: i, j, idx
		real*8 :: h		! x_i+1 - x_i
		real*8, dimension(size(coef, dim=1)+1) :: x, y	!Coordinates of spline reference points
		integer, dimension(size(a)) :: aidx	!indexes of splines to which each data point belongs
		real*8, dimension(size(coef, dim=1)+1) :: Mval !Final M values
		real*8, dimension(size(coef, dim=1)+1) :: cpb
		real*8, dimension(size(coef, dim=1)+1,size(coef, dim=1)+1) :: M  !m = My
		real*8, dimension(size(a),size(coef, dim=1)+1) :: C !C_ji = d/dy_i S(a_j)
		real*8, dimension(size(coef, dim=1)+1,size(a)) :: CT
		real*8, dimension(size(coef, dim=1)+1,size(coef, dim=1)+1) :: Cp !Coefficient matrix for the final equation system

		n = size(coef, dim=1)+1
		np = size( a )
		
		x(1) = minval( a )
		x(n) = maxval( a )
		h = (x(n) - x(1))/real( n-1, kind = 8 )

		do i = 2, n-1
			x(i) = x(1) + h*(i-1)
		end do
		
		aidx = min( int( (a-x(1))/h ) + 1, n-1 ) !Defining, to which spline each data point belongs
		
		M = SplineMatrix( n, h )

		!---Expressing s(a_i) as a linear combination of y---
		!----and putting the coefficients into the matrix----
		C = 0.
		do i = 1, np
			idx = aidx(i)
			do j = 1, n
				C(i,j) = (a(i)-x(idx))**3/(6*h)*(M(idx+1,j)-M(idx,j))+&
				(a(i)-x(idx))**2/2._8*M(idx,j)-&
				(a(i)-x(idx))*h/6._8*(M(idx+1,j)+M(idx,j)*2)
			end do
			C(i,idx) = C(i,idx) + 1-(a(i)-x(idx))/h
			C(i,idx+1) = C(i,idx+1) + (a(i)-x(idx))/h
		end do
		!----------------------------------------------------

		!---Finding y by solving the least squares equation system---
		CT = transpose( C )
		Cp = matmul( CT, C )
		do i = 1, n
			cpb(i) = dot_product( CT(i,:), b )
		end do

		call solve(n, Cp, cpb, y)

		do i = 1, n
			Mval(i) = dot_product( M(i,:), y )
		end do

		coef(:,1) = x(1:n-1)
		coef(:,2) = (Mval(2:n)-Mval(1:n-1))/(6*h)
		coef(:,3) = Mval(1:n-1)/2.
		coef(:,4) = (y(2:n)-y(1:n-1))/h - ((Mval(2:n)+2*Mval(1:n-1))/6.)*h
		coef(:,5) = y(1:n-1)

		!------------------------------------------------------------

	end subroutine CubicSplinesFit
	
	
	subroutine solve(n, A, b, x)
		real*8, dimension(n,n), intent(in) :: A
		real*8, dimension(n), intent(in) :: b
		real*8, dimension(n), intent(out) :: x
		
		integer :: info
		integer, dimension(n) :: ipiv
		real*8, dimension(n,n) :: Mtmp

		Mtmp = A
		call dgetrf(n, n, Mtmp, n, ipiv, info)
		if( info /= 0 ) then
			write(*,*) "dgetrf in solve failed: ", info
			stop
		end if
	
		x = b
		call dgetrs('N', n, 1, Mtmp, n, ipiv, x, n, info)
		if( info /= 0 ) then
			write(*,*) "dgetrs in solve failed: ", info
			stop
		end if

	end subroutine solve

	

	function SplineMatrix( n, h ) 
		
		! Generates the matrix M with coefficients of m expansion as a linear combination
		! of y: m = My, m - (second derivative of spline)/2, all other coefficients can
		! be expressed as a function of m and y.
		
		integer, intent(in) :: n
		real*8, intent(in) :: h
		real*8, dimension(n,n) :: SplineMatrix
		real*8, dimension(n-2,n-2) :: X
		integer, dimension(n-2,n) :: Y
		
		integer :: i, j
		real*8, parameter  :: lambda = 1.31695789692482_8
		
		do i = 1, n-2
			do j = 1, n-2
				X(i,j) = (1-2*iand( i+j, 1 ))*( cosh( (n-1-abs(j-i))*lambda ) - cosh( (n-1-i-j)*lambda ) )/(2*sinh(lambda)*sinh((n-1)*lambda))
			end do
		end do
		
		Y = 0.
		do i = 1, n-2
			Y(i,i) = 1
			Y(i,i+1) = -2
			Y(i,i+2) = 1
		end do
		
		SplineMatrix = 0.
		SplineMatrix(2:n-1,:) = matmul( X, Y )*(6/h**2)
		
		return
		
	end function SplineMatrix

  
  
	function Thomas( n, d, subd, supd, b )
	!Solves a tridiagonal system of linear equation
	!n - number of equations
	!d - main diagonal elements
	!subd - sub-diagonal elements
	!supd - sup-diagonal elements
	!b - right side values
	!WARNING! d and b values are not conserved!
	
		integer, intent(in) :: n
		real*8, dimension(:), intent(inout) :: d, subd, supd, b
		real*8, dimension(n) :: Thomas
		real*8, dimension(:), allocatable :: c
		
		integer :: i
		real*8  :: tmp
		
		do i = 2, n
			tmp = subd(i)/d(i-1)
			d(i) = d(i) - tmp*supd(i-1)
			b(i) = b(i) - tmp*b(i-1)
		end do
		
		Thomas(n) = b(n)/d(n)

		do i = n-1, 1, -1
			Thomas(i) = (b(i)-supd(i)*Thomas(i+1))/d(i)
		end do
	  
	end function Thomas  


end module splines_utilities_mod
