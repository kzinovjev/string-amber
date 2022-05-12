module CV_module_mod

    use prmtop_dat_mod, only : ifbox
    use pbc_mod, only : pbc_box

    public :: CV_bond, CV_angle, CV_dihedral, &
              CV_symdih, CV_asymdih, CV_pointplane, map_periodic
	private
contains



    !-------------------------
    subroutine CV_bond( x, atm, bond, gradient )
    
        real*8, dimension(:), intent(in) :: x
        integer, dimension(2), intent(in) :: atm
        real*8, intent(out) :: bond
        real*8, dimension(:), intent(inout) :: gradient

	call dist( x(atm(1)*3-2:atm(1)*3), x(atm(2)*3-2:atm(2)*3), bond, &
		   gradient(atm(1)*3-2:atm(1)*3), gradient(atm(2)*3-2:atm(2)*3) )

    end subroutine CV_bond
    !-------------------------
    
    
    
    !-------------------------
    subroutine CV_angle( x, atm, angle, gradient )
    
        real*8, dimension(:), intent(in) :: x
        integer, dimension(3), intent(in) :: atm
        real*8, intent(out) :: angle
        real*8, dimension(:), intent(inout) :: gradient
        
        real*8, dimension(3,3) :: xatm, dx
        real*8, dimension(3) :: ab, bc, nab, nbc
        real*8 :: dab, dbc, cosang, d_acos
        integer :: i
        
        do i = 1, 3
            xatm(:,i) = x(atm(i)*3-2:atm(i)*3)
        end do
        
        ab = xatm(:,1) - xatm(:,2)
        bc = xatm(:,3) - xatm(:,2)
        call map_periodic( ab )
        call map_periodic( bc )
        dab = sqrt( sum( ab**2 ) )
        dbc = sqrt( sum( bc**2 ) )
        nab = ab/dab
        nbc = bc/dbc
        
        cosang = dot_product( nab, nbc )
        if( cosang > 1._8  ) cosang =  1._8
        if( cosang < -1._8 ) cosang = -1._8
        
        angle = acos( cosang )*57.295780_8

        d_acos = -57.295780_8/sqrt( 1._8 - cosang**2 )
        dx(:,1) = ( nbc - cosang*nab )/dab
        dx(:,3) = ( nab - cosang*nbc )/dbc
        dx(:,2) = -( dx(:,1) + dx(:,3) )
        
        do i = 1, 3
            gradient(atm(i)*3-2:atm(i)*3) = gradient(atm(i)*3-2:atm(i)*3) + dx(:,i)*d_acos
        end do
        
    end subroutine CV_angle
    !-------------------------



    !-------------------------
    subroutine CV_dihedral( x, atm, dihedral, gradient )

        real*8, dimension(:), intent(in) :: x
        integer, dimension(4), intent(in) :: atm
        real*8, intent(out) :: dihedral
        real*8, dimension(:), intent(inout) :: gradient

        integer :: i, j
        real*8, dimension(3, 4) :: grad

        call CV_dihedral_aux(x, atm, dihedral, grad)

        do i = 1, 4
            j = atm(i)*3
            gradient(j-2:j) = gradient(j-2:j) + grad(:,i)
        end do

    end subroutine CV_dihedral
    !-------------------------


    !-------------------------
    subroutine CV_symdih( x, atm, symdih, gradient )

        real*8, dimension(:), intent(in) :: x
        integer, dimension(8), intent(in) :: atm
        real*8, intent(out) :: symdih
        real*8, dimension(:), intent(inout) :: gradient

        integer :: i, j, k
        real*8 :: dih1, dih2
        real*8, dimension(3, 4) :: grad1, grad2

        call CV_dihedral_aux(x, atm(1:4), dih1, grad1)
        call CV_dihedral_aux(x, atm(5:8), dih2, grad2)

        symdih = dih1 + dih2
        do i = 1, 4
            j = atm(i)*3
            gradient(j-2:j) = gradient(j-2:j) + grad1(:,i)
            k = atm(i+4)*3
            gradient(k-2:k) = gradient(k-2:k) + grad2(:,i)
        end do

    end subroutine CV_symdih
    !-------------------------



    !-------------------------
    subroutine CV_asymdih( x, atm, asymdih, gradient )

        real*8, dimension(:), intent(in) :: x
        integer, dimension(8), intent(in) :: atm
        real*8, intent(out) :: asymdih
        real*8, dimension(:), intent(inout) :: gradient

        integer :: i, j, k
        real*8 :: dih1, dih2
        real*8, dimension(3, 4) :: grad1, grad2

        call CV_dihedral_aux(x, atm(1:4), dih1, grad1)
        call CV_dihedral_aux(x, atm(5:8), dih2, grad2)

        asymdih = dih1 - dih2
        do i = 1, 4
            j = atm(i)*3
            gradient(j-2:j) = gradient(j-2:j) + grad1(:,i)
            k = atm(i+4)*3
            gradient(k-2:k) = gradient(k-2:k) - grad2(:,i)
        end do

    end subroutine CV_asymdih
    !-------------------------
    
    
    
    !-------------------------
    subroutine CV_dihedral_aux( x, atm, dihedral, grad ) !Taken from fdynamo 2.2 with minor changes

        real*8, dimension(:), intent(in) :: x
        integer, dimension(4), intent(in) :: atm
        real*8, intent(out) :: dihedral
        real*8, dimension(3,4), intent(out) :: grad
        real*8, dimension(3,4) :: xatm
        integer :: i

        real*8 :: df, diff, cosang, factij, factkl, msize, nsize, rkj, sgnfac
        real*8, dimension(3) :: drij, drkj, drkl, dti, dtj, dtk, dtl, m, n

        do i = 1, 4
            xatm(:,i) = x(atm(i)*3-2:atm(i)*3)
        end do      

        drij = xatm(:,1) - xatm(:,2)
        drkj = xatm(:,3) - xatm(:,2)
        drkl = xatm(:,3) - xatm(:,4)
        call map_periodic( drij )
        call map_periodic( drkj )
        call map_periodic( drkl )

        m = cross_product( drij, drkj )
        n = cross_product( drkj, drkl )
      
      ! . get the size of m and n.
        msize = sqrt ( dot_product ( m, m ) )
        nsize = sqrt ( dot_product ( n, n ) )
      
      ! . normalize m and n.
        m = m / msize
        n = n / nsize
      
      ! . ensure cosang is between -1 and 1.
        cosang = dot_product ( m, n )
        if( cosang > 1._8  ) cosang =  1._8
        if( cosang < -1._8 ) cosang = -1._8
        
      ! . get the sign of the angle.
        sgnfac = 1.0_8
        if ( dot_product ( drij, n ) < 0.0_8 ) sgnfac = -1.0_8
      
      ! . determine the dihedral.
        dihedral = sgnfac * acos ( cosang )
        dihedral = dihedral*57.295780_8

      ! . calculate rkj.
        rkj = sqrt ( dot_product ( drkj, drkj ) )
      
      ! . calculate dti and dtl.
        dti =   rkj * m / msize
        dtl = - rkj * n / nsize
      
      ! . calculate some additional factors.
        factij = dot_product ( drij, drkj ) / ( rkj * rkj )
        factkl = dot_product ( drkl, drkj ) / ( rkj * rkj )
      
      ! . calculate dtj and dtk.
        dtj = dti * ( factij - 1.0_8 ) - factkl * dtl
        dtk = dtl * ( factkl - 1.0_8 ) - factij * dti
      
      ! . calculate the gradients.

        grad(:,1) = dti
        grad(:,2) = dtj
        grad(:,3) = dtk
        grad(:,4) = dtl
        grad = grad * 57.295780_8

    contains

    
        function cross_product ( a, b )

            real ( kind = 8 ), dimension(1:3) :: cross_product
            real ( kind = 8 ), dimension(1:3), intent(in) :: a, b

            cross_product(1) = a(2) * b(3) - a(3) * b(2)
            cross_product(2) = a(3) * b(1) - a(1) * b(3)
            cross_product(3) = a(1) * b(2) - a(2) * b(1)

        end function cross_product 
    
    
    end subroutine CV_dihedral_aux
    !-------------------------

    
    
    !-------------------------
    subroutine CV_pointplane( x, atm, pointplane, gradient )

        real*8, dimension(:), intent(in) :: x
        integer, dimension(4), intent(in) :: atm
        real*8, intent(out) :: pointplane
        real*8, dimension(:), intent(inout) :: gradient
        real*8, dimension(3,4) :: xatm

        real*8,dimension(3)    :: V, U, W                   !Vectores y puntos definicion
        real*8,dimension(3)    :: P                         !Puntos definicion
        integer                :: i

        !PLANO
        real*8,dimension(4)     :: Pl                        !Ecuacion del plano
        real*8                  :: E, F                      !Escalares auxiliares
        real*8,dimension(3,4,4) :: gPl                       !Gradiente de la Ec del plano
        real*8,dimension(3,4)   :: gE, gF, g                 !Gradientes

        
        do i = 1, 4
            xatm(:,i) = x(atm(i)*3-2:atm(i)*3)
        end do             
 
        !Punto y vectores del plano
        P = xatm(:,2)
        U = xatm(:,3) - xatm(:,2)
        V = xatm(:,4) - xatm(:,2)
        W = xatm(:,4) - xatm(:,3)
        call map_periodic(U)
        call map_periodic(V)
        call map_periodic(W)

        !Ecuacion del Plano
        Pl(1) = U(2)*V(3) - U(3)*V(2)
        Pl(2) = U(3)*V(1) - U(1)*V(3)
        Pl(3) = U(1)*V(2) - U(2)*V(1)
        Pl(4) = -dot_product( Pl(1:3), P )

        !  Pl = Pl * sgn(da)    
        E = Pl(1)*xatm(1,1) + Pl(2)*xatm(2,1) + Pl(3)*xatm(3,1) !Escalar Auxiliar
        F = Pl(1)**2 + Pl(2)**2 + Pl(3)**2


        pointplane =  ( E + Pl(4) ) / sqrt( F )
        
!Gradiente del plano
!A
        gPl(:,1,1) = .0_8     
        gPl(1,2,1) = .0_8     
        gPl(2,2,1) = -W(3)     
        gPl(3,2,1) =  W(2)     
        gPl(1,3,1) = .0_8     
        gPl(2,3,1) =  V(3)     
        gPl(3,3,1) = -V(2)     
        gPl(1,4,1) = .0_8     
        gPl(2,4,1) = -U(3)     
        gPl(3,4,1) =  U(2)     
!B  
        gPl(:,1,2) = .0_8     
        gPl(1,2,2) =  W(3)     
        gPl(2,2,2) = .0_8     
        gPl(3,2,2) = -W(1)     
        gPl(1,3,2) = -V(3)     
        gPl(2,3,2) = .0_8     
        gPl(3,3,2) =  V(1)     
        gPl(1,4,2) =  U(3)     
        gPl(2,4,2) = .0_8     
        gPl(3,4,2) = -U(1)     
!C  
        gPl(:,1,3) = .0_8     
        gPl(1,2,3) = -W(2)     
        gPl(2,2,3) =  W(1)     
        gPl(3,2,3) = .0_8     
        gPl(1,3,3) =  V(2)     
        gPl(2,3,3) = -V(1)     
        gPl(3,3,3) = .0_8     
        gPl(1,4,3) = -U(2)     
        gPl(2,4,3) =  U(1)     
        gPl(3,4,3) = .0_8     

!D  
        gPl(:,:,4) = P(1)*gPl(:,:,1) + P(2)*gPl(:,:,2) + P(3)*gPl(:,:,3) 
        gPl(1,2,4) = gPl(1,2,4) + Pl(1) 
        gPl(2,2,4) = gPl(2,2,4) + Pl(2) 
        gPl(3,2,4) = gPl(3,2,4) + Pl(3) 
        gPl(:,:,4) = -gPl(:,:,4)
     
        !Gradiente de E
        gE = 0._8
        gE(:,1) = Pl(1:3)    

        gE = gE + xatm(1,1)*gPl(:,:,1) + xatm(2,1)*gPl(:,:,2) + xatm(3,1)*gPl(:,:,3)

        !Matriz auxiliar K
        gF = 2*Pl(1)*gPl(:,:,1) + 2*Pl(2)*gPl(:,:,2) + 2*Pl(3)*gPl(:,:,3)

        !Gradiente
        g = ( gE + gPl(:,:,4) ) / sqrt( F ) - ( ( E + Pl(4) )*gF ) / ( 2*F*sqrt(F) ) 

        do i = 1, 4
            gradient(atm(i)*3-2:atm(i)*3) = gradient(atm(i)*3-2:atm(i)*3) + g(:,i)
        end do           
        
    end subroutine CV_pointplane
    !-------------------------
    


!========Local subroutines========

    
    !------------------------
    subroutine dist( A, B, value, A_gradient, B_gradient )

	real*8, dimension(3), intent(in) :: A, B
	real*8, intent(out) :: value
	real*8, dimension(3), intent(inout) :: A_gradient, B_gradient

	real*8, dimension(3) :: dx

	dx = A-B
	call map_periodic( dx )

	value = sqrt( sum( dx**2 ) )

	A_gradient = A_gradient + dx/value
	B_gradient = B_gradient - dx/value

    end subroutine dist
    !------------------------
    
    
    !----------------------------
    subroutine map_periodic( dx )
    !Makes the vector between two points as short as possible
    !by taking the closest periodic cell.

        real*8, dimension(3), intent(inout) :: dx
        integer :: i
       
	if (ifbox==0) return

        do i = 1, 3
            if( abs( dx(i) ) > pbc_box(i)*0.5_8 ) dx(i) = &
                    dx(i) - sign( pbc_box(i), dx(i) )*int( abs(dx(i)/pbc_box(i)) + 0.5 )
        end do    
    
    end subroutine map_periodic
    !----------------------------
    
    
    
end module CV_module_mod
