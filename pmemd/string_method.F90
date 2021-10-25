! <compile=optimized>
!===============================================================================
!                 The String method Module
!===============================================================================
! 
!    The implementation of adaptive string method
!    If you use this software in your work, please cite
!    J. Phys. Chem. A, 2017, 121, 51, 9764-9772
!
!   Public subroutines:
!
!   string_define                Sets up the calculation
!   string_calc                  Launchs the string method calculation
!
!
!    Kirill Zinovjev, 2016, Universitat de Valencia (Spain)
!    ziki@uv.es
!
!===============================================================================
module string_method_mod

#ifdef MPI

    use CV_utilities_mod
    use string_utilities_mod
    use splines_utilities_mod
    
    use gbl_constants_mod, only : KB
    use random_mod, only : amrand
    use mdin_ctrl_dat_mod, only : dt, temp0
    use parallel_dat_mod, only : pmemd_master_comm, mytaskid, &
                                 pmemd_comm_number, numgroups
    use prmtop_dat_mod, only : natom


#include "mpif.h"
    
    private
    
    public :: string_define, string_calc, string_defined, minimize
    
    save

    integer :: node, proc, nnodes !Simulation node, process idx, number of nodes
    integer :: preparation_steps !Number of steps for the preparation stage (default: 1ps)
    integer :: buffer_size !for statistics used to update parameters and calculate PMFs (default: 2ps)
    integer :: string_move_period !string is moved every string_move_period steps (default: 25fs)
    integer :: step !Current step of string dynamics
    integer :: stop_step !step at which to stop the movement of the string

    real*8  :: force_scale !For gradual increase of forces during preparation
    real*8  :: gamma, position_gamma, force_gamma, force_kappa !friction
    real*8  :: Mav_damp !damping coefficient for Mav averaging
    real*8  :: RT !R*temperature
    real*8  :: K_d !Force constant for the d coordinate
    real*8  :: string_length !Arc-length of the string (using cubic splines)
    real*8  :: dpos, dK !displacement of the node position and the force constant
    real*8  :: mean_dx, mean_sigma2 !DEBUG

    logical :: string_defined = .false.
    logical :: server !Whether the current node is the server
    logical :: terminal !Whether the current node is at the end of the string

    logical :: string_move !Whether the string is allowed to move
    logical :: fix_ends    !Whether the terminal nodes are kept fixed
    logical :: write_M     !Whether Minv array should be written to the string file
    logical :: read_M       !Whether Minv are provided. If yes, no Minv update is performed
    logical :: generate_pos !Whether positions of nodes should be generated from the initial guess
    logical :: rescale_forces !Whether the force constants should be rescaled every time the string length changes
    
    !Note: in all the arrays that contain a set of elements for each node (ones that have one
    !of the dimensions == nnodes), the ordering is by processor, not by node. This is done
    !to facilitate mpi broadcasting, since the ordering of nodes changes when replica exchange
    !is used. The only exception is the string array, since the ordering by node in this case
    !is essential for cubic splines interpolation to work. !TODO: change comment
    real*8, dimension(:), allocatable :: K_l, pos !Force constants and positions
    real*8, dimension(:), allocatable :: l_buffer !Buffer with values of l coordinate !TODO: check if really needed
    real*8, dimension(:), allocatable :: weight_buffer !Buffer with weights for reweighting to path CV
    real*8, dimension(:,:), allocatable :: string !Current string
    real*8, dimension(:,:), allocatable :: CV_buffer !Buffer for sampled points in CV space
    real*8, dimension(:,:), allocatable :: n !n vectors (tangents to the string)
    real*8, dimension(:,:), allocatable :: decomp_coef !coefficients for PMF decomposition
    real*8, dimension(:,:), allocatable :: CV !Current values of CVs at all nodes
    real*8, dimension(:), allocatable :: dz !Displacement vector for the string node
    
    !Since all the matrices are symmetric, they are stored as one-dimensional
    !arrays of size msize=n*(n+1)/2 to save memory. Function matmulp is used instead of
    !matmul to work with this representation.
    real*8, dimension(:,:), allocatable :: Mav, Minv, B !Average M, inv(Mav), biasing matrices
    
    real*8, dimension(:,:,:), allocatable :: string_spline !Cubic splines interpolation of the string
    real*8, dimension(:,:,:), allocatable :: decomp_coef_spline !Cubic splines interpolation for PMF decomposition
    real*8, dimension(:,:,:), allocatable :: string_spline_smooth !Same, but for smooth normal vectors

    !---For output---
    integer :: output_period !Data output period (convergence, restart file) (timesteps)
    integer :: dat_unit !Unit for output data writing
    character*20 :: CV_frmt !Format for formatted printing of CV vectors
    character*200 :: dir !Directory for output data and data shared by all jobs
    !----------------
    
    !---For Replica Exchange---
    integer :: REX_period !time between exchange attempts (default: 25fs)
    integer, dimension(:), allocatable :: node_to_proc !Mapping of node indices to replica indices
    integer, dimension(:), allocatable :: proc_to_node !Mapping of replica indices to node indices
    integer, dimension(:), allocatable :: REX_hist !Histogram of exchanges
    !--------------------------
    
    !---For MPI---
    integer :: ierr, commmaster, sanderrank
    integer, dimension(mpi_status_size) :: status
    integer, dimension(:), allocatable :: rcv_size !For mpi_allgatherv
    !-------------
    
    !---For PMF and path CV---
    integer :: points_per_node !How many path points should be interpolated between every two nodes
    integer :: points_extra, bins_extra !Number of extra nodes and extra points for path extrapolation
    integer :: npoints !number of points along the string (without extrapolation)
    integer :: npointstotal !Total number of points (including extrapolation)
    integer :: nbins !Number of bins for Umbrella Integration
    real*8  :: lambda !Inverse distance between two adjacent nodes
    logical :: remove_z_bias !Whether the bias on the orthogonal direction should be removed for path CV
    logical :: reweight_z
    real*8, dimension(:), allocatable :: arc !arc length upto each node (for path CV)
    real*8, dimension(:), allocatable :: force, PMF !mean force and PMF along the path (using path CV)
    real*8, dimension(:,:), allocatable :: path !set of equidistant points along the string (with extrapolation)
    real*8, dimension(:,:), allocatable :: Minv_point !Metric tensor at each point
    !-------------------------

    !---For minimization---
    logical :: minimize !Whether the minimization (not dynamics) is performed
    integer :: nodes !Number of string nodes, node to take bias parameters from
    !----------------------

!==============================================================================
contains
!==============================================================================


    
    !==================================================================
    subroutine string_define(x_)
        
        real*8, dimension(3,natom), intent(in) :: x_
        real*8, dimension(natom*3) :: x

        integer :: i, j, unit
        real*8  :: force_constant, force_constant_d
        character*200 :: CV_file, guess_file, params_file
        logical :: only_PMF
        
    
        namelist /STRINGSETUP/ dir,&
                            step,&
                            output_period,&
                            preparation_steps,&
                            CV_file,&
                            guess_file,&
                            force_constant,&
                            force_constant_d,&
                            fix_ends,&
                            nbins,&
                            points_extra,&
                            REX_period,&
                            buffer_size,&
                            string_move,&
                            only_PMF,&
                            string_move_period,&
                            gamma,&
                            position_gamma,&
                            force_gamma,&
                            force_kappa,&
                            Mav_damp,&
                            write_M,&
                            read_M,&
                            params_file,&
                            remove_z_bias,&
                            reweight_z,&
                            points_per_node,&
                            nbins,&
                            minimize,&
                            nodes,&
                            node

        inquire(file = "STRING", exist=string_defined)
        if (.not. string_defined) return

        x = reshape(x_, (/ natom*3 /))

        ! To minimize changes, just assign to the old sander variables
        ! the corresponding values from pmemd
        commmaster = pmemd_master_comm
        sanderrank = mytaskid
        proc = pmemd_comm_number
        nnodes = numgroups

        if( sanderrank > 0 ) return

        RT = KB*temp0

        !Default values
        step = 0
        output_period = 100
        preparation_steps =  int(1._8/dt) !1 ps
        force_constant = -1
        force_constant_d = -1
        REX_period = int(0.05_8/dt) !50 fs
        buffer_size = int(2._8/dt) !2 ps
        gamma = 2D3
        position_gamma = 2D5
        force_gamma = 5
        force_kappa = 1000
        Mav_damp = 1D-3
        write_M = .false.
        read_M = .false.
        generate_pos = .true.
        rescale_forces = .true.
        only_PMF = .false.
        minimize = .false.
        nodes = 0
        node = 0
        guess_file = ""
        params_file = ""
        CV_file="CVs"

        !---PMF default values---
        nbins = 100
        points_extra = -1
        remove_z_bias = .true.
        reweight_z = .true.
        points_per_node = -1
        !------------------------
        
        !---String movement default values---
        fix_ends = .true.
        string_move = .true.
        string_move_period = 1
        !------------------------------------

        !Open the STRING file and read the parameters
        unit = next_unit()
        open(unit = unit, file = "STRING", status = "old")
        read(unit,NML=STRINGSETUP)
        close(unit)

        gamma = gamma*2.3901D-3 !Conversion from ps-1 to kcal/mol/a.m.u/A^2*ps
        position_gamma = position_gamma*2.3901D-3
    
        !Read CV definitions and calculate the initial values
        call prepare_CVs(trim(adjustl(CV_file)))

        if (minimize) then
            proc = 0
            if (nodes == 0 .or. node == 0) call write_error("nodes and node must be provided for minimization")
            nnodes = nodes
        end if
        proc = proc + 1
        if (.not. minimize) node = proc
        if (points_per_node == -1) points_per_node = 100/nnodes+1
        if (points_extra == -1) points_extra = points_per_node*5

        allocate(string(nCV,nnodes), CV(nCV, nnodes))
        allocate(K_l(nnodes), pos(nnodes), dz(nCV))
        allocate(node_to_proc(nnodes), proc_to_node(nnodes), REX_hist(nnodes-1), rcv_size(nnodes))
        allocate(l_buffer(buffer_size), weight_buffer(buffer_size), CV_buffer(nCV,buffer_size))
        allocate (n(nCV,nnodes), decomp_coef(nCV,nnodes), Mav(msize,nnodes), Minv(msize,nnodes), B(msize, nnodes))
        allocate(decomp_coef_spline(nnodes/2-1,5,nCV), string_spline(nnodes-1,5,nCV), string_spline_smooth(nnodes/2-1,5,nCV))

        REX_hist = 0
        node_to_proc = (/(i, i=0,nnodes-1)/) !MPI numbers processors from 0
        proc_to_node = node_to_proc
        pos = -2
        call update_CV(x)
        
        if (only_PMF .or. minimize) then
            string_move = .false.
            preparation_steps = 0
            read_M = .true.
            rescale_forces = .false.
            if (remove_z_bias) force_constant_d = 0
        end if
        
        if ((preparation_steps > 0) .and. (step == 0)) step = -preparation_steps

        if (.not. read_M) then
            Mav(:,node) = M
            call matinv(Mav(:,node), Minv(:,node))
        end if

        server = proc == 1
        terminal = (node == 1) .or. (node == nnodes)
        
        !Read the initial guess and interpolate nnodes points along the initial guess
        call read_initial_guess

        call string_PMF_prepare
        call reparameterize_linear
        if (.not. string_move) call update_path
        !call reparameterize_splines

        K_d = 0._8
        if (.not. (only_PMF .or. minimize)) then
            K_l = force_constant
            if (force_constant < 0) K_l = RT/((string_length/(nnodes-1))**2*0.25_8)
            K_d = K_l(node)*0.5_8 !By default K_d is half of K_l    
        end if
        if (force_constant_d > 0) K_d = force_constant_d
        force_scale = 1._8
        if (.not. minimize) call update_B

        if (trim(params_file) /= "") call read_params

        dat_unit = -1
        call assign_dat_file

        write(CV_frmt,*) nCV
        CV_frmt = "("//trim(adjustl(CV_frmt))//"E15.5E2)"
        
        if (server) then
            call write_string(0)
            call write_params
        end if

        write(6,*)
        write(6,"(A)")         "---String method------------------"
        write(6,"(A20,I15)")   "Number of CVs    ", nCV
        write(6,"(A20,I15)")   "Number of nodes  ", nnodes
        write(6,"(A20,I15)")   "Node             ", node
        write(6,"(A)")         "----------------------------------"
        write(6,*)
        write(6,*)

    contains
    
        !--------------------------------
        subroutine read_initial_guess
        
            integer :: i, u, ninit
            real*8, dimension(:,:), allocatable :: string_init
            real*8, dimension(msize) :: Mtmp, Mtmpinv
        
            u=next_unit()
            if (only_PMF .or. minimize) then
                open(unit=u, file="0_final.string", status="old")
                read(u,*) string
                read(u,*) Minv
                open(unit=u, file="final_parameters.dat", status="old")
                do i = 1, nnodes
                    read(u,*) pos(i), K_l(i)
                end do
                close(u)
                return
            end if
        
            if (trim(guess_file) == "") then
                allocate(string_init(nCV,nnodes))
                string_init = CV
                do i = 1, nnodes
                    call to_box(string_init(:,i))
                end do
            else
                open(unit=u, file=trim(adjustl(guess_file)), status="old")
                read(u,*) ninit
                allocate(string_init(nCV,ninit))
                read(u,*) string_init
                if (read_M) then
                    read(u,*) Minv
                    call matinv(Minv(:,node), Mav(:,node))
                end if
                close(u)
            end if
            
            if (ninit /= nnodes) then
                call mpi_allreduce(Mav(:,node), Mtmp, msize, mpi_real8, mpi_sum, commmaster, ierr)
                Mtmp = Mtmp/nnodes
                call matinv(Mtmp, Mtmpinv)
                call interpolate_linear(string_init, string, Mtmpinv)
            else
                string = string_init
            end if
            deallocate(string_init)
        
        end subroutine read_initial_guess
        !--------------------------------
        
        
        !--------------------------------
        subroutine interpolate_linear(A, B, M)
    
            real*8, dimension(:,:), intent(inout) :: A
            real*8, dimension(:), intent(in) :: M
            real*8, dimension(:,:), intent(out) :: B
        
            integer :: ninit, nfinal
            real*8  :: Ltotal, Lnew
            real*8, dimension(size(A, dim=2)) :: L
            real*8, dimension(size(A, dim=2),size(A, dim=2)) :: A_cont !continuous string (no jumps because of periodic coordinates)
            
            ninit = size(A, dim=2)
            nfinal = size(B, dim=2)
            
            call to_continuous(A)
        
            B(:,1) = A(:,1)
            B(:,nfinal) = A(:,ninit)
            
            !---Build array of string arc lengths up to each node---
            L(1) = 0._8
            do i = 2, ninit
                L(i) = L(i-1) + CV_distance(A(:,i-1), A(:,i), M)
            end do
            Ltotal = L(ninit)
            !-------------------------------------------------------            
            
            !---Obtain new nodes by simple linear interpolation---
            j = 2
            do i = 2, nfinal-1
                if (pos(1) < -1) then    
                    Lnew = Ltotal*(i-1)/(nfinal-1) !string arc length up to the new node i
                else
                    Lnew = pos(i)*Ltotal/pos(nnodes)
                end if
                do while( L(j) < Lnew )
                  j = j + 1
                end do
                B(:,i) = A(:,j-1) + map_periodic(A(:,j) - A(:,j-1)) * &
                                    (Lnew - L(j-1)) / (L(j) - L(j-1))
            end do
            !-----------------------------------------------------
    
        end subroutine interpolate_linear
        !--------------------------------


        subroutine read_params

            integer :: u
            u = next_unit()
            open(unit=u, file=trim(adjustl(params_file)), status="old")
            read(u,*) pos, K_l
            close(u)
            call reparameterize_linear

        end subroutine read_params
  
  
    end subroutine string_define
    !==================================================================
    
    
    
    !==================================================================
    subroutine assign_dat_file
    
        logical :: exist
        character*80 :: snode, ext
    
        write(snode,*) node
        snode = adjustl(snode)
        
        ext = ".dat"
        if (.not. string_move) then
            ext = "_final.dat"
            inquire(file=trim(dir)//trim(snode)//trim(ext), exist=exist)
        end if
        if (dat_unit > -1) close(dat_unit)
        
        if (.not. minimize) call mpi_barrier(commmaster, ierr)
        dat_unit = next_unit()
        open(unit=dat_unit, file=trim(dir)//trim(snode)//trim(ext), access="append")
        if (.not. string_move .and. .not. exist) write(dat_unit,"(4E15.5E2)") K_l(node), pos(node), K_d, 0._8 !In case WHAM will be used
    
    end subroutine assign_dat_file
    !==================================================================
    
    
    
    !==================================================================
    subroutine string_calc(x_, force, energy)
    
        real*8, dimension(3,natom), intent(in) :: x_
        real*8, dimension(:,:), intent(out) :: force
        real*8, intent(out) :: energy

        real*8, dimension(natom*3) :: x
        real*8 :: dpos_tmp, dK_tmp, pos_target, delta, sigma2, sigma2_target
        real*8 :: s, z
        real*8, dimension(nCV) :: dz_tmp

        logical :: is_move_step, is_REX_step

        energy = 0
        if (.not. string_defined) return
        if( sanderrank > 0 ) return

        x = reshape(x_, (/ natom*3 /))

        step = step + 1

        is_move_step = string_move .and. mod(step, string_move_period) == 0
        is_REX_step = mod(step, REX_period) == 0

        if (is_move_step .or. is_REX_step) then
            call update_CV(x)
        else
            call update_CVs(x)
        end if

        !Forces are applied during both preparation and production
        force = 0.
        if (string_move) then
            call add_force_ld
        else !If string is not moving, the Umbrella Sampling along the path CV is performed
            call add_force_sz
        end if

        if (minimize) then
            dz_tmp = 0
            call write_dat
            return
        end if

        !Update the Mav and Minv matrices
        if (string_move) then
            if (.not. read_M) then
                Mav(:,node) = Mav(:,node)*(1-Mav_damp) + M*Mav_damp
                call matinv(Mav(:,node), Minv(:,node))
                call normalize_M(n(:,node), Minv(:,node))
            end if
        end if        

        if ((step .le. 0) .and. string_move) then !Preparation
        
            !Force is gradually increased from 0 to the target value during the preparation
            force_scale = real(step, kind=8)/preparation_steps+1
            
        else !Production

            call fill_buffers

            dz_tmp = K_d*(map_periodic(CVs-string(:,node)) - &
                          n(:,node)*&
                          dot_product_M(map_periodic(CVs-string(:,node)), n(:,node), Minv(:,node)))
            delta = string_length/(nnodes-1)
            pos_target = delta*(node-1)-pos(node)
            sigma2_target = delta**2*0.25_8
            dpos_tmp = (pos_target-dot_product_M(map_periodic(CVs-string(:,node)),&
                        n(:,node),&
                        Minv(:,node)))

            if (step == 1) then
                dz = 0._8
                dpos = 0._8
                dK = 0._8
                !mean_dx = dpos_tmp
                mean_dx = 0._8
                mean_sigma2 = dpos_tmp**2
            end if
            
            mean_dx = mean_dx*0.99+dpos_tmp*0.01
            mean_sigma2 = mean_sigma2*0.99+dpos_tmp**2*0.01
            sigma2 = dpos_tmp**2
            if (step < 1000) then
                dK_tmp = 0.!Wait for mean_sigma2 and mean_dx to converge
            else    
                dK_tmp = RT/sigma2_target - RT/(mean_sigma2) + force_kappa*mean_dx**2!The last term to make K slightly larger when far
            end if
            dpos_tmp = dpos_tmp*K_l(node)
            
            if (string_move) then
                dz = dz + dz_tmp
                dpos = dpos + dpos_tmp
                dK = dK + dK_tmp
            end if

            if (is_move_step) then
                if (fix_ends .and. terminal) dz = 0._8
                string(:,node) = string(:,node) + dz/gamma*dt/string_move_period
                if (.not. terminal) pos(node) = pos(node) + scale_dpos(dpos/position_gamma*dt/string_move_period)
                K_l(node) = K_l(node) + dK/force_gamma*dt/string_move_period
                rcv_size = 1
                call mpi_allgatherv(K_l(node), 1, mpi_real8, K_l, &
                                    rcv_size, proc_to_node, mpi_real8, commmaster, ierr)
                dz = 0._8
                dpos = 0._8
                dK = 0._8
                call reparameterize_linear
            end if

            call write_dat

            if (is_REX_step) call REX
            if (mod(step, buffer_size) == 0) call string_PMF_calculate

            if (string_move .and. &
                server .and. &
                mod(step, output_period) == 0) then
                call write_string
                call write_convergence
                call write_params
            end if

            if (is_move_step) then
                call stop_string !stops the string movement if requested (by creating file STOP_STRING)
            end if

        end if


        
    contains
    
    
        real*8 function scale_dpos(x) !Ensures that the node positions never cross
        
            real*8, intent(in) :: x
            
            real*8 :: d
            
            d = merge(pos(node+1)-pos(node), pos(node)-pos(node-1), x>0)
            scale_dpos = x*exp(-x**2/d**2*4)
        
        end function scale_dpos
        
    
    
        !--------------------------
        subroutine write_dat
            call to_period(string)
            if (.not. string_move) write(dat_unit,"(2E15.5E2)",advance="no") s, z
            write(dat_unit,CV_frmt,advance="no") CVs
            write(dat_unit,CV_frmt,advance="no") string(:,node)
            write(dat_unit,CV_frmt) dz_tmp/gamma
            flush(dat_unit)
        end subroutine write_dat
        !--------------------------
        
        
        !--------------------------
        subroutine stop_string
        
            integer :: u, tmp
            integer :: av_step
            logical :: exist
        
            if (server) then
                stop_step = -1
                inquire(file=trim(dir)//"STOP_STRING", exist=exist)
                if (exist) then
                    u = next_unit()
                    open(u, file=trim(dir)//"STOP_STRING", status="old")
                    read(u,*,iostat=tmp) av_step, stop_step
                    close(u)
                end if
            end if
            
            call mpi_bcast(stop_step, 1, mpi_integer, 0, commmaster, ierr)
            call mpi_barrier(commmaster, ierr)
            call mpi_bcast(av_step, 1, mpi_integer, 0, commmaster, ierr)
            
            if ((stop_step > -1) .and. (step .ge. stop_step)) then
                string_move = .false.
                step = 0
                call average_string(av_step, stop_step) !average of string and parameters
                                                        !is taken over [av_step:stop_step]
                if (remove_z_bias) K_d = 0._8
                call assign_dat_file
                REX_hist = 0
            end if
        
        end subroutine stop_string
        !--------------------------`
        
        
        !--------------------------
        subroutine average_string(first, last)
        
            integer, intent(in) :: first, last
            
            integer :: i, u, n
            real*8, dimension(nnodes) :: params_tmp
            real*8, dimension(nCV) :: reference
            real*8, dimension(nCV, nnodes) :: string_tmp
            character*200 :: si

            reference = string(:, 1)
            string = 0._8
            n = (last-first)/output_period+1
            
            u = next_unit()
            do i = first, last, output_period
                write(si,*) i
                si = adjustl(si)
                open(unit=u, file=trim(dir)//trim(si)//".string", status="old")
                read(u,*) string_tmp
                string_tmp(:, 1) = reference + &
                        map_periodic(string_tmp(:, 1) - reference)
                call to_continuous(string_tmp)
                string = string + string_tmp
            end do
            string = string/n
            
            pos = 0._8
            open(unit=u, file=trim(dir)//"node_positions.dat", status="old")
            do i = 0, first-1, output_period
                read(u,*)
            end do
            do i = first, last, output_period
                read(u,*) params_tmp
                pos = pos + params_tmp
            end do
            pos = pos/n

            K_l = 0._8
            open(unit=u, file=trim(dir)//"force_constants.dat", status="old")
            do i = 0, first-1, output_period
                read(u,*)
            end do
            do i = first, last, output_period
                read(u,*) params_tmp
                K_l = K_l + params_tmp
            end do
            close(u)
            K_l = K_l/n

            call reparameterize_linear
            call write_string(step, "final", .true.)
            
            u = next_unit()
            open(unit=u, file=trim(dir)//"final_parameters.dat", status="replace")
            do i = 1, nnodes
                write(u,"(2F15.5)") pos(i), K_l(i)
            end do
            close(u)
        
        end subroutine average_string
        !--------------------------
    
        
        !--------------------------
        subroutine add_force_ld !Elliptic potential (for the string method itself)

            integer :: i, j
            real*8, dimension(nCV) :: gradld !Gradient of l+d in CV space

            gradld = matmulp(B(:,node), map_periodic(CVs-string(:,node)))

            do i = 1, n_used_atoms
                j = used_atoms(i)
                force(:,j) = - force_scale * matmul(Jacobian(j*3-2:j*3,:), gradld)
            end do

        end subroutine add_force_ld
        !--------------------------
        
        
        !--------------------------
        subroutine add_force_sz !Potential along path CV (for precise PMF calculation)

            integer :: i, j
            real*8, dimension(3,n_used_atoms) :: grad_s, grad_z, grad
            
            call get_s(CVs, s, z, Jacobian, grad_s, grad_z)
            energy = K_l(node)*0.5*(s-pos(node))**2 + K_d*0.5*merge(z, 0._8, z>0._8)**2

            do i = 1, n_used_atoms
                j = used_atoms(i)
                force(:,j) = - force_scale * &
                               (K_l(node)*(s-pos(node))*grad_s(:,i) + &
                                K_d*merge(z, 0._8, z>0._8)*grad_z(:,i))
            end do
        
        end subroutine add_force_sz
        !--------------------------
        
        
        !--------------------------
        subroutine fill_buffers
        
            integer :: idx
            
            idx = mod(step, buffer_size)+1
            
            l_buffer(idx) = dot_product_M(map_periodic(CVs-string(:,node)),&
                                          n(:,node),&
                                          Minv(:,node))
            weight_buffer(idx) = CV_distance2(CVs, &
                                              string(:,node), &
                                              B(:,node))*0.5_8
            CV_buffer(:,idx) = CVs
            
        end subroutine fill_buffers
        !--------------------------
    
    end subroutine string_calc
    !==================================================================
    
    
    
    !==================================================================
    subroutine write_string(alt_step, postfix, write_s)
    
        integer, intent(in), optional :: alt_step
        character(len=*), intent(in), optional :: postfix
        logical, intent(in), optional :: write_s
    
        character*200 :: filename
        character*20 :: frmt
        
        integer :: u, i, j, k
        
        if (present(alt_step)) then
            write(filename,*) alt_step
        else
            write(filename,*) step
        end if
        filename = adjustl(filename)
        
        if (present(postfix)) filename = trim(filename)//"_"//trim(adjustl(postfix))
        
        u = next_unit()
        open(unit=u, file=trim(dir)//trim(filename)//".string", status="replace")
        call to_continuous(string)
        write(u,CV_frmt) string
        call to_period(string)
        
        if (write_M .or. present(postfix)) then
            write(u,*)
            do i = 1, nnodes
                k = 0
                do j = 1, nCV
                    k = k+j
                    write(frmt,*) j
                    frmt = "("//trim(adjustl(frmt))//"E15.5E2)"
                    write(u,frmt) Minv(k-j+1:k,i)
                end do
                write(u,*)
            end do
        end if
        
         close(u)
         
         if(present(write_s) .and. write_s .eqv. .true.) then
            open(unit=u, file=trim(dir)//trim(filename)//"_CV.string", status="replace")
            write(frmt,*) nCV
            frmt = "(F15.5,"//trim(adjustl(frmt))//"E15.5E2)"
            call to_continuous(string)
            do i = 1, nnodes
                write(u,frmt) pos(i), string(:,i)
            end do
            call to_period(string)
            close(u)
         end if
        
    end subroutine write_string
    !==================================================================    
    
    
    
    !==================================================================    
    subroutine write_params

        integer :: u
        character*80 :: frmt
        
        write(frmt,*) nnodes
        frmt = "("//trim(adjustl(frmt))//"F15.5)"
        u = next_unit()
        open(unit=u, file=trim(dir)//"node_positions.dat", access="append")
        write(u,frmt) pos/pos(nnodes)
        open(unit=u, file=trim(dir)//"force_constants.dat", access="append")
        write(u,frmt) K_l
        close(u)
    
    end subroutine write_params
    !==================================================================    
    
    
    
    !==================================================================
    subroutine write_convergence
        
        integer :: i, u
        character*80 :: si, sstep
        real*8, dimension(nCV,nnodes) :: tmp_string
        
        write(sstep,*) step
        sstep = adjustl(sstep)
        u = next_unit()
        open(unit=u, file=trim(dir)//"convergence.dat", status="replace")

        do i = 0, step, output_period
            write(si,*) i
            si = adjustl(si)
            call read_string(trim(dir)//trim(si)//".string", tmp_string)
            write(u,"(I8,F15.5)") i, strings_distance(string, tmp_string, Minv)
        end do
        
        close(u)
        
    contains
        
        !----------------------------------------
        real*8 function strings_distance(A, B, M_list)
    
            real*8, dimension(:,:), intent(in) :: A, B
            real*8, dimension(:,:), intent(in) :: M_list
            
            integer :: i, nnodes
            
            nnodes = size(A, dim=2)
            strings_distance = 0._8
            do i = 1, nnodes
                strings_distance = strings_distance + &
                                   CV_distance(A(:,i), B(:,i), M_list(:,i))
            end do
            strings_distance = strings_distance/real(nnodes, kind=8)
    
        end function strings_distance
        !----------------------------------------
        
    end subroutine write_convergence
    !==================================================================
    
    
    
    !==================================================================
    subroutine update_B
    
        integer :: i, j, idx
    
        real*8, dimension(nCV) :: Minvn
        real*8, dimension(msize) :: S
    
        call normalize_M(n(:,node), Minv(:,node))
        Minvn = matmulp(Minv(:,node), n(:,node))
        
        idx = 0
        do i = 1, nCV
            do j = 1, i
                idx = idx + 1
                S(idx) = Minvn(i)*Minvn(j) !Outer product
            end do
        end do
        B(:,node) = S*(K_l(node)-K_d)+Minv(:,node)*K_d
        
        rcv_size = msize        
        call mpi_allgatherv(B(:,node), msize, mpi_real8, B, &
                           rcv_size, proc_to_node*msize, mpi_real8, commmaster, ierr)
    
    end subroutine update_B
    !==================================================================    
    
    
    
    !==================================================================    
    subroutine update_CV(x)
    
        real*8, dimension(natom*3), intent(in) :: x
        
        call update_CVs(x)
        CV(:,node) = CVs

        if (minimize) return

        rcv_size = nCV

        call mpi_allgatherv(CV(:,node), nCV, mpi_real8, CV, &
                           rcv_size, proc_to_node*nCV, mpi_real8, commmaster, ierr)        
    
    end subroutine update_CV
    !==================================================================    
    
    
    
    !==================================================================
    subroutine REX
    
        integer :: ex_proc
    
        real*8 :: rnd, dpos_tmp
        real*8, dimension(nCV) :: dz_tmp
        logical, dimension(nnodes) :: exchange
        real*8, dimension(buffer_size) :: buffer_tmp
        real*8, dimension(nCV,buffer_size) :: CV_buffer_tmp
        
        exchange = .false.

        if (mod(step/REX_period, 2) == mod(node, 2)) then !Triggers exchange mode: (1<->2 3<->4) or (1 2<->3 4)
            call amrand(rnd)
            exchange(node) = prob(node) .ge. rnd
            if ((.not. exchange(node)) .and. (node < nnodes) .and. &
                (pos(node+1)-pos(node) < string_length/(nnodes-1)*0.01_8)) then
                exchange = .true.
                write(*,*) "string_method: forced replica exchange at step",&
                           step, "between nodes", node, "and", node+1
                flush(6)
            end if
        end if
    
        rcv_size = 1
        call mpi_allgatherv(exchange(node), 1, mpi_logical, exchange, &
                           rcv_size, proc_to_node, mpi_logical, commmaster, ierr)
        rcv_size = msize
        call mpi_allgatherv(Mav(:,node), msize, mpi_real8, Mav, &
                           rcv_size, proc_to_node*msize, mpi_real8, commmaster, ierr)
                           
        close(dat_unit)
        if (exchange(node)) then

            ex_proc = node_to_proc(node+1)

            call mpi_send(l_buffer, buffer_size, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(l_buffer, buffer_size, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)            
            call mpi_send(weight_buffer, buffer_size, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(weight_buffer, buffer_size, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(CV_buffer, buffer_size*nCV, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(CV_buffer, buffer_size*nCV, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(dz, nCV, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(dz, nCV, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(dpos, 1, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(dpos, 1, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(mean_dx, 1, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(mean_dx, 1, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(mean_sigma2, 1, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(mean_sigma2, 1, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            
            node = node + 1
    
        else if (node > 1 .and. exchange(node-1)) then

            ex_proc = node_to_proc(node-1)

            call mpi_recv(buffer_tmp, buffer_size, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(l_buffer, buffer_size, mpi_real8, ex_proc, 0, commmaster, ierr)            
            l_buffer = buffer_tmp
            call mpi_recv(buffer_tmp, buffer_size, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(weight_buffer, buffer_size, mpi_real8, ex_proc, 0, commmaster, ierr)
            weight_buffer = buffer_tmp
            call mpi_recv(CV_buffer_tmp, buffer_size*nCV, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(CV_buffer, buffer_size*nCV, mpi_real8, ex_proc, 0, commmaster, ierr)
            CV_buffer = CV_buffer_tmp
            call mpi_recv(dz_tmp, nCV, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(dz, nCV, mpi_real8, ex_proc, 0, commmaster, ierr)
            call mpi_recv(dpos_tmp, 1, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(dpos, 1, mpi_real8, ex_proc, 0, commmaster, ierr)
            dpos = dpos_tmp
            call mpi_recv(dpos_tmp, 1, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(mean_dx, 1, mpi_real8, ex_proc, 0, commmaster, ierr)
            mean_dx = dpos_tmp
            call mpi_recv(dpos_tmp, 1, mpi_real8, ex_proc, mpi_any_tag, commmaster, status, ierr)
            call mpi_send(mean_sigma2, 1, mpi_real8, ex_proc, 0, commmaster, ierr)
            mean_sigma2 = dpos_tmp
            
            dz = dz_tmp
            node = node - 1
    
        end if

        terminal = (node == 1) .or. (node == nnodes)
        call update_maps    
        call assign_dat_file
        if (server) call REX_write_statistics

    contains
    
        !---------------------------------------
        real*8 function prob(i) !Calculates the exchange probability between nodes i and i+1
        
            integer, intent(in) :: i
            real*8 :: dE
            
            prob = 0._8
            if (i == nnodes) return
            dE = node_energy(i  , i+1)+&
                 node_energy(i+1, i  )-&
                 node_energy(i  , i  )-&
                 node_energy(i+1, i+1)
            prob = merge(exp(-dE/RT), 1._8, dE>0)
        
        end function prob
        !---------------------------------------
        
        
        !---------------------------------------
        real*8 function node_energy(node_idx, CV_idx)
    
            integer, intent(in) :: node_idx, CV_idx    

            node_energy = CV_distance2(CV(:,CV_idx), string(:,node_idx), B(:,node_idx))*0.5_8
            
        end function node_energy
        !---------------------------------------
        
        
        !---------------------------------------
        subroutine REX_write_statistics
        
            integer :: i, u
            character*200 :: frmt, filename
            
            u = next_unit()
            filename = "exchanges.REX"
            if (.not. string_move) filename = "exchanges_final.REX"
            open(unit=u, file=trim(dir)//trim(adjustl(filename)), access="append")
            do i = 1, nnodes
                if (exchange(i)) then
                    REX_hist(i) = REX_hist(i)+1
                    write(u,"(I10,I5)") step, i 
                end if
            end do
            close(u)
            
            filename = "plot.REX"
            if (.not. string_move) filename = "plot_final.REX"
            open(unit=u, file=trim(dir)//trim(adjustl(filename)), access="append")
            write(frmt,*) nnodes
            frmt = "(I10,"//trim(adjustl(frmt))//"I5)"
            write(u,frmt) step, proc_to_node+1
            close(u)
            
            filename = "histogram.REX"
            if (.not. string_move) filename = "histogram_final.REX"
            open(unit=u, file=trim(dir)//trim(adjustl(filename)), status="replace")
            write(u,"(I10)") REX_hist
            close(u)
        
        end subroutine REX_write_statistics
        !---------------------------------------
        
        
        !---------------------------------------
        subroutine update_maps
        
            integer :: i, tmp
            
            do i = 1, nnodes-1
                if (exchange(i)) then
                    tmp = node_to_proc(i)
                    node_to_proc(i) = node_to_proc(i+1)
                    node_to_proc(i+1) = tmp    
                end if
            end do
            
            do i = 1, nnodes
                proc_to_node(node_to_proc(i)+1) = i-1
            end do
        
        end subroutine update_maps
        !---------------------------------------
        
    end subroutine REX
    !==================================================================
    
    
    
    !==================================================================
    subroutine reparameterize_linear
    
        integer :: i, j
        real*8, dimension(nnodes) :: L
        real*8, dimension(nCV) :: dz

        if (minimize) then
            call to_continuous(string)
            !---Build array of string arc lengths up to each node---
            L(1) = 0._8
            do i = 2, nnodes
                L(i) = L(i-1) + CV_distance(string(:,i-1), string(:,i), (Minv(:,i-1)+Minv(:,i))*0.5_8)
            end do
            string_length = L(nnodes)
            !-------------------------------------------------------
            call CubicSplinesFitND(L, string, string_spline_smooth)
            call to_box(string(:,node))
            return
        end if

        rcv_size = msize
        call mpi_allgatherv(Minv(:,node), msize, mpi_real8, Minv, &
                           rcv_size, proc_to_node*msize, mpi_real8, commmaster, ierr)
        rcv_size = nCV
        call mpi_allgatherv(string(:,node), nCV, mpi_real8, string, rcv_size, &
                            proc_to_node*nCV, mpi_real8, commmaster, ierr)
        rcv_size = 1
        call mpi_allgatherv(pos(node), 1, mpi_real8, pos, rcv_size, &
                            proc_to_node, mpi_real8, commmaster, ierr)

        call to_continuous(string)


        if (rescale_forces) K_l = K_l*string_length**2
        !---Build array of string arc lengths up to each node---
        L(1) = 0._8
        do i = 2, nnodes
            L(i) = L(i-1) + CV_distance(string(:,i-1), string(:,i), (Minv(:,i-1)+Minv(:,i))*0.5_8)
        end do
        string_length = L(nnodes)
        if (generate_pos) then
            pos = L
            generate_pos = .false.
        end if
        !-------------------------------------------------------
        if (rescale_forces) K_l = K_l/string_length**2

        if (pos(1) < -1) then !If interpolating for the first time, take equidistant points
            pos(node) = string_length*(node-1)/(nnodes-1)
        else
            !Rescale positions to the new string length
            pos(node) = pos(node)*string_length/pos(nnodes)
        end if

        !---Obtain new nodes by simple linear interpolation---
        if (.not. terminal) then
            j = 2
            do while( L(j) < pos(node) )
                j = j + 1
            end do
            dz = string(:,j-1) + map_periodic(string(:,j) - string(:,j-1)) * &
                             (pos(node) - L(j-1)) / (L(j) - L(j-1)) - string(:,node)
            string(:,node) = string(:,node) + dz
        end if
        !-----------------------------------------------------

        rcv_size = 1
        call mpi_allgatherv(pos(node), 1, mpi_real8, pos, &
                           rcv_size, proc_to_node, mpi_real8, commmaster, ierr)
        rcv_size = nCV
        call mpi_allgatherv(string(:,node), nCV, mpi_real8, string, rcv_size, &
                            proc_to_node*nCV, mpi_real8, commmaster, ierr)
                        
        !Reinterpolate the string using arc-length as x
        !The interpolation is done with least squares cubic splines to make the normal
        !vectors more stable. Since normal vectors define the orientation of the bias,
        !potential, this is essential for the stability of the simulation.
        call CubicSplinesFitND(L, string, string_spline_smooth)        

        !Extract tangent vectors at reference positions along the string
        n(:,node) = spline_der_ND(pos(node), string_spline_smooth)
        call normalize_M(n(:,node), Minv(:,node))
        
        call update_decomp_coef(L)
        
        call update_path
        call update_B
        call to_box(string(:,node))
        
    end subroutine reparameterize_linear
    !==================================================================
    
    
    
    !==================================================================
    subroutine reparameterize_splines
    
        integer :: i
        real*8, dimension(nnodes) :: x !Initial parameterization (arbitrary, node indices are used here)
        real*8, dimension(nnodes) :: L !arc lengths up to each node
        real*8, dimension(nnodes) :: string_der !|dz/dx|_M - string derivative length with metric Minv
        real*8, dimension(nnodes-1,5) :: string_der_spline !|dz/dx|_M as function of x (to calculate the arc length)
        real*8, dimension(nnodes-1,5) :: xofL_spline !x as function of arc length
        
        rcv_size = msize
        call mpi_allgatherv(Minv(:,node), msize, mpi_real8, Minv, &
                           rcv_size, proc_to_node*msize, mpi_real8, commmaster, ierr)
        rcv_size = nCV
        call mpi_allgatherv(string(:,node), nCV, mpi_real8, string, rcv_size, &
                            proc_to_node*nCV, mpi_real8, commmaster, ierr)        
        rcv_size = 1
        call mpi_allgatherv(pos(node), 1, mpi_real8, pos, rcv_size, &
                            proc_to_node, mpi_real8, commmaster, ierr)        
                            
        call to_continuous(string)
        
        x = (/(i, i=0,nnodes-1)/)
        call CubicSplinesFitND(x, string, string_spline_smooth)
        
        K_l = K_l*string_length**2 !Force constants are adjusted for the stretching/compressing of the string
        
        !---Build array of string arc lengths up to each node---
        do i = 1, nnodes
            string_der(i) = len_M(spline_der_ND(x(i), string_spline_smooth), Minv(:,i))
        end do

        call CubicSplines(x, string_der, string_der_spline)
        L(1) = 0._8
        do i = 2, nnodes
            L(i) = L(i-1) + spline_integrate(x(i-1), x(i), string_der_spline)
        end do
        string_length = L(nnodes)
        !-------------------------------------------------------
        
        if (generate_pos) then !positions are generated from initial positions of point (if provided in the guess file)
            pos = L
            generate_pos = .false.
        end if
        K_l = K_l/string_length**2
        
        if (pos(1) < -1) then !If interpolating for the first time, take equidistant points
            pos(node) = string_length*(node-1)/(nnodes-1)
        else
            !Rescale positions to the new string length
            pos(node) = pos(node)*string_length/pos(nnodes)
        end if
        
        call CubicSplines(L, x, xofL_spline)
        string(:,node) = spline_value_ND(spline_value(pos(node), xofL_spline), string_spline_smooth)
        
        rcv_size = 1
        call mpi_allgatherv(pos(node), 1, mpi_real8, pos, &
                           rcv_size, proc_to_node, mpi_real8, commmaster, ierr)
        rcv_size = nCV
        call mpi_allgatherv(string(:,node), nCV, mpi_real8, string, rcv_size, &
                            proc_to_node*nCV, mpi_real8, commmaster, ierr)
                        
        !Reinterpolate the string using arc-length as x
        !The interpolation is done with least squares cubic splines to make the normal
        !vectos more stable. Since normal vectors define the orientation of the bias,
        !potential, this is essential for the stability of the simulation.
        call CubicSplinesFitND(L, string, string_spline_smooth)        

        !Extract tangent vectors at reference positions along the string
        n(:,node) = spline_der_ND(pos(node), string_spline_smooth)
        call normalize_M(n(:,node), Minv(:,node))
  
        call update_decomp_coef(L)
        
        !call update_path
        call update_B
        call to_box(string(:,node))        
    
    end subroutine reparameterize_splines
    !==================================================================
    
    
    
    !==================================================================
    subroutine update_decomp_coef(L)
    !Create splines interpolation for PMFCV decomposition
    
        real*8, dimension(:), intent(in) :: L !arc lengths up to each node
    
        call mpi_allgatherv(n(:,node), nCV, mpi_real8, n, rcv_size, &
                            proc_to_node*nCV, mpi_real8, commmaster, ierr)  
                            
        decomp_coef(:,node) = matmulp(Minv(:,node), n(:,node))*n(:,node)
        
        call mpi_allgatherv(decomp_coef(:,node), nCV, mpi_real8, decomp_coef, rcv_size, &
                            proc_to_node*nCV, mpi_real8, commmaster, ierr)    
                            
        call CubicSplinesFitND(L, decomp_coef, decomp_coef_spline)  
        
    end subroutine update_decomp_coef
    !==================================================================
    
    
!------------------------------------------------------------------------------!
!               Subroutines for path CV and PMF calculation                    !
!------------------------------------------------------------------------------!


    !==================================================================
    subroutine string_PMF_prepare

        npoints = (nnodes-1)*points_per_node+1
        npointstotal = npoints+points_extra*2
        bins_extra = (nbins-1)*points_extra/(npointstotal-1)
        nbins = nbins+bins_extra*2
    
        allocate(arc(npointstotal), force(nbins), PMF(nbins))
        allocate(path(nCV,npointstotal), Minv_point(msize,npointstotal))
    
    end subroutine string_PMF_prepare
    !==================================================================
    
    
    
    !==================================================================
    subroutine update_path
    
        integer :: i, j
        real*8, dimension(nCV) :: neg_extra, pos_extra !Vectors to extrapolate the path
    
        lambda = (npoints-1)/string_length
        arc = (/(i, i=-points_extra,npoints+points_extra-1)/)/lambda
    
        do i = points_extra+1, npoints+points_extra
            path(:,i) = spline_value_ND(arc(i), string_spline_smooth)
        end do
        
        !---Treat extrapolated points---
        neg_extra = path(:,points_extra+1)-path(:,points_extra+2)
        pos_extra = path(:,npoints+points_extra)-path(:,npoints+points_extra-1)
        do i = 1, points_extra
            path(:,points_extra-i+1) = path(:,points_extra-i+2)+neg_extra
            path(:,npoints+points_extra+i) = path(:,npoints+points_extra+i-1)+pos_extra
        end do
        !-------------------------------

        !---generate Minv at each point---
        do i = 1, nnodes-1
            do j = 0, points_per_node-1
                Minv_point(:,points_extra+(i-1)*points_per_node+j+1) = (Minv(:,i)*(points_per_node-j)+&
                                                                        Minv(:,i+1)*j)/points_per_node
            end do
        end do
        do i = 1, points_extra
            Minv_point(:,i) = Minv(:,1)
        end do
        do i = npoints+points_extra, npointstotal
            Minv_point(:,i) = Minv(:,nnodes)
        end do
        !---------------------------------
    
    end subroutine update_path
    !==================================================================
    


    !==================================================================
    subroutine string_PMF_calculate
    
        integer :: i, nbuf, nnodes, bin, minpos, bin_start, bin_end
        real*8 :: mean, sigma2, bin_width, mean_error, sigma2_error
        real*8 :: K, s0, tmp
        integer, dimension(nbins) :: N_bin !Number of points in each bin
        real*8, dimension(nbins) :: force, t_bin
        real*8, dimension(nbins) :: pop_bin !Bin populations
        real*8, dimension(nbins) :: w_error, force_error, PMF_error
        real*8, dimension(nCV) :: CVtmp
        real*8, dimension(:), allocatable :: s, z, w
        character*200 :: sstep, snode, filename

        if (string_move) then
            call update_path
            nbuf = buffer_size
        else
            close(dat_unit)
            write(snode,*) node
            snode = adjustl(snode)
            filename = trim(dir)//trim(snode)//"_final.dat"
            nbuf = step
            open(unit=dat_unit, file=filename, status="old")
            if (.not. string_move) read(dat_unit,*)
        end if

        allocate(s(nbuf), z(nbuf), w(nbuf))
        bin_width = (npointstotal-1)/lambda/(nbins-1)
        t_bin = (/(i, i=0,nbins-1)/)*bin_width + arc(1) !Lower bound of each bin

        do i = 1, nbuf
            if (string_move) then
                CVtmp = CV_buffer(:,i)
                call get_s(CVtmp, s(i), z(i))
            else
                read(dat_unit,*) s(i), z(i)
            end if
        end do
        if (.not. string_move) close(dat_unit)
        
        K = K_l(node)
        s0 = pos(node)
        do i = 1, nbuf
            !Reweigting the bias to be only along the s coordinate. In case of the string potential
            !the bias should be removed completely and than the bias along s be added. In case of
            !the path CV potential only the bias along z coordinate should be removed.
            w(i) = merge(exp((weight_buffer(i)-K*0.5_8*(s(i)-s0)**2)/RT), &
                         merge(exp(K_d*0.5_8*merge(z(i), 0._8, z(i)>0._8)**2/RT), & 
                               1._8, &
                               reweight_z),&
                         string_move)
        end do
        w = w/sum(w)
                
        !---Umbrella Integration---
        mean = dot_product(s, w)
        sigma2 = dot_product(s**2, w) - mean**2
        N_bin = 0
            
        pop_bin = 0._8    
        do i = 1, nbuf
            bin = int((s(i)-t_bin(1))/bin_width)+1
            bin = max(min(bin, nbins), 1)
            N_bin(bin) = N_bin(bin) + 1
        end do
        t_bin = t_bin+bin_width*0.5_8 !Middle of each bin
        pop_bin = exp(-(t_bin-mean)**2/(sigma2*2))/sqrt(sigma2*2*3.1416_8)*bin_width
        force = RT*(t_bin-mean)/sigma2 - K*(t_bin-s0)
        force = force*pop_bin
        w_error = pop_bin

        call mpi_allreduce(MPI_IN_PLACE, force, nbins, mpi_real8, mpi_sum, commmaster, ierr)
        call mpi_barrier(commmaster, ierr)

        call mpi_allreduce(MPI_IN_PLACE, pop_bin, nbins, mpi_real8, mpi_sum, commmaster, ierr)
        call mpi_barrier(commmaster, ierr)

        call mpi_allreduce(MPI_IN_PLACE, N_bin, nbins, mpi_integer, mpi_sum, commmaster, ierr)
        call mpi_barrier(commmaster, ierr)

        mean_error = sigma2/(nbuf/(0.05/dt))
        sigma2_error = 2*sigma2**2/(nbuf/(0.05/dt)-1)
        w_error = w_error/pop_bin
        force_error = (w_error/sigma2*(force*(1-w_error)*(t_bin-mean)-RT))**2*mean_error + &
                      (w_error/sigma2**2 * &
                        (t_bin-mean)*(force*(1-w_error)*(t_bin-mean)/2-RT))**2*sigma2_error
                        
        call mpi_allreduce(MPI_IN_PLACE, force_error, nbins, mpi_real8, mpi_sum, commmaster, ierr)
        call mpi_barrier(commmaster, ierr)
        
        if (server) then
        
            force = -force/pop_bin
            PMF(1) = 0._8
            do i = 2, nbins
                PMF(i) = PMF(i-1)-force(i-1)*bin_width
            end do
            
            bin_start = 1
            do while (N_bin(bin_start) == 0)
                bin_start = bin_start + 1
            end do
            bin_end = nbins
            do while (N_bin(bin_end) == 0)
                bin_end = bin_end - 1
            end do
            
            minpos = minloc(PMF(bins_extra:nbins/4), dim=1)+bins_extra-1
            PMF = PMF-PMF(minpos)
            PMF_error(minpos) = 0._8
            do i = minpos-1, 1, -1
                PMF_error(i) = PMF_error(i+1)+&
                (force_error(i+1)+force_error(i))*bin_width*0.5_8
            end do
            do i = minpos+1, nbins
                PMF_error(i) = PMF_error(i-1)+&
                (force_error(i-1)+force_error(i))*bin_width*0.5_8
            end do
            
            force_error = sqrt(force_error)*1.96_8
            PMF_error = sqrt(PMF_error)*1.96_8
            write(sstep,*) merge(step, nbuf, string_move)
            sstep = adjustl(sstep)
            if (.not. string_move) sstep = trim(sstep)//"_final"
            call string_PMF_write(trim(dir)//trim(sstep)//".PMF")
            call string_PMFCV_write(trim(dir)//trim(sstep)//"_CV")
        end if
        !--------------------------
        
        if (.not. string_move) call assign_dat_file
        
    contains
    
    
        !--------------------------
        subroutine string_PMFCV_write(filename)
            
                    character(len=*) :: filename
                    
                    integer :: i, uF, uA
                    
                    real*8, dimension(nCV) :: coef !Coefficients for the decomposition (Minv*n)
                    real*8, dimension(nCV) :: tmp
                    real*8, dimension(nCV,nbins) :: forceCV, PMFCV, forceCV_error, PMFCV_error
                    
                    character*80 :: si, frmt
                    
                    write(frmt,*) nCV*2+1
                    frmt = "("//trim(adjustl(frmt))//"F15.5)"

                    !Calculate forceCV and forceCV_error
                    do i = 1, nbins
                        coef = spline_value_ND(t_bin(i)+bin_width*0.5, decomp_coef_spline)
                        forceCV(:,i) = force(i)*coef
                        forceCV_error(:,i) = force_error(i)*coef
                    end do
                    
                    !Calculate PMFCV and PMFCVerror
                    PMFCV(:,1) = 0.
                    do i = 2, nbins
                        PMFCV(:,i) = PMFCV(:,i-1)-forceCV(:,i-1)*bin_width
                    end do
                    tmp = PMFCV(:,minpos)
                    do i = 1, nbins
                        PMFCV(:,i) = PMFCV(:,i)-tmp
                    end do
                    PMFCV_error(:,minpos) = 0._8
                    do i = minpos-1, 1, -1
                        PMFCV_error(:,i) = PMFCV_error(:,i+1)+&
                        (forceCV_error(:,i+1)**2+forceCV_error(:,i)**2)*bin_width*0.5_8
                    end do
                    do i = minpos+1, nbins
                        PMFCV_error(:,i) = PMFCV_error(:,i-1)+&
                        (forceCV_error(:,i-1)**2+forceCV_error(:,i)**2)*bin_width*0.5_8
                    end do
                    PMFCV_error = sqrt(PMFCV_error)*1.96_8
                    
                    !Write everything
                    uF = next_unit()
                    open(unit=uF, file=trim(filename//".force"), status="replace")                    
                    uA = next_unit()
                    open(unit=uA, file=trim(filename//".PMF"), status="replace")
                    write(uF,"(A15)",advance="no") "s"
                    write(uA,"(A15)",advance="no") "s"
                    do i = 1, nCV
                        write(si,*) i
                        si = adjustl(si)
                        write(uF,"(A15)",advance="no") "CV_"//trim(si)
                        write(uA,"(A15)",advance="no") "CV_"//trim(si)
                    end do
                    do i = 1, nCV
                        write(si,*) i
                        si = adjustl(si)
                        write(uF,"(A15)",advance="no") "CV_"//trim(si)//"_CI (95%)"
                        write(uA,"(A15)",advance="no") "CV_"//trim(si)//"_CI (95%)"
                    end do
                    write(uF,*)
                    write(uA,*)
                    do i = 1, nCV*2+1
                        write(uF,"(A15)",advance="no") "---------------"
                        write(uA,"(A15)",advance="no") "---------------"
                    end do
                    write(uF,*)
                    write(uA,*)                    
                    do i = bin_start, bin_end
                        write(uF,frmt) t_bin(i), forceCV(:,i), forceCV_error(:,i)
                        write(uA,frmt) t_bin(i), PMFCV(:,i), PMFCV_error(:,i)
                    end do
                    close(uF)
                    close(uA)
                    
        end subroutine
        !--------------------------
    
        
        !--------------------------
        subroutine string_PMF_write(filename)
        
            character(len=*) :: filename
            
            integer :: i, u
            
            u = next_unit()
            open(unit = u, file = trim(filename), status="replace")
            write(u,"(A15,A10,5A15)") "s", "N_bin", "Pop_bin", "Force", "Force_CI (95%)", "PMF", "PMF_CI (95%)"
            write(u,"(A)") "----------------------------------------------------------------------------------------------------"
            
            do i = bin_start, bin_end
                write(u,"(F15.5,I10,5F15.5)") t_bin(i),&
                                              N_bin(i),&
                                              pop_bin(i),&
                                              force(i),&
                                              force_error(i),&
                                              PMF(i),&
                                              PMF_error(i)
            end do
            
            close(u)
            
        end subroutine string_PMF_write
        !--------------------------
        
    end subroutine string_PMF_calculate
    !==================================================================
    
    
    
    !==================================================================
    subroutine get_s(CV, s, z, Jacobian, grad_s, grad_z)
    !Calculation of a metric-corrected path CV (J. Comput. Chem. 2014, 35, 1672–1681)

        real*8, dimension(nCV), intent(in) :: CV
        real*8, intent(out) :: s
        real*8, intent(out), optional :: z
        real*8, dimension(natom*3,nCV), intent(in), optional :: Jacobian
        real*8, dimension(3,n_used_atoms), intent(out), optional :: grad_s, grad_z

        integer :: i, j
        real*8 :: d, sumw
        
        real*8, dimension(npointstotal) :: w
        real*8, dimension(nCV) :: dCV
        real*8, dimension(nCV,npointstotal) :: grad_tmp !temporary values for gradient calculation

        if (present(grad_s) .neqv. present(Jacobian)) &
            call write_error("get_s: grad_s and Jacobian should be either both present or both not")
        
        do i = 1, npointstotal
            dCV = map_periodic(CV-path(:,i))
            d = len_M(dCV, Minv_point(:,i))
            dCV = matmulp(Minv_point(:,i), dCV)
            w(i) = exp(-lambda*d)
            if (present(grad_s)) grad_tmp(:,i) = dCV/d*w(i)
        end do
        
        sumw = sum(w)
        if (present(z)) z = -log(sumw)/lambda
        w = w/sumw
        grad_tmp = grad_tmp/sumw
        s = dot_product(arc, w)
        
        if (present(grad_z)) then
            do i = 1, n_used_atoms
                j = used_atoms(i)
                grad_z(:,i) = matmul(Jacobian(j*3-2:j*3,:), sum(grad_tmp, dim=2))
            end do
        end if
        
        if (present(grad_s)) then
            do i = 1, npointstotal
                grad_tmp(:,i) = grad_tmp(:,i)*(s-arc(i))
            end do
            do i = 1, n_used_atoms
                j = used_atoms(i)
                grad_s(:,i) = matmul(Jacobian(j*3-2:j*3,:), sum(grad_tmp, dim=2))*lambda
            end do
        end if
    
    end subroutine get_s
    !==================================================================


    
#endif    
    
end module string_method_mod
