#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "nfe-config.h"

!------------------------------------------------------------------------------
! sander: the Molecular Dynamics/NMR Refinement/Modeling Module of AMBER
!         (Assisted Model Building with Energy Refinement).  This main driver
!         routine will manage all operations of the program.
!------------------------------------------------------------------------------
subroutine sander()

  use constants, only : INV_AMBER_ELECTROSTATIC
  use lmod_driver
  use trace

  use state
#if !defined(DISABLE_NFE)
  use nfe_sander_hooks, only : &
      nfe_on_sander_init => on_sander_init, &
      nfe_on_sander_exit => on_sander_exit
  use nfe_sander_proxy, only : infe
#endif /* DISABLE_NFE */

  ! The main qmmm_struct contains all the QMMM variables and arrays
  use qmmm_read_and_alloc, only : read_qmmm_nm_and_alloc
  use qmmm_vsolv_module, only: qmmm_vsolv_store_parameters, new
  use qm2_extern_module, only: qm2_extern_finalize
#ifdef QUICK
  use quick_module, only: quick_finalize
#endif
#ifdef TCPB
  use tcpb_module, only: tcpb_finalize
#endif
#ifdef MPI
  use qmmm_module, only : qmmm_nml, qmmm_struct, deallocate_qmmm, qmmm_mpi, &
                          qm2_struct, qmmm_vsolv, qm2_params, qmmm_mpi_setup
  use sebomd_module, only : sebomd_obj, sebomd_open_files, sebomd_bcast_obj, &
                            sebomd_close_files, sebomd_setup
  use decomp, only : allocate_int_decomp, allocate_real_decomp, &
                     deallocate_int_decomp, deallocate_real_decomp, &
                     synchronize_dec, build_dec_mask, decmask, indx, jgroup, &
                     nat, nrs
#else
  use qmmm_module, only : qmmm_nml, qmmm_struct, deallocate_qmmm, qmmm_mpi, &
                          qm2_struct, qmmm_vsolv, qm2_params
  use sebomd_module, only : sebomd_obj, sebomd_open_files, &
                            sebomd_close_files, sebomd_setup
  use decomp, only : allocate_int_decomp, allocate_real_decomp, &
                     deallocate_int_decomp, deallocate_real_decomp, &
                     nat, nrs
#endif /* MPI */
  use sebomd_arrays, only : init_sebomd_arrays, cleanup_sebomd_arrays
#ifdef LES
  use genbornles
#  ifdef MPI
  use les_data, only : cnum
#  endif /* MPI */
#else
  use genborn
#endif /* LES */
  use les_data, only : temp0les
  use fastwt
  use relax_mat
  use nmr, only: nmrrad, impnum
  use ew_recip, only: deallocate_m1m2m3,first_pme
  use parms
  use molecule, only : mol_info, iwrap_mask_atoms, &
                       allocate_molecule, deallocate_molecule
  use nblist, only: cutoffnb, skinnb, nblist_allocate, nblist_deallocate, &
                    nblist_allreal, nblist_allint, num_calls_nblist, &
                    first_list_flag
  use stack
  use amoeba_runmd, only : AM_RUNMD_get_coords, AM_RUNMD
  use amoeba_mdin, only : beeman_integrator, iamoeba, am_nbead
  use amoeba_interface, only : AMOEBA_deallocate, AMOEBA_readparm

  use pol_gauss_mdin, only : ipgm
  use pol_gauss_interface, only : POL_GAUSS_readparm, POL_GAUSS_deallocate

#ifdef RISMSANDER
  use sander_rism_interface, only: rism_setparam, rism_init, rism_finalize
#endif

#ifdef PUPIL_SUPPORT
  use pupildata
#endif /* PUPIL */

#ifdef APBS
  use apbs
#endif /* APBS */

#ifdef MPI /* SOFT CORE */
  use softcore, only: setup_sc, cleanup_sc, ifsc, extra_atoms, sc_sync_x, &
                      summarize_ti_changes, sc_check_perturbed_molecules, &
                      ti_check_neutral, tishake
  use mbar, only: setup_mbar, cleanup_mbar, ifmbar
#endif

  ! For Linear Interaction Energy calculations
  use linear_response, only: ilrt, setup_linear_response, &
                              cleanup_linear_response

#if defined(MPI)
  use evb_parm, only: xch_type
#  if defined(LES)
  use evb_pimd, only: evb_pimd_init
#  endif /* LES */
  ! Replica Exchange Molecular Dynamics
  use remd, only : rem, mdloop, remd1d_setup, remd_exchange, reservoir_remd_exchange, &
                   remd_cleanup, hremd_exchange, ph_remd_exchange, e_remd_exchange, &
                   multid_remd_setup, multid_remd_exchange, setup_pv_correction, &
                   rremd, stagid
#  ifdef VERBOSE_REMD
  use remd, only : repnum
#  endif
  use bintraj, only: setup_remd_indices
#else
#  define rem 0
#endif /* MPI */

  use pimd_vars, only: ipimd
  use neb_vars, only: ineb
  use trajenemod, only: trajene

  ! CHARMM support
  use charmm_mod, only : charmm_active, charmm_deallocate_arrays, &
                         charmm_filter_out_qm_atoms
  use ff11_mod, only : cmap_active, deallocate_cmap_arrays
  use memory_module

  ! Self-Guided molecular/Langevin Dynamics (SGLD)
  use sgld, only: isgld, psgld

  use nbips, only: ipssys,ips
  use crg_reloc, only: ifcr, cr_backup_charge, cr_cleanup, cr_allocate, &
                       cr_read_input, cr_check_input, cr_print_info
  use emap,only: temap,pemap,qemap

  use file_io_dat
  use constantph, only: cnstph_finalize
  use constante, only: cnste_finalize
  use external_module, only : external_cleanup
  use barostats, only: mcbar_setup
  use random, only: amrset

  ! Accelerated MD
  use amd_mod

  ! scaledMD
  use scaledMD_mod

  ! Adaptive Buffered Force QM/MM
  use abfqmmm_module

  use music_module, only: read_music_nml, print_music_settings

  use commandline_module, only: cpein_specified

  ! ASM
#ifdef MPI
  use string_method, only : string_define, string_defined
#endif

  implicit none

  logical belly, erstop
  integer ier,ncalls,xmin_iter
  logical ok
  logical newstyle
#include "nmr.h"
#include "box.h"
#include "../include/md.h"
#include "extra.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "parallel.h"

  ! AMBER/MPI
#ifdef MPI
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
  include 'mpif.h'
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
#  ifdef MPI_BUFFER_SIZE
  integer*4 mpibuf(mpi_buffer_size)
#  endif
!  REMD: loop is the current exchange. runmd is called numexchg times.
  integer loop

  _REAL_ emtmd
#endif /* MPI */
  ! End of declarations and inclues for AMBER/MPI

#include "ew_pme_recip.h"
#include "ew_frc.h"
#include "ew_erfc_spline.h"
#ifdef MPI
#  include "ew_parallel.h"
#endif
#include "ew_mpole.h"
#include "ew_cntrl.h"
#include "def_time.h"

  type(state_rec) ::  ene
  integer native,nr3,nr

  ! nmrcal vars
  _REAL_ enmr(6), devdis(4), devang(4), devtor(4), &
         devplpt(4), devpln(4), devgendis(4)
  _REAL_ ag(1), bg(1), cg(1)

  ! Updated 9/2007 by Matthew Seetin to enable
  ! plane-point and plane-plane restraints
  integer numphi, nhb

  ! runmin/trajene var
  _REAL_ carrms

  character(len=8) initial_date, setup_end_date, final_date
  character(len=10) initial_time, setup_end_time, final_time

  ! Needed for final time printout
  integer nstlim_total
  integer, dimension(:), allocatable :: dummy

  _REAL_ time0, time1

  integer i

  integer, dimension(:), allocatable :: ipairs
  logical qsetup
#ifdef MPI
  integer :: n_force_calls
  logical :: do_list_update=.false.
#endif
  _REAL_ :: box_center(3)
#ifdef MPI_DEBUGGER
  integer, volatile :: release_debug

  ! So only the master master thread has release_debug = 0
  release_debug = worldrank

  ! Lock us into an infinite loop while release_debug == 0 on any thread (only
  ! the master here). This allows you to connect a debugger to any running
  ! process without having to 'race' program execution.  A debugger MUST be
  ! attached to the master thread (typically the thread with the lowest PID),
  ! and have release_debug set to NOT 0 (e.g., via "set release_debug=1").
  do
    if (release_debug .ne. 0) then
      exit
    end if
  end do

  ! Prevent any other threads from progressing past this point until all
  ! threads you want to watch with a debugger are watched and all those
  ! threads are continued.
  call mpi_barrier(mpi_comm_world, ier)
#endif

  ! ASM
#ifdef MPI
  inquire(file = "STRING", exist=string_defined)
#endif

  ! Here begin the executable statements.
  call Trace_enter( 'sander' )
  ! Uncomment this to compare master vs nonmaster code paths:
  ! call Trace_set_output_focus( 1 )
  ! Start by initializing the cpu timer.
  ! Needed for machines where returned cpu times are relative.
  call date_and_time( initial_date, initial_time )
  call wallclock( time0 )
  call init_timers()

  ! Initialize the printing of ongoing time and performance summaries.
  call print_ongoing_time_summary(0,0,0.0d0,7)

  ! Initialize the number of copies -- always assume 1 until we know otherwise
  ncopy = 1

  ! BPR - original location of PUPIL interface. I moved it further down
  ! because, if it's here, it can't print stuff; write(6,...) statements
  ! assume mdread1() has already been invoked. However, moving this down
  ! may break other things.

  ! Flag to tell list builder to print size of list on first call
  first_list_flag = .true.

  ! Flag to tell recip space routines to allocate on first call
  first_pme = .true.

  ! Initialize first_call flags for QMMM
  qmmm_struct%qm_mm_first_call = .true.
  qmmm_struct%fock_first_call = .true.
  qmmm_struct%fock2_2atm_first_call = .true.
  qmmm_struct%qm2_allocate_e_repul_first_call = .true.
  qmmm_struct%qm2_calc_rij_eqns_first_call = .true.
  qmmm_struct%qm2_scf_first_call = .true.
  qmmm_struct%zero_link_charges_first_call = .true.
  qmmm_struct%adj_mm_link_pair_crd_first_call = .true.
  qmmm_struct%num_qmmm_calls = 0

#ifdef MPI
  ! Parallel initialization (setup is done in multisander.F90).
  ! Make PE 0 the master.
  master = (mytaskid == 0)
  master_master = (masterrank == 0)
  if (master .and. numtasks > MPI_MAX_PROCESSORS) then
    write(0, '(a,i4,a,i4)') &
         'Error: the number of processors must not be greater than ', &
         MPI_MAX_PROCESSORS, ', but is ', numtasks
    call mexit(6,1)
  end if
#  ifdef MPI_BUFFER_SIZE
  call mpi_buffer_attach(mpibuf, mpi_buffer_size*4, ier)
#  endif
#else /* not MPI follows */
  ! In the single-threaded version, the one process is master
  master = .true.
#endif /* MPI */
  temp0les = 0.d0
  erstop = .false.
  qsetup = .true.

  ! Generic packing scheme
  nwdvar = 1
  native = 32
#ifdef ISTAR2

  ! Int*2 packing scheme
  nwdvar = 2
#endif  /*ISTAR2*/
  numpk = nwdvar
  nbit = native/numpk

  ! Only the master node (only node when single-process)
  ! performs the initial setup and reading/writing
  call timer_start(TIME_TOTAL)
  call abfqmmm_init_param()

  ! In abfQM/MM simulations the QM region changes dynamically, which
  ! requires re-initialization of related data structures.  This is
  ! implemented by running multiple MD simulations in sequence.  For
  ! regular MM or QM/MM MD simulations, we go through this loop only once.
  abfqmmmmasterdoloop: &
  do while ((abfqmmm_param%qmstep <= abfqmmm_param%maxqmstep) .or. &
            (abfqmmm_param%maxqmstep == 0 .and. abfqmmm_param%system == 2))
    masterwork: if (master) then
      if (abfqmmm_param%abfqmmm == 0) then

        ! First, initial reads to determine memory sizes
        call mdread1()
        call amopen(8,parm,'O','F','R')
        call rdparm1(8)
        if (mtmd /= 'mtmd' .or. itgtmd == 2) then
          call mtmdlx(natom)
        end if

        ! Now, we can allocate memory
        call locmem()
        write(6,'(/,a,5x,a)') '|','Memory Use     Allocated'
        write(6,'(a,5x,a,i14)') '|', 'Real      ', lastr
        write(6,'(a,5x,a,i14)') '|', 'Hollerith ', lasth
        write(6,'(a,5x,a,i14)') '|', 'Integer   ', lasti
        write(6,'(a,5x,a,i14)') '|', 'Max Pairs ', lastpr

        ! Dynamic memory allocation: allocate space for
        ! module molecule in the master node
        mol_info%natom = natom
        mol_info%nres  = nres
        call allocate_molecule()

        ! Allocate all global arrays
        allocate( x(lastr), ix(lasti), ipairs(lastpr), ih(lasth), stat=ier)
        REQUIRE(ier == 0)
        ix(1:lasti) = 0

        ! This sets up pointer arrays in MEMORY_MODULE to match array-offsets
        ! into the shared X, IX, and IH arrays. Eventually, LOCMEM code should
        ! be merged with MEMORY_MODULE to allocate individual allocatable
        ! arrays, but that will also require updating the MPI code to handle
        ! individual arrays.
        call memory_init()

        ! Allocate the parm arrays
        call allocate_parms()

        if ((igb .ne. 0 .and. igb .ne. 10 .and. ipb == 0) .or. &
            hybridgb > 0 .or. icnstph > 1 .or. icnste > 1) then
          call allocate_gb( natom, ncopy )
        end if
        if (idecomp > 0) then
#ifdef MPI
          if (ifsc > 0) then
            call synchronize_dec(natom, nres)
          else
            nat = natom
            nrs = nres
          end if
#else
          nat = natom
          nrs = nres
#endif
          call allocate_int_decomp(natom)
        else
          call allocate_int_decomp(1)
        end if

        write(6,'(a,5x,a,i14)'  ) '|', 'nblistReal', nblist_allreal
        write(6,'(a,5x,a,i14)'  ) '|', 'nblist Int', nblist_allint
        write(6,'(a,5x,a,i14,a)') '|', '  Total   ', &
              (8*(lastr+lastrst+nblist_allreal)  &
              + 4*(lasth+lasti+lastpr+lastist+nblist_allint))/1024, &
              ' kbytes'

        ! Finish reading the prmtop file and other user input:
        call rdparm2(x, ix, ih, 8)

        ! ntf and ntc get reset if this is an amoeba prmtop
        call AMOEBA_readparm(8, ntf, ntc, natom, x(lmass))
        if ( ipgm /= 0 ) call POL_GAUSS_readparm(8, natom)
      end if

      call sebomd_setup
      ! Branch for QM/MM systems
      if (qmmm_nml%ifqnt .or. abfqmmm_param%abfqmmm == 1) then
        if (abfqmmm_param%abfqmmm == 0) then
          call read_qmmm_nm_and_alloc(igb, ih, ix, x, cut, use_pme, ntb, 0, &
                                      dummy, 0, .true.)
          if (qmmm_nml%qmtheory%SEBOMD) then

            ! Don't do QM/MM
            qmmm_nml%ifqnt = .false.
            sebomd_obj%do_sebomd = .true.
          end if
        end if
        if (qmmm_struct%abfqmmm == 1 .and. abfqmmm_param%abfqmmm == 0) then
          call abfqmmm_setup(natom, nres, ix(i02), ih(m04), ih(m02), &
                             x(lmass), nbonh, nbona, ix(iibh), ix(ijbh), &
                             ix(iiba), ix(ijba))
          nr = natom
          call AMOEBA_check_newstyle_inpcrd(inpcrd, newstyle)
          if (newstyle) then
            call AM_RUNMD_get_coords(natom, t, irest, ntb, x(lcrd), x(lvel))
          else
#ifdef MPI
            call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, &
                        box, irest, t, temp0, .FALSE., solvph, solve, stagid)
#else
            call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, &
                        box, irest, t, .FALSE.)
#endif
          end if
          abfqmmm_param%maxqmstep = nstlim
        end if
        if (abfqmmm_param%abfqmmm == 1) then
          if (abfqmmm_param%system == 1) then
            call abfqmmm_update_qmatoms(x(lcrd))
            if (abfqmmm_param%ntwpdb < 0) then
              call abfqmmm_write_pdb(x(lcrd), ix(i70))
              close(6)
              call mexit(6, 1)
            end if
          end if
          call abfqmmm_select_system_qmatoms(natom)
          if (qmmm_nml%ifqnt) then
            call read_qmmm_nm_and_alloc(igb, ih, ix, x, cut, use_pme, ntb, &
                                        abfqmmm_param%qmstep, &
                                        abfqmmm_param%isqm, &
                                        abfqmmm_param%abfcharge, .true.)
          endif
        endif
      endif

      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
        call mdread2(x, ix, ih)
        call read_music_nml()
        call print_music_settings()
      endif
#ifdef MPI
      if (ifmbar .ne. 0) then
        call setup_mbar(nstlim)
      end if
#endif

#if defined(RISMSANDER)
      call rism_setparam(mdin, commsander, natom, ntypes, x(L15:L15+natom-1), &
                         x(LMASS:LMASS+natom-1), cn1, cn2, &
                         ix(i04:i04+ntypes**2-1), ix(i06:i06+natom-1))
#endif /*RISMSANDER*/
      if (ifcr .ne. 0) then
        call cr_read_input(natom)
        call cr_check_input(ips)
        call cr_backup_charge( x(l15), natom )
      end if

      ! Allocate memory for energy decomposition
      ! module (needs info from mdread2)
      if (idecomp == 1 .or. idecomp == 2) then
        call allocate_real_decomp(nrs)
#ifdef MPI
        ! Thermodynamic Integration decomposition
        if (ifsc > 0) then
          ! following lines don't really seem to make sense(?)
          ! partner = ieor(masterrank,1)
          ! if (nat == natom) then
          !    nrank = masterrank
          ! else
          !    nrank = partner
          ! end if
          call mpi_bcast(jgroup, nat, MPI_INTEGER, 0, commmaster, ier)
        end if
#endif
      else if( idecomp == 3 .or. idecomp == 4 ) then
        call allocate_real_decomp(npdec*npdec)
      end if

      ! Evaluate constants frommderead settings
      nr = nrp
      nr3 = 3*nr
      belly = (ibelly > 0)

      ! The PUPIL interface was moved down here so that
      ! write() statements work as advertised.
#ifdef PUPIL_SUPPORT
      ! Initialise the CORBA interface
      puperror = 0

      ! We set only one quantum domain by default
      pupnumdomains = 1

      call fixport()
      call inicorbaintfcmd(puperror)
      if (puperror .ne. 0) then
        write(6,*) 'Error creating PUPIL CORBA interface.'
        call mexit(6, 1)
      end if
      pupactive = .true.
      write(6,*) 'PUPIL CORBA interface initialized.'

      ! Allocation of memory and initialization
      pupStep = 0
      puperror = 0
      allocate (qcell   (12         ),stat=puperror)
      allocate (atmqmdomains(2*natom),stat=puperror)
      allocate (pupqlist(natom      ),stat=puperror)
      allocate (pupatm  (natom      ),stat=puperror)
      allocate (pupchg  (natom      ),stat=puperror)
      allocate (qfpup   (natom*3    ),stat=puperror)
      allocate (qcdata  (natom*9    ),stat=puperror)
      allocate (keyMM   (natom      ),stat=puperror)
      allocate (pupres  (nres       ),stat=puperror)
      allocate (keyres  (nres       ),stat=puperror)

      if (puperror .ne. 0) then
        write(6,*) 'Error allocating PUPIL interface memory.'
        call mexit(6, 1)
      end if

      ! Initialise the "atomic numbers" and "quantum forces" vectors
      pupparticles = 0
      iresPup = 1
      pupres(1) = 1
      do iPup = 1, natom
        bs1 = (iPup - 1)*3
        call get_atomic_number_pupil(ih(iPup+m06-1), x(lmass+iPup-1), &
                                     pupatm(iPup))
        if (iresPup .lt. nres) then
          if (iPup .ge. ix(iresPup+i02)) then
            iresPup = iresPup + 1
            pupres(iresPup) = iPup
          end if
        end if
        write (strAux, "(A4,'.',A4)") trim(ih(iresPup+m02-1)), &
                                      adjustl(ih(iPup+m04-1))
        keyres(iresPup) = trim(ih(iresPup+m02-1))
        keyMM(iPup) = trim(strAux)

        ! Retrieve the initial charges
        pupchg(iPup) = x(L15+iPup-1)
        do jPup = 1, 3
          qfpup(bs1+jPup) = 0.0d0
        end do
      end do
      write(6,*) 'Got all atomic numbers.'

      ! Initialise the PUPIL cell
      do iPup = 1, 12
        qcell(iPup) = 0.0d0
      end do

      ! Submit the KeyMM particles and their respective atomic numbers to PUPIL
      puperror = 0
      call putatomtypes(natom, puperror, pupatm, keyMM)
      if (puperror .ne. 0) then
        write(6,*) 'Error sending MM atom types to PUPIL.'
        call mexit(6, 1)
      end if

      ! Submit the Residue Pointer vector to PUPIL
      write(6, "(a20,1x,i6,3x,a17,1x,i6)") 'Number of residues =', nres, &
                                           'Number of atoms =', natom
      puperror = 0
      call putresiduetypes(nres, puperror, pupres, keyres)
      if (puperror .ne. 0) then
        write(6,*) 'Error sending MM residue types to PUPIL.'
        call mexit(6, 1)
      end if
      write(6,*) 'Sent system data to PUPIL.'
      write(*,*) 'PUPIL structure initialized.'
#endif
      ! End of PUPIL interface

      ! Seed the random number generator (master node only)
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
#ifdef MPI
        if (rem == 0) then
          call amrset(ig)
        else

          ! Set the random seed different for different replicas,
          ! but keep same seed for CPUs in the same replica since
          ! we want data from diff numbers of cpus to match.
          !
          ! The variable nodeid is declared through md.h
          ! and is equal to repnum-1
          call amrset(ig + (17 * nodeid))
        end if
#else
        call amrset(ig)
#endif /* MPI */
        if (nbit < 32 .and. nr > 32767) then
          write(6, *) '  Too many atoms for 16 bit pairlist -'
          write(6, *) '    Recompile without ISTAR2'
          call mexit(6, 1)
        end if
        if (ntp > 0.and.iabs(ntb) /= 2) then
          write(6,*) 'Input of NTP/NTB inconsistent'
          call mexit(6, 1)
        end if
      end if

      ! Read coordinates and velocities.  This begins the qmstep = 1 block
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
        call timer_start(TIME_RDCRD)
#ifdef LES
#  ifdef MPI
        call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, t, &
                    temp0les, .TRUE., solvph, solve, stagid)
#  else
        call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, t, &
                    .TRUE.)
#  endif
#else
        call AMOEBA_check_newstyle_inpcrd(inpcrd, newstyle)
        if (newstyle) then
          call AM_RUNMD_get_coords(natom, t, irest, ntb, x(lcrd), x(lvel))
        else
#  ifdef MPI
          call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, &
                      t, temp0, .TRUE., solvph, solve, stagid)
#  else
          call getcor(nr, x(lcrd), x(lvel), x(lforce), ntx, box, irest, t, &
                      .TRUE.)
#  endif
        end if
#endif /* LES */
        if (iamoeba > 0) then
          natom = natom*am_nbead
          nrp = nrp*am_nbead
          nr = nr*am_nbead
          nr3 = nr3*am_nbead
          ncopy = am_nbead
        end if

        ! If this is a polarizable model, read input dipole information
        if (igb == 0 .and. ipb == 0 .and. induced > 0) then
          call get_dips(x, nr)
        end if

#ifdef APBS
        ! APBS initialization
        if (mdin_apbs) then

          ! IN: natom, coords, charge and radii (from prmtop)
          ! OUT: pb charges and pb radii (via apbs_vars module)
          call apbs_init(natom, x(lcrd), x(l15), x(l97))
        end if
#endif /* APBS */

        ! Set the initial velocities
        if (ntx <= 3) then
          call setvel(nr, x(lvel), x(lwinv), tempi, iscale, scalm)

          ! Random numbers may have been "used up" in setting the intial
          ! velocities; re-set the generator so that all nodes are back in
          ! sync.  This again happens only on the master node.
#ifdef MPI
          if (rem == 0) then
            call amrset(ig)
          end if
#else
          call amrset(ig)
#endif
        end if
        if (belly) then
          call bellyf(natom, ix(ibellygp), x(lvel))
        end if
        call timer_stop(TIME_RDCRD)
        if (abfqmmm_param%abfqmmm == 1 .and. ntb > 0) then
          call iwrap2(abfqmmm_param%n_user_qm, abfqmmm_param%user_qm, &
                      x(lcrd), box_center)
        end if

        ! If we are reading NMR restraints/weight changes, read them now:
        if (nmropt >= 1) then
          call nmrcal(x(lcrd), x(lforce), ih(m04), ih(m02), ix(i02), &
                      x(lwinv), enmr, devdis, devang, devtor, devplpt, &
                      devpln, devgendis, temp0, tautp, cut, x(lnmr01), &
                      ix(inmr02), x(l95), 5, 6, rk, tk, pk, cn1, cn2, ag, &
                      bg, cg, numbnd, numang, numphi, nimprp, nhb, natom, &
                      natom, ntypes, nres, rad, wel, radhb, welhb, rwell, &
                      tgtrmsd, temp0les, -1, 'READ')

          ! Updated 9/2007 by Matthew Seetin to enable plane-point and
          ! plane-plane restraints.  Determine how many of the torsional
          ! parameters are impropers
          call impnum(ix(i46), ix(i56), ix(i48), ix(i58), nphih, nphia, &
                      0, nptra, nimprp)
        end if

        ! Set up info related to weight changes for the non-bonds:
        call nmrrad(rad, wel, cn1, cn2, ntypes, 0, 0.0d0)
        call decnvh(asol, bsol, nphb, radhb, welhb)
        if (iredir(4) > 0) then
          call noeread(x, ix, ih)
        end if
        if (iredir(8) > 0) then
          call alignread(natom, x(lcrd))
        end if
        if (iredir(9) > 0) then
          call csaread
        end if
      end if
      ! End of the qmstep = 1 block

      ! Call the fastwat subroutine to tag those bonds which are part
      ! of 3-point water molecules. Constraints will be performed for
      ! these waters using a fast analytic routine.
      call timer_start(TIME_FASTWT)
      call fastwat(ih(m04), nres, ix(i02), ih(m02), nbonh, nbona, ix(iibh), &
                   ix(ijbh), ibelly, ix(ibellygp), iwtnm, iowtnm, ihwtnm, &
                   jfastw, ix(iifstwt), ix(iifstwr), ibgwat, ienwat, ibgion, &
                   ienion, iorwat, 6, natom)
      call timer_stop(TIME_FASTWT)
      call getwds(ih(m04), nres, ix(i02), ih(m02), nbonh, nbona, 0, ix(iibh), &
                  ix(ijbh), iwtnm, iowtnm, ihwtnm, jfastw, ix(iicbh), req, &
                  x(lwinv), rbtarg, ibelly, ix(ibellygp), 6)

      ! Assign link atoms between quantum mechanical and molecular mechanical
      ! atoms if quantum atoms are present.  After assigning the link atoms,
      ! delete all connectivity between the QM atoms.
      if (qmmm_nml%ifqnt) then
        if (abfqmmm_param%abfqmmm .ne. 1) then
          call identify_link_atoms(nbona, ix(iiba), ix(ijba))
        end if

        ! Variable QM solvent:
        ! Store the original bond parameters since we will need to rebuild
        ! the QM region (delete bonded terms etc) repeatedly
        if (qmmm_nml%vsolv > 0) then
          call new(qmmm_vsolv, nbonh, nbona, ntheth, ntheta, nphih, nphia)
          call qmmm_vsolv_store_parameters(qmmm_vsolv, numbnd, ix(iibh), &
                                           ix(ijbh), ix(iicbh), ix(iiba), &
                                           ix(ijba), ix(iicba), ix(i24), &
                                           ix(i26), ix(i28), ix(i30), &
                                           ix(i32), ix(i34), ix(i36), &
                                           ix(i38), ix(i40), ix(i42), &
                                           ix(i44), ix(i46), ix(i48), &
                                           ix(i50), ix(i52), ix(i54), &
                                           ix(i56), ix(i58))
        end if
        if (abfqmmm_param%abfqmmm == 1) then
          if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
            call abfqmmm_allocate_arrays_of_parameters(numbnd, nbonh, nbona, &
                                                       ntheth, ntheta, nphih, &
                                                       nphia)
            call abfqmmm_store_parameters(ix(iibh), ix(ijbh), ix(iicbh), &
                                          ix(iiba), ix(ijba), ix(iicba), &
                                          ix(i24), ix(i26), ix(i28), ix(i30), &
                                          ix(i32), ix(i34), ix(i36), ix(i38), &
                                          ix(i40), ix(i42), ix(i44), ix(i46), &
                                          ix(i48), ix(i50), ix(i52), ix(i54), &
                                          ix(i56), ix(i58), x(l15), rk, req)
          else
            call abfqmmm_set_parameters(numbnd, nbonh, nbona, ntheth, ntheta, &
                                        nphih, nphia, ix(iibh), ix(ijbh), &
                                        ix(iicbh), ix(iiba), ix(ijba), &
                                        ix(iicba), ix(i24), ix(i26), ix(i28), &
                                        ix(i30), ix(i32), ix(i34), ix(i36), &
                                        ix(i38), ix(i40), ix(i42), ix(i44), &
                                        ix(i46), ix(i48), ix(i50), ix(i52), &
                                        ix(i54), ix(i56), ix(i58), x(l15), &
                                        rk, req)
            call init_extra_pts(ix(iibh), ix(ijbh), ix(iicbh), ix(iiba), &
                                ix(ijba), ix(iicba), ix(i24), ix(i26), &
                                ix(i28), ix(i30), ix(i32), ix(i34), ix(i36), &
                                ix(i38), ix(i40), ix(i42), ix(i44), ix(i46), &
                                ix(i48), ix(i50), ix(i52), ix(i54), ix(i56), &
                                ix(i58), ih(m06), ix, x, ix(i08), ix(i10), &
                                nspm, ix(i70), x(l75), tmass, tmassinv, &
                                x(lmass), x(lwinv), req)
          end if
          call identify_link_atoms(nbona, ix(iiba), ix(ijba))
        end if

        ! Remove bonds between QM atoms from list (Hydrogen)
        if (nbonh .gt. 0) then
          call setbon(nbonh, ix(iibh), ix(ijbh), ix(iicbh), ix(ibellygp))
        end if

        ! Remove bonds between QM atoms from list (Heavy)
        if (nbona .gt. 0) then
          call setbon(nbona, ix(iiba), ix(ijba), ix(iicba), ix(ibellygp))
        end if

        ! Remove angles between QM atoms from list (Hydrogen)
        if (ntheth .gt. 0) then
          call setang(ntheth, ix(i24), ix(i26), ix(i28), ix(i30), ix(ibellygp))
        end if

        ! Remove angles between QM atoms from list (Heavy)
        if (ntheta .gt. 0) then
          call setang(ntheta, ix(i32), ix(i34), ix(i36), ix(i38), ix(ibellygp))
        end if

        ! Remove dihedrals between QM atoms from list (Hydrogen)
        if (nphih .gt. 0) then
          call setdih(nphih, ix(i40), ix(i42), ix(i44), ix(i46), &
                      ix(i48), ix(ibellygp))
        end if

        ! Remove dihedrals between QM atoms from list (Heavy)
        if (nphia .gt. 0) then
          call setdih(nphia, ix(i50), ix(i52), ix(i54), ix(i56), &
                      ix(i58), ix(ibellygp))
        end if

        ! Remove CHARMM energy terms from QM region
        call charmm_filter_out_qm_atoms()

        ! Now we should work out the type of each quantum atom present.
        ! This is used for our arrays of pre-computed parameters. It is
        ! essentially a re-basing of the atomic numbers and is done to save
        ! memory. Note: qm_assign_atom_types will allocate the qm_atom_type
        ! array for us. Only the master calls this routine. All other
        ! threads get this allocated and broadcast to them by the mpi setup
        ! routine.
        call qm_assign_atom_types

        ! Set default QMMM MPI parameters - for single cpu operation.
        ! These will get overwritten by qmmm_mpi_setup if MPI is on.
        qmmm_mpi%commqmmm_master = master
        qmmm_mpi%numthreads = 1
        qmmm_mpi%mytaskid = 0
        qmmm_mpi%natom_start = 1
        qmmm_mpi%natom_end = natom
        qmmm_mpi%nquant_nlink_start = 1
        qmmm_mpi%nquant_nlink_end = qmmm_struct%nquant_nlink

        ! Now we know how many link atoms we can allocate the scf_mchg array...
        allocate(qm2_struct%scf_mchg(qmmm_struct%nquant_nlink), stat=ier)

        ! qm2_struct%scf_mchg will be deallocated in deallocate qmmm,
        ! so the allocation must have been successful here.
        REQUIRE(ier == 0)

        ! We can also allocate ewald_memory
        if (qmmm_nml%qm_ewald > 0) then
          call allocate_qmewald(natom)
        end if
        if (qmmm_nml%qmgb == 2) then
          call allocate_qmgb(qmmm_struct%nquant_nlink)
        end if
        allocate(qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat=ier)

        ! Again, qmmm_struct%dxyzqm will be deallocated in deallocate qmmm
        REQUIRE(ier == 0)
      else if(abfqmmm_param%abfqmmm == 1) then
        call abfqmmm_set_parameters(numbnd, nbonh, nbona, ntheth, ntheta, &
                                    nphih, nphia, ix(iibh), ix(ijbh), &
                                    ix(iicbh), ix(iiba), ix(ijba), ix(iicba), &
                                    ix(i24), ix(i26), ix(i28), ix(i30), &
                                    ix(i32), ix(i34), ix(i36), ix(i38), &
                                    ix(i40), ix(i42), ix(i44), ix(i46), &
                                    ix(i48), ix(i50), ix(i52), ix(i54), &
                                    ix(i56), ix(i58), x(l15), rk, req)
        call init_extra_pts(ix(iibh), ix(ijbh), ix(iicbh), ix(iiba), &
                            ix(ijba), ix(iicba), ix(i24), ix(i26), ix(i28), &
                            ix(i30), ix(i32), ix(i34), ix(i36), ix(i38), &
                            ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), &
                            ix(i50), ix(i52), ix(i54), ix(i56), ix(i58), &
                            ih(m06), ix, x, ix(i08), ix(i10), nspm, ix(i70), &
                            x(l75), tmass, tmassinv, x(lmass), x(lwinv), req)
      end if
      ! End if branch based on the presence of quantum atoms.  This block
      ! assigned link atoms and replaced MM interactions along the QM/MM
      ! boundary if necessary.

      ! Open the data dumping files and position them
      ! depending on the type of run:
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
#ifdef MPI
        ! Adaptive QM/MM (qmmm_nml%vsolv=2) via multisander.
        ! All groups have identical coords and velocities.
        ! Only the master process needs to dump results.
        if ((qmmm_nml%vsolv < 2) .or. (worldrank == 0)) then
#endif
          call open_dump_files
#ifdef MPI
        end if
#endif
      end if
      if (master) then
        call amflsh(6)
      end if

    end if masterwork
    ! End of master process setup

#if defined(RISMSANDER)
    call rism_init(commsander)
#endif /* RISMSANDER */

#ifdef MPI
    call mpi_barrier(commsander,ier)

    ! AMBER/MPI
    !
    ! NOTE: in the current AMBER/MPI implementation, two means of
    ! running in parallel within sander are supported. The value
    ! of mpi_orig determines which approach is used.
    ! This is turned on when minimization (imin .ne. 0) is requested,
    ! and is otherwise off.  When it is on the master and the other
    ! processes follow substantially different code paths, and in that
    ! case collective MPI operations must not be used in those code paths.
    !
    ! When running the mpi_orig case, a variable notdone is now
    ! set by the master and determines when to exit the force()
    ! loop.  When the master has finished calling force, the
    ! master changes notdone to 0 and broadcasts the data one more
    ! time to signal end of the loop.  force() is modified so that
    ! in the mpi_orig case, an initial broadcast is done to receive
    ! the value from the master to decide whether to do the work or
    ! simply exit.
    !
    ! Set up initial data and send all needed data to other nodes,
    ! now that the master has it
    !
    ! First, broadcast parameters in memory.h, so that all processors
    ! will know how much memory to allocate:
    call mpi_bcast(natom, BC_MEMORY, mpi_integer, 0, commsander, ier)

    ! Thermodynamic Integration decomposition
    call mpi_bcast(idecomp, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(nat, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(nrs, 1, mpi_integer, 0, commsander, ier)

    ! Set up integer stack initial size
    call mpi_bcast(lastist, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(lastrst, 1, mpi_integer, 0, commsander, ier)

    call mpi_bcast(qmmm_struct%abfqmmm, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(abfqmmm_param%abfqmmm, 1, mpi_integer, 0, commsander, ier)
    if (abfqmmm_param%abfqmmm == 1) then
      call mpi_bcast(abfqmmm_param%natom, 1, mpi_integer, 0, commsander, ier)
      call mpi_bcast(abfqmmm_param%qmstep, 1, mpi_integer, 0, commsander, ier)
      call mpi_bcast(abfqmmm_param%maxqmstep, 1, mpi_integer, 0, &
                     commsander, ier)
      call mpi_bcast(abfqmmm_param%system, 1, mpi_integer, 0, commsander,ier)
    end if

    if (abfqmmm_param%abfqmmm == 1) then
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
        call mpi_bcast(abfqmmm_param%r_core_in, 1, mpi_double_precision, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%r_core_out, 1, mpi_double_precision, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%r_qm_in, 1, mpi_double_precision, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%r_qm_out, 1, mpi_double_precision, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%r_buffer_in, 1, mpi_double_precision, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%r_buffer_out, 1, mpi_double_precision, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%mom_cons_type, 1, mpi_integer, &
                       0, commsander, ier)
        call mpi_bcast(abfqmmm_param%mom_cons_region, 1, mpi_integer, &
                       0, commsander, ier)
        if (.not. master) then
          allocate(abfqmmm_param%x(3*natom), stat=ier)
          REQUIRE(ier == 0)
          allocate(abfqmmm_param%id(natom), stat=ier)
          REQUIRE(ier == 0)
          allocate(abfqmmm_param%mass(natom), stat=ier)
          REQUIRE(ier == 0)
        end if
        allocate(abfqmmm_param%v(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
        allocate(abfqmmm_param%f(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
        allocate(abfqmmm_param%f1(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
        allocate(abfqmmm_param%f2(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
      end if
      call mpi_bcast(abfqmmm_param%x, 3*natom, mpi_double_precision, &
                     0, commsander, ier)
      call mpi_bcast(abfqmmm_param%id, natom, mpi_integer, &
                     0, commsander, ier)
      call mpi_bcast(abfqmmm_param%mass, natom, mpi_double_precision, &
                     0, commsander, ier)
    end if
    if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      call stack_setup()
    else
      call deallocate_stacks
      call stack_setup()
    end if
    call mpi_bcast(plumed, 1, MPI_INTEGER, 0, commsander, ier)
    if (plumed .eq. 1) then
      call mpi_bcast(plumedfile, MAX_FN_LEN, MPI_CHARACTER, 0, commsander, ier)
    endif

    ! GMS: Broadcast parameters from module 'molecule'
    call mpi_bcast(mol_info%natom, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(mol_info%nres, 1, mpi_integer, 0, commsander, ier)
    call mpi_barrier(commsander, ier)

    ! Allocate memory on the non-master nodes
    notmasterallocation: if (.not. master) then

      if (abfqmmm_param%qmstep .ne. 1 .or. abfqmmm_param%system .ne. 1) then
        deallocate(x, ix, ipairs, ih)
      end if

      allocate(x(1:lastr), stat=ier)
      REQUIRE(ier == 0)

      allocate(ix(1:lasti), stat=ier)
      REQUIRE(ier == 0)

      allocate(ipairs(1:lastpr), stat=ier)
      REQUIRE(ier == 0)

      allocate(ih(1:lasth), stat=ier)
      REQUIRE(ier == 0)

      ! AWG Set up pointer arrays also on non-master nodes
      call memory_init()

      ! Allocate space for molecule module arrays in the other nodes
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
        call allocate_molecule()
      else
        call deallocate_molecule()
        call allocate_molecule()
      end if

      ! Thermodynamic Integration decomposition
      if (idecomp > 0) then
        call allocate_int_decomp(natom)
        if (idecomp == 1 .or. idecomp == 2) then
          call allocate_real_decomp(nrs)
        else if( idecomp == 3 .or. idecomp == 4 ) then
          call allocate_real_decomp(npdec*npdec)
        end if
      end if
    end if notmasterallocation
    ! End memory allocation on non-master nodes

    call Trace_note( 'After end if notmasterallocation' )

    ! Broadcast arrays from module 'molecule'
    call mpi_bcast(mol_info%natom_res, mol_info%nres, MPI_INTEGER, &
                   0, commsander, ier)
    call mpi_bcast(mol_info%atom_to_resid_map, mol_info%natom, MPI_INTEGER, &
                   0, commsander, ier)
    call mpi_bcast(mol_info%atom_mass, mol_info%natom, MPI_DOUBLE_PRECISION, &
                   0, commsander, ier)
    if (idecomp == 1 .or. idecomp == 2) then
      call mpi_bcast(jgroup, nat, MPI_INTEGER, 0, commsander, ier)
    end if
    if (icfe == 0 .and. (idecomp ==3 .or. idecomp == 4)) then
      call mpi_bcast(jgroup, natom, MPI_INTEGER, 0, commsander, ier)
      call mpi_bcast(indx, nres, MPI_INTEGER, 0, commsander, ier)
    end if
    call startup_groups(ier)
    call startup(x, ix, ih)

    ! Broadcast Empirical Valence Bond / Path Integral MD inputs/parameters to
    ! all PEs.  Note: the masters have all required EVB/PIMD inputs/parameters
    ! via call to mdread2 (evb_input, evb_init, evb_pimd_init).  For EVB/PIMD,
    ! all PEs need the inputs/parameters ... so we perform this initialization
    ! again for all PEs besides the masters.  The alternative is to use
    ! MPI_BCAST.
    call mpi_bcast(ievb , 1, MPI_INTEGER, 0, commworld, ier)
    call mpi_bcast(ipimd, 1, MPI_INTEGER, 0, commworld, ier)
    if (ievb .ne. 0) then
      call mpi_bcast(evbin, MAX_FN_LEN, MPI_CHARACTER, 0, commworld, ier)
      if (.not. master) then
        call evb_input
        call evb_init
      end if
    end if

#  if defined(LES)
    call mpi_bcast (ncopy, 1, MPI_INTEGER, 0, commworld, ier)
    call mpi_bcast (cnum(1:natom), natom, MPI_INTEGER, 0, commworld, ier)
    call mpi_bcast (evbin, MAX_FN_LEN, MPI_CHARACTER, 0, commworld, ier)
#  endif /* LES */
    if (ievb .ne. 0) then
#  if defined(LES)
      if (ipimd > 0 .and. .not. master) then
        call evb_input
        call evb_init
      end if
      call evb_pimd_init
#  endif /* LES */
      call evb_bcast
      call evb_alloc
    end if

    ! Obtain B vector for Schlegel's distributed Gaussian method
    if (trim(adjustl(xch_type)) == "dist_gauss") then
      call schlegel_dg
    end if
    if (ifsc .ne. 0) then

      ! Multi-CPU minimization does not work with soft core
      if (imin > 0 .and. numtasks > 1) then
        call sander_bomb('imin > 0 and numtasks > 1', &
                     'TI minimizations cannot be performed with > 2 CPUs','')
      end if
      call setup_sc(natom, nres, ih(m04), ih(m06), &
                    ix(i02), ih(m02), x(lcrd), ntypes, clambda, nstlim)
      if (ntp > 0 .and. master) then

        ! Check which molecules are perturbed in NPT runs
        call sc_check_perturbed_molecules(nspm, ix(i70))
      end if

      ! Thermodynamic Integration decomposition
      if (idecomp > 0) then
        if (sanderrank == 0) then
          call build_dec_mask
        end if
        call mpi_bcast(decmask, natom, MPI_INTEGER, 0, commsander, ier)
      end if

      ! Make sure all common atoms have the same v (that of V0) in TI runs
      if (ifsc .ne. 2) then
        if (master) then
          call sc_sync_x(x(lvel), nr3)
        end if
        if (numtasks > 1) then
          call mpi_bcast(nr3, 1, MPI_INTEGER, 0, commsander, ier)
          call mpi_bcast(x(lvel), nr3, MPI_DOUBLE_PRECISION, &
                         0, commsander, ier)
        end if
      end if
      if (tishake .ne. 0) then
        call setnoshake_sc(ix, ntc, num_noshake, master)
      end if
    else
      extra_atoms=0
    end if
    if (.not. master .and. igb == 0 .and. ipb == 0) then
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
        call nblist_allocate(natom,ntypes,num_direct,numtasks)
      else
        call nblist_deallocate
        call nblist_allocate(natom, ntypes, num_direct, numtasks)
      end if
    end if

    ! Allocate memory for GB on the non-master nodes:
    if (.not. master) then
      if ((igb /= 0 .and. igb /= 10 .and. ipb == 0) .or. &
          hybridgb > 0 .or. icnstph.gt.1 .or. icnste.gt.1) then
        if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
          call allocate_gb( natom, ncopy )
        else
          call deallocate_gb
          call allocate_gb(natom, ncopy)
        end if
      end if
    end if
    ! End non-msater process GB memory allocation

    nr = nrp
    nr3 = 3 * nr
    belly = (ibelly > 0)

    ! Do setup for QMMM in parallel if ifqnt is on - this is essentially
    ! everything in the QMMM namelist and a few other bits and pieces.
    ! Note, currently only master node knows if qmmm_nml%ifqnt is
    ! on so we need to broadcast this first and then make decisions
    ! based on this.
    call mpi_bcast(qmmm_nml%ifqnt, 1, mpi_logical, 0, commsander, ier)
    if (qmmm_nml%ifqnt) then

      ! Broadcast all of the stuff in qmmm_nml and allocate the relevant
      ! arrays on all processors. This will also set up information for
      ! openmp on the master processor if it is in use.
      call qmmm_mpi_setup(master, natom)
      if (qmmm_nml%qm_ewald > 0 .and. .not. master) then
        call allocate_qmewald(natom)
      end if
      if (qmmm_nml%qmgb==2 .and. .not. master) then
        call allocate_qmgb(qmmm_struct%nquant_nlink)
      end if
      if (.not. master) then
        allocate(qmmm_struct%dxyzqm(3, qmmm_struct%nquant_nlink), stat=ier)

        ! qmmm_struct%dxyzqm will be deallocated in deallocate qmmm
        REQUIRE(ier == 0)
      end if
    end if
    ! End of QM/MM MPI setup

    ! All nodes are calling this. amrset(ig) has already been called
    ! by the master node (twice if initial velocities are set, ntx <= 3).
    ! Replica Exchange MD now requires need to call amrset on all child
    ! threads, masters have called it above before initial coord read.
    if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      if (rem == 0) then
        call amrset(ig+1)
      else
        if (.not. master) then
          call amrset(ig + 17*nodeid)
        end if
      endif
      if (nmropt >= 1) then

        ! Updated 9/2007 by Matthew Seetin to enable
        ! plane-point and plane-plane restraints
        call nmrcal(x(lcrd), x(lforce), ih(m04), ih(m02), ix(i02), x(lwinv), &
                    enmr, devdis, devang, devtor, devplpt, devpln, devgendis, &
                    temp0, tautp, cut, x(lnmr01), ix(inmr02), x(l95), 5, 6, &
                    rk, tk, pk, cn1, cn2, ag, bg, cg, numbnd, numang, numphi, &
                    nimprp, nhb, natom, natom, ntypes, nres, rad, wel, radhb, &
                    welhb, rwell, tgtrmsd, temp0les, -1, 'MPI ')
      end if
    end if
    call mpi_bcast(lmtmd01, 1, mpi_integer, 0, commsander, ier)
    call mpi_bcast(imtmd02, 1, mpi_integer, 0, commsander, ier)
    if (itgtmd == 2) then
      call mtmdcall(emtmd, x(lmtmd01), ix(imtmd02), x(lcrd), x(lforce), &
                    ih(m04), ih(m02), ix(i02), ih(m06), x(lmass), natom, &
                    nres, 'MPI ')
    end if

    ! Check that the system is neutral and print warning message
    ! if not.  Adjust charges for roundoff error.
    if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then
      if (icfe == 0) then
        call check_neutral(x(l15),natom)
      else
        call ti_check_neutral(x(l15),natom)
      end if
    end if

    ! Tell all nodes if this is a Semi-Empirical Born-Oppenheimer run:
    ! transfer SEBOMD information to the nodes, open necessary files,
    ! and initialize SEBOMD arrays.
    call mpi_bcast(sebomd_obj%do_sebomd, 1, MPI_LOGICAL, 0, commsander, ier)
    if (sebomd_obj%do_sebomd) then
      call sebomd_bcast_obj()
      if (master) then
        call sebomd_open_files
      end if
      call init_sebomd_arrays(natom)
    end if

    ! Use old parallelism for energy minimization.
    ! This is turned on when minimization (imin .ne. 0) is requested,
    ! and is otherwise off.  When it is on the master and the other
    ! processes follow substantially different code paths, and in that
    ! case collective MPI operations must not be used
    if (imin .ne. 0) then
      mpi_orig = .true.
      notdone = 1
    else
      mpi_orig = .false.
    end if
    call Trace_logical( 'mpi_orig is ', mpi_orig )

    notmastermpi_orig: if (mpi_orig .and. .not. master) then

      ! All nodes only do the force calculations.  Minimization
      ! so only master gets past the loop below
      ! hence need to zero QM charges on non-master threads here.
      if (qmmm_nml%ifqnt) then

        ! Apply charge correction if required.
        if (qmmm_nml%adjust_q > 0) then
          call qmmm_adjust_q(qmmm_nml%adjust_q, natom, qmmm_struct%nquant, &
                             qmmm_struct%nquant_nlink, qmmm_struct%nlink, &
                             x(L15), qmmm_struct%iqmatoms, qmmm_nml%qmcharge,&
                             qmmm_struct%atom_mask,qmmm_struct%mm_link_mask, &
                             master,x(LCRD), qmmm_nml%vsolv)
        end if

        ! At this point we can also fill the qmmm_struct%scaled_mm_charges
        ! array - we only need to do this once as the charges are constant
        ! during a run. Having a separate array of scaled charges saves us
        ! having to do it on every qmmm routine call. Do this BEFORE zeroing
        ! the QM charges since that routine take care of these values as well.
        do i = 1, natom

          ! Charge scaling factor for Free Energy Perturbation
          qmmm_struct%scaled_mm_charges(i) = x(L15+(i-1)) * &
                                             INV_AMBER_ELECTROSTATIC * &
                                             qmmm_nml%chg_lambda
        end do

        ! Zero out the charges on the quantum mechanical atoms
        call qm_zero_charges(x(L15), qmmm_struct%scaled_mm_charges, .true.)
        if (qmmm_struct%nlink > 0) then

          ! We need to exclude all electrostatic interactions with
          ! MM link pairs, both QM-MM and MM-MM. Do this by zeroing
          ! the MM link pair charges in the main charge array.  These
          ! charges are stored in qmmm_struct%mm_link_pair_resp_charges
          ! in case they are later needed.
          call qm_zero_mm_link_pair_main_chg(qmmm_struct%nlink, &
                                             qmmm_struct%link_pairs, x(L15), &
                                             qmmm_struct%scaled_mm_charges, &
                                             .true.)
        end if
      end if
      if (igb == 7 .or. igb == 8 .or. hybridgb == 7 .or. hybridgb == 8) then

        ! x(l97) is rborn(), the Generalized Born radii
        call igb7_init(natom, x(l97))

        ! TODO: add igb == 8 here
      end if
      if (ifcr .ne. 0) then
        call cr_allocate(master, natom)
      end if
      n_force_calls = 0
      do while(notdone == 1)
        n_force_calls = n_force_calls+1
        call force(x, ix, ih, ipairs, x(lcrd), x(lforce), ene,ene%vir, &
                   x(l96), x(l97), x(l98), x(l99), qsetup, &
                   do_list_update, n_force_calls)
      end do

      ! Deallocate and return
      goto 999
    end if notmastermpi_orig

    call Trace_note( 'After end if notmastermpi_orig' )

    ! Report the parallel configuration
    if (master) then
      if (abfqmmm_param%abfqmmm .ne. 1) then
        write(6, '(a,i4,a,/)') '|  Running AMBER/MPI version on ', &
                               numtasks, ' nodes'
      end if

      ! BTREE is selected by default if cpu is a power of two.
      ! The number of processes is required to be a power of two for Btree
      ! Print a warning about inefficiency with cpus not being a power of 2.
      if (numtasks > 1 .and. logtwo(numtasks) <= 0) then
        if (abfqmmm_param%abfqmmm .ne. 1) then
          write(6, '(a,i4,a,/)') &
                '|  WARNING: The number of processors is not a power of 2'
        end if
        if (abfqmmm_param%abfqmmm /=1) then
          write(6, '(a,i4,a,/)') &
                '|           this may be inefficient on some systems.'
        end if
      end if
    end if
    ! End reporting of parallel configuration

    if (master .and. numgroup > 1) then
      write(6, '(a,i4,a,i4,a,i4,a)') '|  MULTISANDER: ', numgroup, &
           ' groups. ',  numtasks, ' processors out of ', worldsize, &
           ' total.'
    end if

    if (master) then
      call amflsh(6)
    end if
    ! End of AMBER/MPI work in this block

#else
    if (abfqmmm_param%abfqmmm == 1) then
      if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
        allocate(abfqmmm_param%v(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
        allocate(abfqmmm_param%f(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
        allocate(abfqmmm_param%f1(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
        allocate(abfqmmm_param%f2(3*natom+iscale), stat=ier)
        REQUIRE(ier == 0)
      end if
    end if

    ! For debugging, the charges must be copied at the start so that
    ! they can't change later.  Check that the system is neutral and
    ! print warning message adjust charges for roundoff error.
    if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then
      call check_neutral(x(l15),natom)
    end if
    if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      call amrset(ig+1)
    end if

    if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1) then
      call stack_setup()
    else
      call deallocate_stacks
      call stack_setup()
    end if

    ! Open necessary files and initialize SEBOMD arrays
    if (sebomd_obj%do_sebomd) then
      call sebomd_open_files
      call init_sebomd_arrays(natom)
    endif
#endif /* MPI */

#ifdef OPENMP

    ! If -openmp was specified to configure_amber then -DOPENMP is defined and
    ! the threaded version of MKL will have been linked in. It is important
    ! here that we set the default number of openmp threads for MKL to be 1 to
    ! stop conflicts with threaded vectorization routines when running in
    ! parallel etc.  Individual calls to MKL from routines that know what they
    ! are doing - e.g. QMMM calls to diagonalizers etc can increase this limit
    ! as long as they put it back afterwards.
    call omp_set_num_threads(1)

    ! If we are using openmp for matrix diagonalization print some information.
    if (qmmm_nml%ifqnt .and. master) then
      call qm_print_omp_info()
    end if
#endif

    ! Allocate memory for crg relocation
    if (ifcr /= 0) then
      call cr_allocate( master, natom )
    end if

    ! Initialize LIE module if used
    if ( ilrt /= 0 ) then
      call setup_linear_response(natom, nres, ih(m04), ih(m06), ix(i02), &
                                 ih(m02), x(lcrd), x(l15), ntypes, ix(i04), &
                                 ix(i06), cn1, cn2, master)
    end if
    call date_and_time(setup_end_date, setup_end_time)

    ! Initialize the printing of ongoing time and performance summaries.
    ! We do this quite late here to avoid including all the setup time.
    if (master) then
      call print_ongoing_time_summary(0, 0, 0.0d0, 7)
    end if

    ! Now do the dynamics or minimization.
    if (igb == 7 .or. igb == 8 .or. hybridgb == 7 .or. hybridgb == 8) then
      call igb7_init(natom, x(l97)) !x(l97) is rborn()

      ! Add igb = 8 here
    end if

    if (qmmm_nml%ifqnt) then

      ! Apply charge correction if required.
      if (qmmm_nml%adjust_q > 0) then
        call qmmm_adjust_q(qmmm_nml%adjust_q, natom, qmmm_struct%nquant, &
                           qmmm_struct%nquant_nlink, qmmm_struct%nlink, &
                           x(L15), qmmm_struct%iqmatoms, qmmm_nml%qmcharge, &
                           qmmm_struct%atom_mask, qmmm_struct%mm_link_mask, &
                           master, x(LCRD), qmmm_nml%vsolv)
      end if

      ! At this point we can also fill the qmmm_struct%scaled_mm_charges
      ! array - we only need to do this once as the charges are constant
      ! during a run. Having a separate array of scaled charges saves us
      ! having to do it on every qmmm routine call. Do this BEFORE zeroing
      ! the QM charges since that routine take care of these values as well.
      do i = 1, natom

        ! Charge scaling factors for Free Energy Perturbation
        qmmm_struct%scaled_mm_charges(i) = x(L15+(i-1)) * &
                                           INV_AMBER_ELECTROSTATIC * &
                                           qmmm_nml%chg_lambda
      end do

      ! Zeroing of QM charges MUST be done AFTER call to check_neutral.
      ! Zero out the charges on the quantum mechanical atoms.
      call qm_zero_charges(x(L15), qmmm_struct%scaled_mm_charges, .true.)
      if (qmmm_struct%nlink > 0) then

        ! We need to exclude all electrostatic interactions with
        ! MM link pairs, both in QM-MM and MM-MM.  Do this by zeroing
        ! the MM link pair charges in the main charge array.  These
        ! charges are stored in qmmm_struct%mm_link_pair_resp_charges
        ! in case they are later needed.
        call qm_zero_mm_link_pair_main_chg(qmmm_struct%nlink, &
                                           qmmm_struct%link_pairs, x(L15), &
                                           qmmm_struct%scaled_mm_charges, &
                                           .true.)
      end if
    end if

    ! Use the debugf namelist to activate
    call debug_frc(x, ix, ih, ipairs, x(lcrd), x(lforce), cn1, cn2, qsetup)

    ! Prepare for SGLD simulation
    if (isgld > 0) then
      call psgld(natom,ix(i08), ix(i10), x(lmass),x(lcrd),x(lvel), rem)
    end if

    ! Prepare for EMAP constraints
    if (temap) then
      call pemap(dt,temp0,x,ix,ih)
    end if

    ! Prepare for Isotropic periodic sum of nonbonded interaction
    if (ips .gt. 0) then
      call ipssys(natom, ntypes, ntb, x(l15), cut, cn1, cn2, ix(i04), &
                  ix(i06), x(lcrd))
    end if

    if (master .and. (.not. qmmm_nml%ifqnt) .and. &
        (abfqmmm_param%abfqmmm .ne. 1)) then
      write(6,'(/80(''-'')/,''   4.  RESULTS'',/80(''-'')/)')
    end if

    ! Set up the MC barostat if requested
    if (ntp > 0 .and. barostat == 2) then
      call mcbar_setup(ig)
    end if

    ! ASM
#ifdef MPI
    if (string_defined) call string_define(x(lcrd:lcrd+natom*3-1))
#endif

    ! Input flag imin determines the type of calculation: MD, minimization, ...
    call Trace_integer( 'At label imincase; imin is ', imin )
    imincase: select case (imin)
      case (0)

        ! Dynamics:
        call timer_start(TIME_RUNMD)

        ! Set up Accelerated Molecular Dynamics
        if (iamd .gt. 0) then
          call amd_setup(ntwx)
        endif

        ! Set up scaledMD
        if (scaledMD .gt. 0) then
          call scaledMD_setup(ntwx)
        endif

        ! Path Integral Molecular Dynamics
        if (ipimd > 0) then
          call pimd_init(natom, x(lmass), x(lwinv), x(lvel), ipimd)
        end if

        ! Nudged Elastic Band simulations
        if (ineb > 0) then
          call neb_init()
        end if

#ifdef MPI
        ! Replica Exchange Molecular Dynamics.  If this is not a REMD run,
        ! runmd is called only once.  If this is a REMD run, runmd is
        ! called 0 to numexchg times, where the 0th runmd is just for getting
        ! initial PEs (no dynamics).
        if (rem == 0) then

          ! Not a REMD run. runmd will be called once.
          loop = 0
        else

          ! This is a REMD run. runmd will be called numexchg times.
          loop = numexchg
          if (rem < 0) then

            ! Multi-D REMD run
            call multid_remd_setup(numexchg, numwatkeep, temp0, mxvar, &
                                   natom, ig, solvph, solve, irest)
          else

            ! 1D REMD. Set up temptable, open remlog, etc.
            call remd1d_setup(numexchg, hybridgb, numwatkeep, &
                        temp0, mxvar, natom, ig, solvph, solve, stagid)
          end if
          call setup_pv_correction(6, ntp, sanderrank)

          ! Now set up REMD indices for traj/restart writes. Only do this on
          ! master since only master handles writes.
          if (master) then
            call setup_remd_indices
          end if
        end if ! Replica run setup

        ! Loop over REMD exchanges
        do mdloop = 0, loop

          ! REMD exchange handling.  mdloop == 0 is just used
          ! to get initial energies for the first exchange.
          if (rem < 0 .and. mdloop > 0) then
            call multid_remd_exchange(x, ix, ih, ipairs, qsetup, &
                                      do_list_update, temp0, solvph, solve, &
                                      reservoir_exchange_step)
          else if ((rem == 1 .or. rem == 2) .and. mdloop > 0) then
             if(rremd>0) then
                call reservoir_remd_exchange(1, rem, x(lcrd), x(lvel), x(lmass), &
                               nr3, natom, nr, temp0, reservoir_exchange_step)
             else
                call remd_exchange(1, rem, x(lcrd), x(lvel), x(lmass), &
                               nr3, natom, nr, temp0)
             end if
          else if (rem == 3 .and. mdloop > 0) then
            call hremd_exchange(1, x, ix, ih, ipairs, qsetup, do_list_update)

            ! Force was called inside hremd_exchange, so call nmrdcp
            ! to decrement the NMR counter, since this should not count
            ! as a real step. This is OK, since the counter got
            ! incremented at the _very_ end of nmrcal, so we haven't already
            ! printed an unwanted value (JMS 2/12)
            if (nmropt .ne. 0) then
              call nmrdcp
            end if
          else if (rem == 4 .and. mdloop > 0) then
            call ph_remd_exchange(1, solvph)
          else if (rem == 5 .and. mdloop > 0) then
            call e_remd_exchange(1, temp0, solve)
          end if
          ! End of block for replica exchange handling

#  ifdef VERBOSE_REMD
          if (rem > 0 .and. mdloop .eq. 0 .and. master) then
            write (6,'(a,i4)') "| REMD: Getting initial energy for replica ", &
                               repnum
          end if
#  endif /* VERBOSE_REMD */
#endif  /* MPI */

#if !defined(DISABLE_NFE)
#  ifdef MPI
          call mpi_bcast(infe, 1, mpi_integer, 0, commsander, ier)
#  endif
          if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1 .and. &
              infe == 1) then
            call nfe_on_sander_init(ih, x(lmass), x(lcrd), rem)
          end if
#endif /* DISABLE_NFE */

          ! Branch for calling the Beeman integrator
          if (beeman_integrator == 1) then
            call AM_RUNMD(ix,ih,ipairs, x(lwinv), x(lmass), x, &
                          x(lcrd), x(lvel), x(lforce), qsetup)
          else
            call runmd(x, ix, ih, ipairs, x(lcrd), x(lwinv), x(lmass), &
                       x(lforce), x(lvel), x(lvel2), x(l45), x(lcrdr), &
                       x(l50), x(l95), ix(i70), x(l75), erstop, qsetup)
          end if
          ! End branch for Beeman integrator

#if !defined(DISABLE_NFE)
          if (abfqmmm_param%qmstep == 1 .and. abfqmmm_param%system == 1 .and. &
              infe == 1) then
            call nfe_on_sander_exit()
          end if
#endif /* DISABLE_NFE */
          ! End of Replica Exchange MD block

#ifdef MPI
        end do
        ! End of loop over REMD exchanges

        ! Cleanup REMD files.
        if (rem .ne. 0) then
          call remd_cleanup()
        end if
#endif
        call timer_stop(TIME_RUNMD)
        if (master) then
          call amflsh(6)
        end if

        ! The erstop error condition stems from subroutine shake;
        ! furthermore, it seems that erstop can never be true since shake
        ! can never return with its third last argument, niter, equal to 0.
        if (erstop) then
          if (master) then
            write(6, *) 'FATAL ERROR'
          end if
          call mexit(6,1)
        end if
      case (1)

        ! Minimization: input flag ntmin determines the method of minimization
        call Trace_integer( 'At label ntmincase; ntmin is ', ntmin )
        ntmincase: select case (ntmin)
          case (0, 1, 2)
            call runmin(x, ix, ih, ipairs, x(lcrd), x(lforce), x(lvel), &
                        ix(iibh), ix(ijbh), x(l50), x(lwinv), ix(ibellygp), &
                        x(l95), ene, carrms, qsetup)

            ! If a conventional minimisation is being done,
            ! the restart file is written inside the runmin routine.
          case (LMOD_NTMIN_XMIN)
            write(6,'(a)') '  LMOD XMIN Minimization.'
            write(6,'(a)') ''
            write(6,'(a)') '  Note: Owing to the behaviour of the XMIN &
                            &algorithm,'
            write(6,'(a)') '        coordinates in the trajectory and &
                            &intermediate'
            write(6,'(a)') '        restart files will not match up with &
                            &energies'
            write(6,'(a)') '        in the mdout and mdinfo files. The final &
                            &energy'
            write(6,'(a)') '        and final coordinates do match.'
            write(6,'(a)') ''
            xmin_iter = 0
            call run_xmin(x, ix, ih, ipairs, x(lcrd), x(lforce), &
                          ene, qsetup, xmin_iter, ntpr)
            if (master) then

              ! Write the restart file
              call minrit(0,nrp,ntxo,x(lcrd))
            end if
          case (LMOD_NTMIN_LMOD)
            write(6,'(a)') '  LMOD LMOD Minimization.'
            write(6,'(a)') ''
            write(6,'(a)') '  Note: Owing to the behaviour of the XMIN &
                            &algorithm,'
            write(6,'(a)') '        coordinates in the trajectory and &
                            &intermediate'
            write(6,'(a)') '        restart files will not match up with &
                            &energies'
            write(6,'(a)') '        in the mdout and mdinfo files. The final &
                            &energy'
            write(6,'(a)') '        and final coordinates do match.'
            write(6,'(a)') ''
            call run_lmod(x, ix, ih, ipairs, x(lcrd), x(lforce), ene, qsetup)
            if (master) then

              ! Write the restart file
              call minrit(0, nrp, ntxo, x(lcrd))
            end if
          case default

            ! invalid ntmin: input validation occurs in mdread.f
            ASSERT(.false.)
        end select ntmincase
      case (5)

        ! Modified for reading trajectories (trajene option)
        write (6,*) "POST-PROCESSING OF TRAJECTORY ENERGIES"

        ! Read trajectories and calculate energies for each frame
        call trajene(x, ix, ih, ipairs, ene, ok, qsetup)
        if (.not. ok) then
          write (6,*) 'error in trajene()'
          call mexit(6, 1)
        end if
      case default

        ! Invalid imin: input validation should be transferred to mdread.f
        write(6,'(/2x,a,i3,a)') 'Error: Invalid IMIN (', imin, ').'
        ASSERT(.false.)
    end select imincase

#ifdef MPI /* SOFT CORE */
    if (master) then
      if (icfe .ne. 0 .and. ifsc == 1) then
        call summarize_ti_changes(natom, resat)
      end if
    end if
#endif

    ! Finish up EMAP
    if (temap) then
      call qemap()
    end if
    if (abfqmmm_param%abfqmmm .ne. 1) then
      exit
    else
#ifdef MPI
    call mpi_barrier(commsander,ier)
      if (.not. master) then
        if (charmm_active) then
          call charmm_deallocate_arrays()
        end if
        if (cmap_active) then
          call deallocate_cmap_arrays()
        end if
      end if
#endif /* MPI */
      if (abfqmmm_param%qmstep .ne. abfqmmm_param%maxqmstep .or. &
          abfqmmm_param%system .ne. 2) then
        if (qmmm_nml%ifqnt) then
          call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
        end if
      end if
      call deallocate_m1m2m3()
      if (abfqmmm_param%system == 2 .and. master) then
        call abfqmmm_write_idrst()
        call abfqmmm_write_pdb(x(lcrd),ix(i70))
      end if
      call abfqmmm_next_step()
    end if
  end do abfqmmmmasterdoloop

  if (abfqmmm_param%abfqmmm == 1) then
    deallocate(abfqmmm_param%id, stat=ier)
    REQUIRE(ier == 0)
    deallocate(abfqmmm_param%v, stat=ier)
    REQUIRE(ier == 0)
    deallocate(abfqmmm_param%f, stat=ier)
    REQUIRE(ier == 0)
    deallocate(abfqmmm_param%f1, stat=ier)
    REQUIRE(ier == 0)
    deallocate(abfqmmm_param%f2, stat=ier)
    REQUIRE(ier == 0)
    if (master) then
      deallocate(abfqmmm_param%isqm, stat=ier)
      REQUIRE(ier == 0)
    endif
  end if

  ! Calc time spent running vs setup
  call timer_stop(TIME_TOTAL)
  call wallclock( time1 )
  call date_and_time( final_date, final_time )
#ifdef MPI
  call profile_time(time1 - time0, num_calls_nblist, profile_mpi)
#else
  call profile_time(time1 - time0, num_calls_nblist)
#endif

#ifdef MPI
  ! Set and broadcast notdone in mpi_orig case to inform
  ! other nodes that we are finished calling force().
  if (mpi_orig) then
    notdone = 0
    call mpi_bcast(notdone, 1, mpi_integer, 0, commsander, ier)
  end if
#endif

  ! PUPIL interface: block to finalize the CORBA work
#ifdef PUPIL_SUPPORT
  ! Finalize Corba Interface
  puperror = 0
  call killcorbaintfc(puperror)
  if (puperror /= 0) then
     write(6,*) 'Error ending PUPIL CORBA interface.'
  end if
  write(6,'(a)') 'PUPIL CORBA interface finalized.'
  pupactive = .false.
#endif
  ! End of this PUPIL interface block

  ! Finalize X-ray refinement work
  call amflsh(6)

  if (master) then
#ifdef MPI
    ! Adaptive QM/MM (qmmm_nml%vsolv=2) via multisander:
    ! all groups have identical coordinates and velocities.
    ! Only the master process needs to dump results.
    if ((qmmm_nml%vsolv < 2) .or. (worldrank == 0)) then
#endif
      call close_dump_files
#ifdef MPI
    end if
#endif

    ! Close out external library
    if (iextpot .ne. 0) then
      call external_cleanup
    end if

    ! Close out constant pH work
    if (icnstph .ne. 0 .or. (icnste .ne. 0 .and. cpein_specified)) then
      call cnstph_finalize
    end if

    ! Close out constant Redox potential work
    if (icnste .ne. 0 .and. .not. cpein_specified) then
      call cnste_finalize
    end if

    ! Write out the final timings, taking Replica Exchange MD into account
#ifdef MPI
    if (rem .ne. 0) then
      nstlim_total = nstlim * numexchg
    else
      nstlim_total = nstlim
    end if
#else
    nstlim_total = nstlim
#endif
    if (imin == 0) then
      call print_ongoing_time_summary(nstlim_total, nstlim_total, dt, 6)
    end if
    write(6,'(12(a))') '|           Job began  at ', initial_time(1:2), &
          ':', initial_time(3:4), ':', initial_time(5:10), '  on ',&
          initial_date(5:6), '/', initial_date(7:8), '/', initial_date(1:4)
    write(6,'(12(a))') '|           Setup done at ', setup_end_time(1:2),  &
          ':', setup_end_time(3:4), ':', setup_end_time(5:10), '  on ', &
          setup_end_date(5:6), '/', setup_end_date(7:8), '/', &
          setup_end_date(1:4)
    write(6,'(12(a))') '|           Run   done at ', final_time(1:2),  &
          ':', final_time(3:4), ':', final_time(5:10), '  on ', &
          final_date(5:6), '/', final_date(7:8), '/', final_date(1:4)
    call nwallclock( ncalls )
    write(6, '(''|'',5x,''wallclock() was called'',I8,'' times'')') ncalls
    call amflsh(6)

    if (iesp > 0) then
      call esp(natom, x(lcrd) ,x(linddip))
    end if
  end if
  call amflsh(6)

#ifdef MPI
  999 continue
  ! end if notmastermpi_orig
  ! This is the effective end of the notmastermpi_orig if statement
  ! because the last statement inside that if is a jump here.
#endif

  ! --- dynamic memory deallocation:
  if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%EXTERN .and. master) then
    call qm2_extern_finalize()
  endif
#ifdef QUICK
  if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%ISQUICK) then
    call quick_finalize()
  endif
#endif
#ifdef TCPB
  if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%ISTCPB) then
    call tcpb_finalize()
  endif
#endif
  if (qmmm_nml%ifqnt .and. .not. qmmm_struct%qm_mm_first_call) then

    ! If first_call is still true, this thread never really
    ! called the QMMM routine. E.g. more threads than PIMD replicates
    call deallocate_qmmm(qmmm_nml, qmmm_struct, qmmm_vsolv, qm2_params)
  end if

  if (ipimd > 0) then
    call pimd_finalize(ipimd)
  end if
  if (ineb > 0) then
    call neb_finalize()
  end if
  if (idecomp > 0) then
    call deallocate_real_decomp()
    call deallocate_int_decomp()
  end if
  if (master .and. idecomp == 0) then
    call deallocate_int_decomp()
  end if

#ifdef MPI /* SOFT CORE */
  if (ifsc .ne. 0) then
    call cleanup_sc()
  end if
  if (ifmbar .ne. 0) then
    call cleanup_mbar()
  end if
#endif

  ! Finalize Linear Interaction Energy module if initiated above
  if (ilrt .ne. 0) then
    call cleanup_linear_response(master)
  end if

#ifdef RISMSANDER
  call rism_finalize()
#endif

  if (ifcr .ne. 0) then
    call cr_cleanup()
  end if
  if (sebomd_obj%do_sebomd) then
    if (master) then
      call sebomd_close_files
    end if
    call cleanup_sebomd_arrays
  end if

  if (master .and. iwrap == 2) then
    deallocate(iwrap_mask_atoms, stat=ier)
    REQUIRE(ier == 0)
  end if
  call nblist_deallocate()
  call deallocate_stacks()
  if ((igb /= 0 .and. igb /= 10 .and. ipb == 0) .or. &
      hybridgb > 0 .or. icnstph > 1 .or. icnste > 1) then
    call deallocate_gb()
  end if
  if (master) then
    if (igb == 10 .or. ipb .ne. 0) then
      call pb_free()
    end if
  end if
  deallocate(ih, stat=ier)
  REQUIRE(ier == 0)
  deallocate(ipairs, stat=ier)
  REQUIRE(ier == 0)
  deallocate(ix, stat=ier)
  REQUIRE(ier == 0)
  deallocate(x, stat=ier)
  REQUIRE(ier == 0)
  if (ntb > 0 .and. ifbox == 1 .and. ew_type == 0 .and. mpoltype == 0) then
    call deallocate_m1m2m3()
  end if
  call AMOEBA_deallocate
  call POL_GAUSS_deallocate
  call deallocate_molecule()

  if (abfqmmm_param%abfqmmm .ne. 1) then
    if (charmm_active) then
      call charmm_deallocate_arrays()
    end if
    if (cmap_active) then
      call deallocate_cmap_arrays()
    end if
  end if
  call Trace_exit( 'sander' )
  if (master .and. mdout .ne. 'stdout') then
    close(6)
  end if

  return

end subroutine sander

!------------------------------------------------------------------------------
! esp: Calculate the ElectroStatic Potential for a system.
!
! Arguments:
!   natom:     the number of atoms in the system
!   x:         array of coordinates for the atoms
!   mom_ind:   array of induced moments on each atom
!------------------------------------------------------------------------------
subroutine esp(natom, x, mom_ind)

  ! Routine to calculate the ESP due to the induced moments (and only
  ! those induced moments) at the same spatial points as the reference QM.
  use constants, only : zero, BOHRS_TO_A, INV_AMBER_ELECTROSTATIC
  use file_io_dat
  implicit none
  integer natom
  _REAL_  x(3,*)
  _REAL_  mom_ind(3,*)

#  include "ew_mpole.h"

  integer dat_unit, new_unit, minus_new_unit
  parameter(dat_unit=30, new_unit=31, minus_new_unit=33)

  integer inat, nesp, idum
  _REAL_  xin, yin, zin
  integer jn, kn
  _REAL_  esp_qm, xb_esp, yb_esp, zb_esp
  _REAL_  x_esp, y_esp, z_esp
  _REAL_  e_x, e_y, e_z, e_q, esp_new
  _REAL_  dist, dist3

  call amopen(dat_unit,"esp.dat",'O','F','R')
  call amopen(new_unit,"esp.induced",owrite,'F','W')
  call amopen(minus_new_unit,"esp.qm-induced",owrite,'F','W')
  read (dat_unit,'(3i5)')inat,nesp,idum
  write(6,'(t2,''inat = '',i5)')inat
  write(6,'(t2,''nesp = '',i5)')nesp

  write(new_unit,'(2i5)')inat,nesp
  write(minus_new_unit,'(2i5)')inat,nesp

  if (inat /= natom) then
    write(6,'(t2,''natom mismatch with esp file'')')
    call mexit(6,1)
  end if

  do jn = 1,inat
    read (dat_unit,'(17x,3e16.7)')xin,yin,zin
    write(new_unit,'(17x,3e16.7)')xin,yin,zin
    write(minus_new_unit,'(17x,3e16.7)')xin,yin,zin
  end do

  do jn = 1,nesp
    e_x = zero
    e_y = zero
    e_z = zero
    e_q = zero
    read(dat_unit,'(1x,4e16.0)') esp_qm, xb_esp, yb_esp, zb_esp
    x_esp = xb_esp * BOHRS_TO_A
    y_esp = yb_esp * BOHRS_TO_A
    z_esp = zb_esp * BOHRS_TO_A
    do kn = 1, natom
      dist = sqrt((x(1,kn)-x_esp)**2 + (x(2,kn)-y_esp)**2 + &
                  (x(3,kn)-z_esp)**2)
      dist3 = dist**3
      e_x = e_x - mom_ind(1,kn)*(x(1,kn) - x_esp)/dist3
      e_y = e_y - mom_ind(2,kn)*(x(2,kn) - y_esp)/dist3
      e_z = e_z - mom_ind(3,kn)*(x(3,kn) - z_esp)/dist3
    end do

    e_x = e_x * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
    e_y = e_y * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
    e_z = e_z * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
    e_q = e_q * BOHRS_TO_A * INV_AMBER_ELECTROSTATIC
    esp_new = e_x + e_y + e_z
    write(new_unit, '(1x,4e16.7)') esp_new, xb_esp, yb_esp, zb_esp
    write(minus_new_unit, '(1x,4e16.7)') esp_qm-esp_new, &
                                           xb_esp, yb_esp, zb_esp
  end do

  close(dat_unit)
  close(new_unit)
  close(minus_new_unit)

  return

end subroutine esp
