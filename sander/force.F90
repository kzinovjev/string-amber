! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"
#include "nfe-config.h"

!------------------------------------------------------------------------------
! force: The main driver routine to compute energies and forces in sander.
!        This is a simple nexus of many different modules, each with their own
!        contributions to how the atoms will move next.
!
! Arguments:
!   xx:             global real array (holding e.g. coordinates at position
!                   l15 etc., see locmem.F90)
!   ix:             global array holding all integer arrays (see locmem.F90)
!   ih:             Hollerith array containing atom names, types,
!                   residue names, and more
!   ipairs:         ?? Global pairlist ?? -- needs explanation according to JMS
!   x:              Coordinates of all atoms
!   f:              Forces on all atoms
!   ener:           Energy with components (see the state_rec type in the
!                   state.F90 module)
!   vir:            Virial (four element _REAL_ vector)
!   fs:
!   rborn:          The Generalized Born radii (they will be recalculated
!                   within this subroutine)
!   reff:           The effective Born radii (recalculated in this subroutine)
!   onereff:        The inverse effective Born radii, 1.0/reff
!   qsetup:         Flag to activate setup of multiple components, .false. on
!                   first call
!   do_list_update: flag to have the non-bonded list updated (returns .TRUE. or
!                   .FALSE. to the calling subroutine following a call to
!                   nonbond_list)
!   nstep:          The step number
!------------------------------------------------------------------------------
subroutine force(xx, ix, ih, ipairs, x, f, ener, vir, fs, rborn, reff, &
                 onereff, qsetup, do_list_update, nstep)

#if !defined(DISABLE_NFE)
  use nfe_sander_hooks, only: nfe_on_force => on_force
  use nfe_sander_proxy, only: infe
#endif /* DISABLE_NFE */
  use file_io_dat
#ifdef LES
  use genbornles
  use pimd_vars, only: nrg_all
#  ifdef MPI
  use les_data, only: elesa, elesb, elesd, elesp
#  endif /* MPI */
#else
  use genborn
  use qmmm_module, only: qmmm_struct, qm2_struct
#endif /* LES */
#ifdef MPI
  use neb_vars, only: ineb, neb_force
  use full_pimd_vars, only: totener
  use softcore, only: sc_ener
#endif /* MPI */
#if defined(LES) && defined(MPI)
  use evb_data, only: nrg_frc
  use pimd_vars, only: equal_part
  use miller, only: dlnQ_dl
  use remd, only: rem ! wasn't used for LES above
#endif /* LES && MPI */
  use poisson_boltzmann, only: pb_force
  use dispersion_cavity, only: npopt, np_force
  use pbtimer_module, only: pbtimer_init, pbtimer_summary
#ifdef RISMSANDER
  use sander_rism_interface, only: rismprm, rism_force
#endif /* RISMSANDER */
#ifdef APBS
  use apbs
#endif /* APBS */
  use trace
  use stack
  use pimd_vars, only: ipimd, nbead, bnd_vir, Epot_spring, Epot_deriv, &
                       real_mass, itimass
  use qmmm_module, only : qmmm_nml
  use constants, only: zero, one
  use relax_mat
  use ew_recip
  use parms, only: cn1, cn2, cn6, asol, bsol, pk, rk, tk, numbnd, numang, &
                   nptra, nphb, nimprp, cn3, cn4, cn5 ! for another vdw model
#ifdef PUPIL_SUPPORT
  use nblist, only: nonbond_list, a, b, c, alpha, beta, gamma, ucell
#else
  use nblist, only: nonbond_list, a, b, c, alpha, beta, gamma
#endif /*PUPIL_SUPPORT*/
#ifdef DSSP
  use dssp, only: fdssp, edssp, idssp
#endif /* DSSP */

  ! AMOEBA modifications
  use amoeba_interface, only: AM_VAL_eval, AM_NonBond_eval
  use amoeba_mdin, only : iamoeba,am_nbead

  ! pGM additions
  use pol_gauss_interface, only: pGM_NonBond_eval
  use pol_gauss_mdin, only : ipgm

  use amd_mod
  use scaledMD_mod
  use nbips, only: ips, eexips
  use emap, only: temap, emapforce

#ifdef PUPIL_SUPPORT
  use pupildata
#endif /*PUPIL_SUPPORT*/

  use linear_response, only: ilrt, ee_linear_response, energy_m0, energy_w0, &
                             energy_vdw0, cn1_lrt, cn2_lrt, crg_m0, crg_w0, &
                             do_lrt, f_scratch, lrt_solute_sasa
  use cns_xref
  use xray_interface_module, only: xray_get_derivative, xray_active
#ifdef USE_ISCALE
  use xray_globals_module, only: atom_bfactor
#endif

  ! CHARMM Force Field Support
  use charmm_mod, only: charmm_active, charmm_calc_impropers, &
                        charmm_calc_cmap, charmm_calc_urey_bradley,&
                        charmm_dump_gold, &
                        do_charmm_dump_gold
  use ff11_mod, only: cmap_active, calc_cmap
#ifdef MPI
  use decomp, only: collect_dec2, init_dec
#else
  use decomp, only: init_dec
#endif /* MPI */
  use state
  use crg_reloc, only: ifcr, cr_reassign_charge, cr_calc_force
  use sebomd_module, only: sebomd_obj, sebomd_save_forces
  use abfqmmm_module
  use les_data, only: temp0les
  use music_module, only : music_force
  use external_module, only : pme_external, gb_external

  implicit none

  integer, intent(in) :: nstep

#ifdef PUPIL_SUPPORT
  character(kind=1,len=5) :: routine="force"
#endif
#if defined(PUPIL_SUPPORT) || defined(MPI)
  integer ierr
#endif
  integer   ipairs(*)
  _REAL_ xx(*)
  integer   ix(*)
  character(len=4) ih(*)
  _REAL_ fs(*), rborn(*), reff(*), dvdl
  _REAL_, intent(out) :: onereff(*)
#include "def_time.h"
#include "ew_frc.h"
#include "ew_cntrl.h"
#include "extra_pts.h"
#include "parallel.h"

#ifdef MPI
#  include "ew_parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
  include 'mpif.h'
  integer gb_rad_mpistart, j3, j, i3
  _REAL_ :: temp_amd_totdih
#  ifdef LES
  _REAL_ :: vel0_nrg_sum
#  endif /* LES */
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
#endif /* MPI */

  logical belly
#include "../include/md.h"
#include "../pbsa/pb_md.h"
#include "box.h"
#include "nmr.h"
#include "../include/memory.h"
#include "extra.h"
#include "tgtmd.h"
#include "multitmd.h"
#include "flocntrl.h"
  integer istart, iend
  _REAL_ evdwex, eelex
  _REAL_ enemap
  _REAL_ xray_e

  logical, intent(inout) :: qsetup
  logical, intent(out) :: do_list_update

  _REAL_  enmr(6), devdis(4), devang(4), devtor(4), devpln(4), devplpt(4), &
          devgendis(4), entr, ecap, enfe
  _REAL_  x(*), f(*), vir(4)
  type(state_rec)  ener

  ! Local
  _REAL_                     :: ene(30)    !Used locally ONLY
  type(potential_energy_rec) :: pot        !Used locally ONLY
#ifdef USE_ISCALE
  logical, save :: first=.true.
#endif

#if defined(LES) && defined(MPI)
  _REAL_  :: nrg_bead(nbead)
#endif /* LES && MPI */

#ifndef LES
  _REAL_ escf
#endif /* LES */

  integer i
  _REAL_  virvsene, eelt, epol, esurf, edisp
#ifdef APBS
  _REAL_ enpol
#endif /* APBS */
#ifdef RISMSANDER
  _REAL_ erism
#endif /*RISMSANDER*/

  ! Charge transfer
  _REAL_ ect

  _REAL_ epolar, aveper, aveind, avetot, emtot, dipiter, dipole_temp
  integer, save :: newbalance

  ! Accelerated MD variables
  _REAL_ amd_totdih

  ! SEBOMD Gradient test variables
  _REAL_ :: em2, em1, e0, ep1, ep2
  _REAL_ :: fx1, fx2, fxx, xdx, dx

  ! MuSiC
  _REAL_ :: music_vdisp, music_vang, music_vgauss, music_spohr89

  ect = 0.0

  call trace_enter( 'force' )

  call timer_start(TIME_FORCE)
  if ( idecomp /= 0 .and. icfe == 0) then
    call init_dec
  end if
  ene(:) = ZERO
  call zero_pot_energy(pot)

  belly = ibelly > 0

#ifdef MPI
  if (mpi_orig) then

    ! Check to see if we are done yet in mpi_orig case (tec3).
    ! This is done by monitoring the status of an integer notdone.
    ! If notdone .eq. 1 then we keep going.  notdone is set to zero
    ! when we no longer want to call force().  This perhaps is not the
    ! most efficient means to implement this check...
    call mpi_bcast(notdone,1,mpi_integer,0,commsander,ierr)
    if (notdone /= 1) return

    ! Send copies of xyz coords, setbox common block, vir array
    ! and NTNB value to all nodes from master with a broadcast.
    if (numtasks > 1) then
      call mpi_bcast(box, BC_BOXR, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(ntb, BC_BOXI, mpi_integer, 0, commsander, ierr)
      call mpi_bcast(vir, 3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
      call mpi_bcast(xx(lcrd), 3*natom, MPI_DOUBLE_PRECISION, 0, &
                     commsander, ierr)
      call mpi_bcast(ntnb, 1, mpi_integer, 0, commsander, ierr)
      if (iabs(ntb) >= 2) then
        call mpi_bcast(xx(l45), 3*natom, MPI_DOUBLE_PRECISION, &
                       0, commsander, ierr)
      end if
    end if
  end if
  istart = iparpt(mytaskid) + 1
  iend = iparpt(mytaskid+1)
#else
  istart = 1
  iend = natom
#endif /* MPI */

  ! QM/MM variable QM solvent scheme
  if (qmmm_nml%ifqnt .and. qmmm_nml%vsolv > 0) then
    call timer_start(TIME_QMMM)
    call timer_start(TIME_QMMMVARIABLESOLVCALC)

    ! For the moment this is NOT parallel: all threads need to call this.
    call qmmm_vsolv_update(nstep, natom, qsetup, xx(l15), x,  &
                           nbonh, nbona, ntheth, ntheta, nphih, nphia, &
                           ix(iibh), ix(ijbh), ix(iicbh), &
                           ix(iiba), ix(ijba), ix(iicba), &
                           ix(i24), ix(i26), ix(i28), ix(i30), &
                           ix(i32), ix(i34), ix(i36), ix(i38), &
                           ix(i40), ix(i42), ix(i44), ix(i46), &
                           ix(i48), ix(i50), ix(i52), ix(i54), &
                           ix(i56), ix(i58), ix(ibellygp))
    call timer_stop(TIME_QMMMVARIABLESOLVCALC)
    call timer_stop(TIME_QMMM)
  end if
  ! End QM/MM variable QM solvent scheme

  if (iamoeba .eq. 1) then
    REQUIRE(am_nbead .eq. ncopy)
  end if

  ! Zero out the energies and forces
  enoe = 0.d0
  aveper = 0.d0
  aveind = 0.d0
  avetot = 0.d0
  dipiter = 0.d0
  dvdl = 0.d0
  dipole_temp = 0.d0
  enmr(1:6) = 0.d0
  enfe = 0.d0
  vir(1:4) = 0.d0
  virvsene = 0.d0
  f(1:3*natom+iscale) = 0.d0
#ifdef LES
   if (ipimd > 0) then
     nrg_all(1:nbead) = 0.d0
   end if
#endif

  if (sebomd_obj%do_sebomd) then
    ! step 0: initialize temporary array for saving restraint
    !         (required since SEBOMD forces are stored in a different array)
    call sebomd_save_forces(0, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
  end if

#ifndef PUPIL_SUPPORT
  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then

    ! For GB: do all nonbondeds together below
    call timer_start(TIME_NONBON)
    call timer_start(TIME_LIST)
    if (abfqmmm_param%abfqmmm == 1) then
      qsetup = .true.
    end if
    call nonbond_list(x, ix(i04), ix(i06), ix(i08), ix(i10), ntypes, &
                      natom/am_nbead, xx, ix, ipairs, ntnb, ix(ibellygp), &
                      belly, newbalance, qsetup, do_list_update)
    call timer_stop(TIME_LIST)
    call timer_stop(TIME_NONBON)
  end if
  ! charge reassign here !
  if (ifcr /= 0) then
    call cr_reassign_charge(x, f, pot%ct, xx(l15), natom)
  end if
#endif
  if (sebomd_obj%do_sebomd) then
    ! step 1/2 to save forces from nfe_on_force
    call sebomd_save_forces(1, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
  end if

#if !defined(DISABLE_NFE)
  if (infe == 1) then
    call nfe_on_force(x, f, enfe)
  end if
#endif

  ! Semi-Empirical Born-Oppenheimer MD
  if (sebomd_obj%do_sebomd) then
#include "sebomd_force.inc"
  endif

  ! Do weight changes, if requested.  Updated 9/2007 by Matthew Seetin
  ! to enable plane-point and plane-plane restraints.
  if (nmropt > 0) then
    call nmrcal(x, f, ih(m04), ih(m02), ix(i02), xx(lwinv), enmr, devdis, &
                devang, devtor, devplpt, devpln, devgendis, temp0, tautp, &
                cut, xx(lnmr01), ix(inmr02), xx(l95), 31, 6, rk, tk, pk, cn1, &
                cn2, asol, bsol, xx(l15), numbnd, numang, nptra-nimprp, &
                nimprp, nphb, natom, natom, ntypes, nres, rad, wel, radhb, &
                welhb, rwell, tgtrmsd, temp0les, -1, 'WEIT')
  end if

! If calling an External library
if (iextpot .gt. 0) then
#ifdef MPI
  if (sanderrank .eq. 0) then
#endif /* MPI */
    if (igb == 0) then
      call pme_external(x, f, ener%pot%tot)
    else
      call gb_external(x, f, ener%pot%tot)
    endif
    if (nmropt > 0) then
      call nmrcal(x, f, ih(m04), ih(m02), ix(i02), xx(lwinv), enmr, devdis, &
                  devang, devtor, devplpt, devpln, devgendis, temp0, tautp, &
                  cut, xx(lnmr01), ix(inmr02), xx(l95), 31, 6, rk, tk, pk, cn1, &
                  cn2, asol, bsol, xx(l15), numbnd, numang, nptra-nimprp, &
                  nimprp, nphb, natom, natom, ntypes, nres, rad, wel, radhb, &
                  welhb, rwell, tgtrmsd, temp0les, -1, 'CALC')
    end if
    if (natc > 0 .and. ntr==1) then   ! ntr=1 (positional restraints)
      call xconst(natc, entr, ix(icnstrgp), x, f, xx(lcrdr), xx(l60))
    end if
    ener%pot%constraint = + sum(enmr(1:6)) + entr + enfe
    ener%pot%tot = ener%pot%tot + sum(enmr(1:6)) + entr + enfe
#ifdef MPI
  endif

  call mpi_bcast(f, 3*natom, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
#endif /* MPI */
! If not calling an External library
else

  epolar = 0.d0

  ! EGB: if Generalized Born is in effect (there is a GB solvent,
  ! not just electrostatics in a vacuum or some Poisson-Boltzmann
  ! solvent), then we need to calculate the GB radii for this
  ! structure.  If we are doing GB in the context of QM (that is,
  ! qm_gb = 2) then we need to calculate the GB radii before
  ! calling qm_mm.
  if (igb > 0 .and. igb /= 6 .and. igb /= 10 .and. ipb == 0 .and. &
      (irespa < 2 .or. mod(irespa,nrespai) == 0)) then
#ifdef MPI
    gb_rad_mpistart = mytaskid+1
#endif
    call timer_start(TIME_EGB)
    call timer_start(TIME_GBRAD1)

    ! In QM/MM, this will calculate us the radii for the link atoms
    ! and not the MM link pair atoms.  The initial radii used are
    ! those of the mm link pair's atom type though.  The MM link
    ! pair atoms must be the coordinates of the link atoms here.
    if (qmmm_nml%ifqnt) then
      call adj_mm_link_pair_crd(x)
    end if

#ifdef LES
#  ifdef MPI
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        gbalpha, gbbeta, gbgamma, gb_rad_mpistart)
#  else
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        gbalpha, gbbeta, gbgamma)
#  endif
#else
#  ifdef MPI
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        xx(l2402), xx(l2403), xx(l2404), gb_rad_mpistart)
#  else
    call egb_calc_radii(igb, natom, x, fs, reff, onereff, fsmax, rgbmax, &
                        rborn, offset, rbornstat, xx(l188), xx(l189), &
                        xx(l186), xx(l187), gbneckscale, ncopy, rdt, &
                        xx(l2402), xx(l2403), xx(l2404))
#  endif
#endif
    if (qmmm_nml%ifqnt) then
      call rst_mm_link_pair_crd(x)
    end if
    call timer_stop(TIME_GBRAD1)
    call timer_stop(TIME_EGB)
  end if
  ! End EGB

  ! QM/MM Contributions are now calculated before the NON-Bond info.
  if (qmmm_nml%ifqnt) then

    ! If we are doing periodic boundaries with QM/MM PME then we need to
    ! do the PME calculation twice. First here to get the potential at
    ! each QM atom due to the PME and then again after all the QM is done
    ! to get the MM-MM potential and all of the gradients.
    if (qmmm_nml%qm_pme) then
      ! Ewald force will put the potential into the qm_ewald%mmpot array.
      call timer_start(TIME_EWALD)
      call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1, cn2, cn6, &
                       eelt, epolar, f, xx, ix, ipairs, xx(l45), virvsene, &
                       xx(lpol), &
#ifdef HAS_10_12
                       xx(lpol2), .true., cn3, cn4, cn5, asol, bsol)
#else
                       xx(lpol2), .true., cn3, cn4, cn5 )
#endif
      call timer_stop(TIME_EWALD)
    endif

    call timer_start(TIME_QMMM)
#ifndef LES
    call qm_mm(x, natom, qmmm_struct%scaled_mm_charges, f, escf, periodic, &
               reff, onereff, intdiel, extdiel, Arad, cut, &
               qm2_struct%scf_mchg, ntypes, ih(m04), ih(m06), xx(lmass), &
               ix(i04), nstep)
    pot%scf = escf
#endif
    call timer_stop(TIME_QMMM)
  end if
  !END qm/mm contributions, triggered by ifqnt in the qmmm_nml namelist

#include "pupil_force.inc"

  if (iamoeba .eq. 1) vir(1:4)=0.0

  ! Calculate the non-bonded contributions
  call timer_start(TIME_NONBON)

  ! EEXIPS: three dimensional Isotropic Periodic Sums.  This routine
  ! initializes electrostatic and vdW energies in IPS and updates
  ! the IPS parameters.
  call timer_start(TIME_EEXIPS)
  if (ips > 0) then
    if (iamoeba == 0) then
      call eexips(evdwex, eelex, istart, iend, ntb, ntypes, ix(i04), &
                  ix(i06), ix(i08), ix(i10), xx(l15), cn1, cn2, f, x)
    else
      evdwex = 0.0d0
      eelex = 0.0d0
    endif
  endif
  call timer_stop(TIME_EEXIPS)

  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then

    ! (for GB: do all nonbondeds together below)
    call timer_start(TIME_EWALD)
    if (iamoeba == 1) then
      call AM_NonBond_eval(natom, x, f, vir, xx, ipairs, evdw, eelt, epolar, &
                           enb14, ee14, diprms, dipiter)
    else if (ipgm == 1) then
       call pGM_NonBond_eval(natom, x, f, vir, xx, ipairs, &
                             ntypes, ix(i04), ix(i06), cn1, cn2, & ! R. Luo: amber vdw
                             evdw, eelt, epolar, enb14, ee14, &
                             diprms, dipiter)
    else
      if (induced > 0) then
        call handle_induced(x, natom, ix(i04), ix(i06), xx(l15), cn1, cn2, &
                            cn6, eelt, epolar, f, xx, ix, ipairs, xx(lpol), &
                            xx(lpol2), xx(lpolbnd), xx(l45), virvsene, &
                            ix(i02), ibgwat, nres, aveper, aveind, avetot, &
                            emtot, diprms, dipiter, dipole_temp, &
#ifdef HAS_10_12
                            cn3, cn4, cn5, asol, bsol)
#else
                            cn3, cn4, cn5)
#endif
      else
        if (ilrt /= 0) then

          ! Modifications for computing interaction energy
          ! according to the Linear Response Theory, LIE module
          if (do_lrt) then

            ! call with molecule charges set to zero
            call ewald_force(x, natom, ix(i04), ix(i06), crg_m0, cn1, cn2, &
                             cn6, energy_m0, epolar, f_scratch, xx, ix, &
                             ipairs, xx(l45), virvsene, xx(lpol), &
#ifdef HAS_10_12
                             xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                             xx(lpol2), .false. , cn3, cn4, cn5)
#endif

            ! call with water charges set to zero
            call ewald_force(x, natom, ix(i04), ix(i06), crg_w0, cn1, cn2, &
                             cn6, energy_w0, epolar, f_scratch, xx, ix, &
                             ipairs, xx(l45), virvsene, xx(lpol), &
#ifdef HAS_10_12
                             xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                             xx(lpol2), .false. , cn3, cn4, cn5)
#endif
            ! call with full charges but no vdw interaction
            ! between solute and solvent
            call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1_lrt, &
                             cn2_lrt, cn6, eelt, epolar, f_scratch, xx, ix, &
                             ipairs, xx(l45), virvsene, xx(lpol), &
#ifdef HAS_10_12
                             xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                             xx(lpol2), .false. , cn3, cn4, cn5)
#endif
            energy_vdw0 = evdw
            call lrt_solute_sasa(x,natom, xx(l165))
          end if

          ! call normal_ewald force this will overwrite everything
          ! computed above except energy_m0 and energy_w0
          call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1, cn2, &
                           cn6, eelt, epolar, f, xx, ix, ipairs, xx(l45), &
                           virvsene, xx(lpol), &
#ifdef HAS_10_12
                           xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                           xx(lpol2), .false. , cn3, cn4, cn5)
#endif
          energy_vdw0 = evdw - energy_vdw0

          ! count call to ltr, maybe calculate Eee and print it
          call ee_linear_response(eelt, master)
        else ! just call ewald_force normally
          call ewald_force(x, natom, ix(i04), ix(i06), xx(l15), cn1, cn2, &
                           cn6, eelt, epolar, f, xx, ix, ipairs, xx(l45), &
                           virvsene, xx(lpol), &
#ifdef HAS_10_12
                           xx(lpol2), .false. , cn3, cn4, cn5, asol, bsol)
#else
                           xx(lpol2), .false. , cn3, cn4, cn5)
#endif
        end if ! ilrt /= 0
      end if ! induced > 0
    end if ! iamoeba == 1

    call timer_stop(TIME_EWALD)
#ifdef MPI
      if (mytaskid == 0) then
#endif
        pot%vdw   = evdw
        pot%elec  = eelt
        pot%hbond = ehb  !whereis ehb?
#ifdef MPI
      else
        ! energies have already been reduced to the master
        ! node in ewald_force, so here we zero out elements
        ! on non-master nodes:
        pot%vdw      = 0.d0
        pot%elec     = 0.d0
        pot%hbond    = 0.d0
      end if
#endif
      if( ips > 0 )then
         pot%vdw   = pot%vdw   + evdwex
         pot%elec  = pot%elec  + eelex
      endif
   end if  ! ( igb == 0 .and. ipb == 0 .and. iyammp == 0 )

   ! End of non-bonded computations
   call timer_stop(TIME_NONBON)

   ! Calculate other contributions to the forces
   !   - When igb==10, a Poisson-Boltzmann solvent is active.  All nonbonds
   !     are done in the subroutine pb_force, and all nonpolar interactions
   !     are done in np_force:
   !   - Holger Gohlka put this part here such that "outflag" is known from
   !     a call of pb_force; outflag is needed in the "bond" routine in the
   !     case of ifcap == 2,5 (i.e., ivcap == 1,5)
#ifdef MPI
   if(mytaskid == 0)then
#endif
     if (igb == 10 .or. ipb /= 0 ) then
       call timer_start(TIME_PBFORCE)
       call pbtimer_init
       call pb_force(natom, nres, ntypes, npdec, ix(i02), ix(i04), ix(i06), &
                     ix(i10), cn1, cn2, xx(l15), x, f, evdw, eelt, epol)
       if (pbgrid) pbgrid = .false.
       if (pbinit) pbinit = .false.
       pot%vdw  = evdw
       pot%elec = eelt
       pot%pb   = epol
       call timer_stop(TIME_PBFORCE)

       call timer_start(TIME_NPFORCE)
       esurf = 0.0d0
       edisp = 0.0d0
       if (ifcap == 0  .and. npopt /= 0) then
         call np_force(natom, nres, ntypes, ix(i02), ix(i04), ix(i06), &
                       cn1, cn2, x, f, esurf, edisp)
       end if
       if (pbprint) pbprint = .false.
       pot%surf = esurf
       pot%disp = edisp
       call pbtimer_summary
       call timer_stop(TIME_NPFORCE)
     end if  ! ( igb == 10 .or. ipb /= 0 )
#ifdef MPI
   end if
#endif

  ! Bonds with H
  call timer_start(TIME_BOND)

  ! initialize bond virial
  if (ipimd > 0) bnd_vir = zero

#ifdef MPI /* SOFT CORE */
  ! zero only once, sc bond energy is sum of H and non-H terms
  sc_ener(1) = 0.0d0
#endif
  if (ntf < 2) then
    ebdev = 0.d0
    call bond(nbonh, ix(iibh), ix(ijbh), ix(iicbh), x, f, ene(6))
    pot%bond = pot%bond + ene(6)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesb
    endif
#  endif
#endif
  end if

  ! Bonds without H
  if (ntf < 3) then
    call bond(nbona+nbper, ix(iiba), ix(ijba), ix(iicba), x, f, ene(7))
    pot%bond = pot%bond + ene(7)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesb
    endif
#  endif
#endif
    if (nbonh+nbona > 0) then
      ebdev = sqrt(ebdev/(nbonh+nbona))
    end if
  end if

  ! Angles with H
  if (ntf < 4) then

#ifdef MPI /* SOFT CORE */
    ! zero only once, sc bond energy is sum of H and non-H terms
    sc_ener(2) = 0.0d0
#endif
    eadev = 0.d0

    call angl(ntheth, ix(i24), ix(i26), ix(i28), ix(i30), x, f, ene(8))
    pot%angle = pot%angle + ene(8)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesa
    endif
#  endif
#endif
  end if

  ! Angles without H
  if (ntf < 5) then
    call angl(ntheta+ngper, ix(i32), ix(i34), ix(i36), ix(i38), x, f, &
              ene(9))
    pot%angle = pot%angle + ene(9)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesa
    end if
#  endif
#endif
    if (ntheth+ntheta > 0) then
      eadev = 57.296*sqrt( eadev/(ntheth+ntheta) )
    end if
  end if

  ! Added by Romelia Salomon: for Accelerated MD, calculate
  ! dihedral energy first to estimate dihedral weight, then
  ! use it in the regular ephi function.  Fix this later, as
  ! AMD dihedral weight does NOT support AMOEBA
  if (iamd .gt. 1) then
    ! Dihedrals with H
    if (ntf < 6) then
      call ephi_ene_amd(nphih, ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), &
                        x, amd_dih_noH)
    endif
    ! Dihedrals without H
    if (ntf < 7) then
      call ephi_ene_amd(nphia+ndper, ix(i50), ix(i52), ix(i54), ix(i56), &
                        ix(i58), x, amd_dih_H)
    endif
#ifdef MPI
    temp_amd_totdih = amd_dih_noH + amd_dih_H
#  ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE, temp_amd_totdih, 1, &
                       MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
    amd_totdih = temp_amd_totdih
#  else
    call mpi_allreduce(temp_amd_totdih, amd_totdih, 1, MPI_DOUBLE_PRECISION, &
                       mpi_sum, commsander, ierr)
#  endif /* USE_MPI_IN_PLACE */
#else
    amd_totdih = amd_dih_noH + amd_dih_H
#endif /* MPI */
    call calculate_amd_dih_weights(amd_totdih)
  endif

  ! Dihedrals with H
  ! initialize 14 nb energy virial
  if (ipimd>0) then
    e14vir = zero
  end if
  if (ntf < 6) then

#ifdef MPI /* SOFT CORE */
    ! zero only once, sc bond energy is sum of H and non-H terms
    sc_ener(3) = 0.0d0
#endif
    call ephi(nphih, ix(i40), ix(i42), ix(i44), ix(i46), ix(i48), xx(l15), &
              ix(i04), x, f, dvdl, ene(10), ene(11), ene(12), xx(l190))

    ! Combine contributions from dihedrals with H, including
    ! the torsions themselves and 1:4 vdW / elec interactions.
    pot%dihedral = pot%dihedral + ene(10)
    pot%vdw_14   = pot%vdw_14   + ene(11)
    pot%elec_14  = pot%elec_14  + ene(12)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesd
    endif
#  endif
#endif
  end if

  ! Dihedrals without H
  if (ntf < 7) then
    call ephi(nphia+ndper, ix(i50), ix(i52), ix(i54), ix(i56), ix(i58), &
              xx(l15), ix(i04), x, f, dvdl, ene(13), ene(14), ene(15), &
              xx(l190))

    ! Combine contributions from dihedrals without H, including
    ! the torsions themselves and 1:4 vdW / elec interactions.
    pot%dihedral = pot%dihedral + ene(13)
    pot%vdw_14   = pot%vdw_14   + ene(14)
    pot%elec_14  = pot%elec_14  + ene(15)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
       pot%les = pot%les + elesd
    endif
#  endif
#endif

    ! Do CHARMM impropers, if the CHARMM force field is in use
    ! Note: CHARMM does not distinguish between impropers with and
    ! without hydrogen.  Hence, it is not possible to strictly
    ! conform to the ntf options of sander. Here CHARMM impropers
    ! are calculated as long as ntf < 7 - so CHARMM impropers are
    ! essentially considered to be in the same set as dihedrals
    ! NOT involving hydrogen.
    if (charmm_active) then
      call charmm_calc_urey_bradley(x, pot%angle_ub, f)
      call charmm_calc_impropers(x, pot%imp, f)
      call charmm_calc_cmap(x, pot%cmap, f)
    end if
    ! End CHARMM impropers computation

    if (cmap_active) then
      call calc_cmap(x,pot%cmap,f)
    end if
  end if
  ! End computations for ntf < 7

  ! AMOEBA valence terms: bond, angle, torsion
  if (iamoeba == 1) then
    call AM_VAL_eval(x, f, vir, ene(6), ene(8), ene(10))
    pot%bond     = pot%bond + ene(6)
    pot%angle    = pot%angle + ene(8)
    pot%dihedral = pot%dihedral + ene(10)
  end if

  ! End of valence contribution computations
  call timer_stop(TIME_BOND)

  ! Step 1/2 to save emap + entr + ecap forces
  if (sebomd_obj%do_sebomd) then
     call sebomd_save_forces(1, natom, f, xx(gradsebomd), xx(grad1tmp), &
                             xx(grad2tmp))
  end if

  ! Calculate the EMAP constraint energy
  if (temap) then   ! ntr=1 (positional restraints)
    call emapforce(natom, enemap, xx(lmass), x, f)
    pot%emap = enemap
  end if

  ! Calculate the position constraint energy
#ifdef MPI /* SOFT CORE */
  ! Zero all restraint/constraint energies
  sc_ener(14:19) = 0.d0
#endif
  if (natc > 0 .and. ntr==1) then   ! ntr=1 (positional restraints)
    call xconst(natc, entr, ix(icnstrgp), x, f, xx(lcrdr), xx(l60))
    pot%constraint = entr
  end if
  if (itgtmd==1 .and. (nattgtfit > 0 .or. nattgtrms > 0)) then

    ! Calculate rmsd for targeted md (or minimization) if requested.
    ! All nodes do rms fit, could just be master then broadcast.
    ! All nodes need all coordinates for this.
    call rmsfit(xx(lcrdr), x, xx(lmass), ix(itgtfitgp), ix(itgtrmsgp), &
                rmsdvalue, nattgtrms, nattgtfit, rmsok)
    if (.not. rmsok) then
      if (master) then
        write (6,*) 'Fatal error: Error calculating RMSD!'
      end if
      call mexit(6, 1)
    end if

    call xtgtmd(entr, ix(itgtrmsgp), x, f, xx(lcrdr), xx(lmass), tgtrmsd, &
                tgtmdfrc, rmsdvalue, nattgtrms)
    pot%constraint = entr
  else if (itgtmd == 2) then
    call mtmdcall(entr, xx(lmtmd01), ix(imtmd02), x, f, ih(m04), ih(m02), &
                  ix(i02), ih(m06), xx(lmass), natom, nres, 'CALC')
    pot%constraint = entr
  end if

  if (ifcap == 1 .or. ifcap == 2) then
    call capwat(natom, x, f, ecap)
    pot%constraint = pot%constraint + ecap
  else if (ifcap == 3) then
    write(6,*) 'No energy expression for spherical boundary known yet'
    call mexit(6,1)
  else if(ifcap == 4) then
    write(6,*) 'No energy expression for orthorhombic boundary known yet'
    call mexit(6,1)
    !call orth(natom,ix(ibellygp),x,f,eorth)
    !ene(20) = ene(20) + eorth
  end if

  if (sebomd_obj%do_sebomd) then
    ! step 2/2 to save emap + entr + ecap forces
    call sebomd_save_forces(2, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
  end if

  ! No energy expression for ifcap == 5 given because only
  !    one step of minimization is allowed with this.

  !  (this seems very weird: we have already done an allreduce on molvir
  !  in ewald_force(); this just collects it on processor 0 (with zeroes
  !  on all slave nodes), then later does an allreduce...)
  if (mytaskid == 0 .and. iamoeba == 0) then
    vir(1) = vir(1)+0.5d0*molvir(1,1)
    vir(2) = vir(2)+0.5d0*molvir(2,2)
    vir(3) = vir(3)+0.5d0*molvir(3,3)
  end if

  if (igb == 0 .and. ipb == 0 .and. iyammp == 0) then
    ener%virvsene = virvsene
    ener%diprms = diprms
    ener%dipiter = dipiter
    ener%dipole_temp = dipole_temp
  end if

  ! Get the noesy volume penalty energy
  pot%noe = 0.d0
  if (iredir(4) /= 0) then
    call timer_start(TIME_NOE)
    call noecalc(x, f, xx, ix)
    call timer_stop(TIME_NOE)
  end if
  ! Do we need a pot%noe here?  mjw TODO

  ! When any form of Generalized Born is active but Poisson-Boltzman is
  ! not, all nonbonded interactions are done in the subroutine egb:
  esurf = 0.d0
  if (igb /= 0 .and. igb /= 10 .and. ipb == 0) then
    call timer_start(TIME_EGB)
    call egb(x, f, rborn, fs, reff, onereff, xx(l15), ix(i04), ix(i06), &
             ix(i08), ix(i10), xx(l190), cut, ntypes, natom, natbel, epol, &
             eelt, evdw, esurf, dvdl, xx(l165), ix(i82), xx(l170), xx(l175), &
             xx(l180), xx(l185), ncopy &
#ifndef LES
             , xx(l2402),xx(l2403),xx(l2404) &
#endif
             )
    pot%vdw  = evdw
    pot%elec = eelt
    pot%gb   = epol
    pot%surf = esurf
    pot%dvdl = dvdl
    call timer_stop(TIME_EGB)
#ifdef MPI
#  ifdef LES
    if (rem == 2) then
      pot%les = pot%les + elesp
    endif
#  endif
#endif
  end if  ! ( igb /= 0 .and. igb /= 10 .and. ipb == 0 )
  ! End handoff of nonbonded computations to subroutine egb

#ifdef RISMSANDER
  ! Force and energy computations by the Reference Interaction Site Model
  if (rismprm%rism == 1) then
    call timer_start(TIME_RISM)
    call rism_force(x, f, erism, irespa, imin)
    pot%rism = erism
    call timer_stop(TIME_RISM)
  endif
#endif

#ifdef APBS
  ! Force computations from the Adaptive Poisson Boltzmann solver
  if (mdin_apbs) then
    if (igb /= 6) then
      write(6, '(a)') '&apbs keyword requires igb=6.'
      call mexit(6,1)
    end if
    call timer_start(TIME_PBFORCE)

    ! Input:  coordinates, radii, charges
    ! Output: updated forces (via apbs_params) and solvation energy
    !         (pol + apolar)
    if (sp_apbs) then
      call apbs_spenergy(natom, x, f, eelt, enpol)
    else
      call apbs_force(natom, x, f, pot%vdw, eelt, enpol)
    end if
    pot%pb   = eelt
    pot%surf = enpol
    call timer_stop(TIME_PBFORCE)

  end if
  ! End APBS force computations
#endif /* APBS */

  if (sebomd_obj%do_sebomd) then
    ! step 1/2 to save eshf + epcshf + ealign + ecsa forces
    call sebomd_save_forces(1, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
  end if
  if (master) then
    !  These parts of the NMR energies are not parallelized, so only
    !  are done on the master node:
    eshf = 0.d0
    epcshf = 0.d0
    ealign = 0.d0
    ecsa = 0.d0
    if (iredir(5) /= 0) call cshf(natom,x,f)
    if (iredir(7) /= 0) call pcshift(natom,x,f)
    if (iredir(9) /= 0) call csa1(natom,x,f)
    if (iredir(8) /= 0) call align1(natom,x,f,xx(lmass))
  end if

  if (sebomd_obj%do_sebomd) then
    ! step 2/2 to save eshf + epcshf + ealign + ecsa forces
    call sebomd_save_forces(2,natom,f,xx(gradsebomd),xx(grad1tmp),xx(grad2tmp))
  end if

  ! additional force due to charge relocation
  if (ifcr /= 0) then
    call cr_calc_force(f)
  end if

  ! MuSiC - GAL17 force field
  call music_force(ipairs, music_vdisp, music_vang, music_vgauss, music_spohr89)

#ifdef MPI
  call timer_barrier( commsander )
  call timer_start(TIME_COLLFRC)
  !     add force, ene, vir, copies from all nodes
  !            also add up newbalance for nonperiodic.

  ! Remember to work on the local instance of the
  ! potential energy array, i.e. pot and NOT the global one,
  ! i.e. ener%pot
  call fdist(f,xx(lfrctmp),pot,vir,newbalance)
  call timer_stop(TIME_COLLFRC)
#endif

  ! ---- at this point, the parallel part of the force calculation is
  !      finished, and the forces have been distributed to their needed
  !      locations.  All forces below here are computed redundantly on
  !      all processors, and added into the force vector.  Hence, below
  !      is the place to put any component of the force calculation that
  !      has not (yet) been parallelized.
  if (sebomd_obj%do_sebomd) then
    ! step 1/2 to save enmr + edssp forces
    call sebomd_save_forces(1, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
  end if

  ! Calculate the NMR restraint energy contributions, if requested.
  ! (Even though this is not parallelized, it needs to be run on all
  ! threads, since this code is needed for weight changes as well as
  ! for NMR restraint energy analysis.  The whole thing could stand a
  ! major re-write....)
  if (nmropt > 0) then
    call nmrcal(x, f, ih(m04), ih(m02), ix(i02), xx(lwinv), enmr, devdis, &
                devang, devtor, devplpt, devpln, devgendis, temp0, tautp, &
                cut, xx(lnmr01), ix(inmr02), xx(l95), 31, 6, rk, tk, pk, cn1, &
                cn2, asol, bsol, xx(l15), numbnd, numang, nptra-nimprp, &
                nimprp, nphb, natom, natom, ntypes, nres, rad, wel, radhb, &
                welhb, rwell, tgtrmsd, temp0les, -1, 'CALC')
  end if
#ifdef MPI
  call mpi_reduce(enoe, pot%noe, 1, MPI_DOUBLE_PRECISION, mpi_sum, 0, &
                  commsander, ierr)
  enoe = pot%noe ! so all processors now have the full enoe value
#else
  pot%noe = enoe
#endif

#ifdef DSSP
  if (idssp > 0) then
    call fdssp(natom, x, f, edssp)
  else
    edssp = 0.d0
  end if
#endif
  if (sebomd_obj%do_sebomd) then
    ! step 2/2 to save enmr + edssp forces
    call sebomd_save_forces(2, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))
  end if

  ! X-ray refinement section: add the target function and gradient
  call cns_xref_run(natom,ih(m04), x,f,ener)

  ! Built-in X-ray target function and gradient
  xray_e = 0.d0
  if( xray_active ) then
#ifdef USE_ISCALE
     if (first) then
        ! set coordinates to current bfactors:
        x(3*natom+1:4*natom) = atom_bfactor(1:natom)
        first = .false.
     else
        ! get current bfactors from the end of the coordinate array:
        atom_bfactor(1:natom) = x(3*natom+1:4*natom)
     endif
     call xray_get_derivative(x,f,xray_e,dB=f(3*natom+1))
#else
     call xray_get_derivative(x,f,xray_e)
#endif
  endif

    ! Calculate the total energy and group the components
#ifndef LES
  if (igb == 0 .and. ipb == 0) then
    pot%vdw_14   = pot%vdw_14   + enb14
    pot%elec_14  = pot%elec_14  + ee14
  endif
#endif
  pot%constraint = pot%constraint + eshf + epcshf + pot%noe + &
                   sum(enmr(1:6)) + ealign + ecsa + pot%emap + xray_e + enfe
#ifdef DSSP
  pot%constraint = pot%constraint + edssp
#endif
  pot%polar = epolar
  pot%tot   = pot%vdw + pot%elec + pot%gb + pot%pb + pot%bond + pot%angle + &
              pot%dihedral + pot%vdw_14 + pot%elec_14 + pot%hbond + &
              pot%constraint + pot%rism + pot%ct
  pot%tot = pot%tot + pot%polar + pot%surf + pot%scf + pot%disp

  !Charmm related
  pot%tot = pot%tot + pot%angle_ub + pot%imp + pot%cmap

  ! MuSiC - GAL17 force field
  pot%tot = pot%tot + music_vdisp + music_vang + music_vgauss + music_spohr89

  ! SEBOMD: apply lambda term to forces and energy between full QM
  ! calculation and full MM calculation.  If lambda .eq. 1.0 (default),
  ! the MM calculations are zeroed.
  if (sebomd_obj%do_sebomd) then

    ! step 3: apply restraint forces to SEBOMD forces
    call sebomd_save_forces(3, natom, f, xx(gradsebomd), xx(grad1tmp), &
                            xx(grad2tmp))

    ! new potential energy
    ! (to have full restraint, we add pot%constraint to the SEBOMD energy)
    ! (see sebomd_save_forces subroutine for explanation)
    pot%tot = (one - sebomd_obj%lambda) * pot%tot + &
               sebomd_obj%lambda * (sebomd_obj%esebomd + pot%constraint)

    ! new forces
    do i = 1, 3*natom
      f(i) = (one - sebomd_obj%lambda)*f(i) + &
             sebomd_obj%lambda * xx(gradsebomd + i - 1)
    end do
  end if

  ! Accelerate MD: calculate the total potential energy weight,
  ! then apply it to all the force elements f=f*fwgt.  Update
  ! the total energy.  Added by Romelia Salomon
  if (iamd .gt. 0) then
    call calculate_amd_total_weights(natom, pot%tot, pot%dihedral, &
                                     pot%amd_boost, f, temp0)
    pot%tot = pot%tot + pot%amd_boost
  end if

  ! scaledMD: scale the total potential and forces by the
  ! factor scaledMD_lambda.  Added by Romelia Salomon
  if (scaledMD .gt. 0) then
    call scaledMD_scale_frc(natom, pot%tot, f)
    pot%tot = pot%tot * scaledMD_lambda
  end if

  ! The handover
  ener%pot = pot
  ener%aveper = aveper
  ener%aveind = aveind
  ener%avetot = avetot

  ! This is now historical; MJW Feb 2010
  !
  !    Here is a summary of how the ene array is used.  For parallel runs,
  !    these values get summed then rebroadcast to all nodes (via
  !    mpi_allreduce).

  !    ene(1):      total energy
  !    ene(2):      van der Waals
  !    ene(3):      electrostatic energy
  !    ene(4):      10-12 (hb) energy, or GB energy when igb.gt.0
  !    ene(5):      bond energy
  !    ene(6):      angle energy
  !    ene(7):      torsion angle energy
  !    ene(8):      1-4 nonbonds
  !    ene(9):      1-4 electrostatics
  !    ene(10):     constraint energy
  !    ene(11-19):  used as scratch, but not needed further below
  !    ene(20):     position constraint energy + cap energy
  !    ene(21):     charging free energy result
  !    ene(22):     noe volume penalty
  !    ene(23):     surface-area dependent energy, or cavity energy
  !    ene(24):     potential energy for a subset of atoms
  !    ene(25):     SCF Energy when doing QMMM
  !    ene(26):     implicit solvation dispersion energy

#ifdef PUPIL_SUPPORT
  ! QM/MM structural considerations are now dealt with.
  ! Add the quantum forces from last QM calculation.
  do iPup = 1, pupparticles
    bs1 = (abs(pupqlist(iPup))-1)*3
    do jPup = 1, 3
      bs2 = bs1 + jPup
      f(bs2) = f(bs2) + qfpup(bs2)
    enddo
  enddo

  ! If there are more that one QM Domain add vdw interaction
  ! among qm particles from different QM Domains
  if (pupnumdomains .gt. 1) then
    call add_vdwqmqm(r_stack(l_puptmp),f,ener,ntypes,ih(m04),ih(m06),ix(i04))
  endif

  ! Deallocate temporary stack
  call free_stack(l_puptmp,routine)

  ! Disconnect QM/MM interactions
  qmmm_nml%ifqnt = .false.
#endif

  ! If freezemol has been set, zero out all of the forces for
  ! the real atoms. (It is no longer necessary to set ibelly.)
  if (ifreeze > 0) then
    do i = 1, 3*natom
      f(i) = 0.d0
    end do
  end if

  ! If a bellymask is being used, set the belly atom forces to zero.
  if (belly) call bellyf(natom,ix(ibellygp),f)

  ! Interface to EVB
#ifdef MPI
#  ifdef LES
  if (nbead > 0) then
    call mpi_allreduce (nrg_all, nrg_bead, nbead, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, commsander, ierr)
    nrg_all(:) = nrg_bead(:)
  end if
  if (ievb /= 0) then
    call evb_ntrfc(x, f, ener, ix, ipairs, vel0_nrg_sum)
  end if
#  else
  if (ievb /= 0) then
    call evb_ntrfc(x, f, ener, xx(lmass), ix, ipairs)
  end if
#  endif /* LES */
#endif /* MPI */

  if (ipimd > 0) then
#ifdef LES
    call pimd_part_spring_force(x, f, real_mass, Epot_spring, Epot_deriv, dvdl)
#  ifdef MPI
    if (ievb /= 0) then
      nrg_frc(3) = vel0_nrg_sum
      nrg_frc(2) = equal_part + Epot_deriv
      nrg_frc(1) = nrg_frc(3) + nrg_frc(2)
      dlnQ_dl = dvdl
    endif
#  endif /* MPI */
#else /* NOT LES below */
   ener = ener/nbead
   f(1:natom*3) = f(1:natom*3)/nbead
   vir(1:3) = vir(1:3)/nbead
   atvir = atvir/nbead
   e14vir = e14vir/nbead
   bnd_vir = bnd_vir/nbead
   call pimd_full_spring_force(x, f, real_mass, Epot_spring, Epot_deriv, dvdl)
#  ifdef MPI
   if (master) then
     call mpi_reduce(ener, totener, state_rec_len, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, 0, commmaster, ierr)
   end if
#  endif /* MPI */
#endif /* LES */
    ! Pass dvdl = dV/dl for TI w.r.t. mass.
    if (itimass > 0) then
      ener%pot%dvdl = dvdl
    end if
  end if

! Ending if of External library
end if

#ifdef MPI
  ! Nudged Elastic Band (NEB) is only supported with MPI
  if (ineb > 0) then

    ! Only the master thread will do NEB
    if (sanderrank .eq. 0) then
      call full_neb_forces(xx(lmass), x, f, ener%pot%tot, ix(itgtfitgp), &
                           ix(itgtrmsgp))
    endif

    ! Now, the master thread will broadcast the NEB forces for the rmsgp
    ! atoms.  All nodes will add this to their current force totals.
    ! If we wanted all nodes to calculate the NEB forces, they would all
    ! need access to the neighbor bead coordinates, which is probably
    ! more expensive than having the master send out the modified forces.
    call mpi_bcast(neb_force, nattgtrms*3, MPI_DOUBLE_PRECISION, 0, &
                   commsander, ierr)

    do i = 1, nattgtrms

      ! Make j3 a pointer into the packed NEB force array.  The
      ! actual atom number for an atom in the rms group can then
      ! be looked up in the ix array.  i3 is set as a pointer
      ! into the force array f.
      j3 = 3*(i - 1)
      j = ix(itgtrmsgp+i-1)
      i3 = 3*(j - 1)
      f(i3+1) = f(i3+1) + neb_force(j3+1)
      f(i3+2) = f(i3+2) + neb_force(j3+2)
      f(i3+3) = f(i3+3) + neb_force(j3+3)
    enddo

    ! CARLOS: what is this doing? ener(27) is neb energy
    ! looks like it reduces entire energy array
    ! needs MUCH better documentation on details (such as 28)
    if (sanderrank .eq. 0) then
      call mpi_reduce(ener, totener, state_rec_len, MPI_DOUBLE_PRECISION, &
                      MPI_SUM, 0, commmaster, ierr)
    end if
  end if
  ! End NEB segment
#endif /*MPI*/

  ! Dump forces in CHARMM format if appropriate and requested
  if (charmm_active) then
    if (do_charmm_dump_gold == 1) then
      call charmm_dump_gold(f, natom, ener)
    endif
  end if

#ifdef MPI
  ! Gather the data in dec (decomposition) arrays of slave processes
  ! back to the master for Thermodynamic Integration.
  if (icfe == 0) then
    if (idecomp == 1 .or. idecomp == 2) then
      call collect_dec2(nres)
    end if
    if (idecomp >= 3) then
      call collect_dec2(npdec*npdec)
    end if
  end if
#endif

  ! End force computations and exit
  call timer_stop(TIME_FORCE)
  call trace_exit('force')
  return

end subroutine force
