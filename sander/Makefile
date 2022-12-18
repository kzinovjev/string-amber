#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
# Makefile for Amber 21 SANDER
#
include ../config.h

# PLUMED-related include file
ifdef PLUMED_INCLUDE_FILE
include $(PLUMED_INCLUDE_FILE)
endif

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA

SHELL=/bin/sh
VPATH=

# Sander Fortran source files are free format
# Add also search paths for .mod files
LOCALFLAGS = $(FREEFORMAT_FLAG) -I../pbsa -I../sqm -I../rism -I../../../include


EMILOBJ= \
    emil.o mdin_emil_dat.o

NFE_OBJECTS = \
    nfe-abmd-hooks.o nfe-bbmd-ctxt.o nfe-bbmd-hooks.o  \
    nfe-colvar.o nfe-colvar-math.o nfe-colvar-type.o nfe-colvar-utils.o \
    nfe-constants.o nfe-cv-ANGLE.o nfe-cv-COM_ANGLE.o \
    nfe-cv-COM_DISTANCE.o nfe-cv-COM_TORSION.o nfe-cv-COS_OF_DIHEDRAL.o \
    nfe-cv-DISTANCE.o nfe-cv-HANDEDNESS.o nfe-cv-LCOD.o \
    nfe-cv-MULTI_RMSD.o nfe-cv-N_OF_BONDS.o nfe-cv-N_OF_STRUCTURES.o \
    nfe-cv-R_OF_GYRATION.o nfe-cv-TORSION.o  \
    nfe-cv-ORIENTATION_ANGLE.o nfe-cv-ORIENTATION_PROJ.o\
    nfe-cv-SPINANGLE.o nfe-cv-TILT.o\
    nfe-cv-QUATERNION0.o nfe-cv-QUATERNION1.o\
    nfe-cv-QUATERNION2.o nfe-cv-QUATERNION3.o\
    nfe-pmd-hooks.o nfe-rmsd.o nfe-sander-hooks.o nfe-sander-proxy.o \
    nfe-stsm-hooks.o nfe-cv-PCA.o nfe-read-pca.o nfe-cv-SIN_OF_DIHEDRAL.o \
    nfe-cv-PAIR_DIHEDRAL.o nfe-cv-PATTERN_DIHEDRAL.o nfe-cv-DF_COM_DISTANCE.o \
    nfe-smd-hooks.o nfe-umbrella.o nfe-utils.o nfe-value.o mt19937.o

EVB_BASE = evb_vars.o evb_mathf.o evb_init.o evb_bcast.o \
    evb_dealloc.o check_input.o evb_io.o evb_keywrd.o \
    exchange_gauss.o exchange_warshel.o react_flux_init.o \
    react_path.o react_flux.o kappa_keck.o \
    schlegel_dg.o schlegel_full.o schlegel_diis.o schlegel_gmres.o \
    schlegel_gbasis.o schlegel_pgbasis.o \
    bnd_strtch.o ang_bend.o dihed_tors.o torinternal.o \
    evb_gfchk.o bmat_inv.o bmat_grad.o \
    unitv.o crossv.o cart2internal.o wdc_bmat.o select_coord.o \
    schlegel_uff.o schlegel_dihed.o schlegel_angle.o \
    schlegel_poh_uff.o schlegel_poh_irc.o schlegel_hcn_cart.o \
    schlegel_hcn_int_irc.o schlegel_hcn_cart_irc.o evb_outoptions.o

EVB_DEBUG = evb_2stdebug.o egap_umb_2stdebug.o dg_debug.o \
    evb_io_debug.o bmat_debug.o bmat_inv_debug.o evb_io_debug2.o \
    morse_anal2num.o xwarshel_anal2num.o xgauss_anal2num.o \
    vdw_anal2num.o rcdbonds_anal2num.o rcbond_anal2num.o \
    dg_grad_anal2num.o

EVBOBJ =  $(EVB_BASE) $(EVB_DEBUG) evb_input.o evb_alloc.o evb_matrix.o \
    morsify.o mod_vdw.o evb_ntrfc.o evb_force.o out_evb.o \
    dbonds_anal2num.o bond_anal2num.o evb_umb.o evb_umb_primitive.o

EVBPIMD =  $(EVB_BASE) $(EVB_DEBUG) evb_input.LES.o evb_alloc.LES.o evb_matrix.LES.o \
    morsify.LES.o mod_vdw.LES.o evb_ntrfc.LES.o evb_force.LES.o \
    evb_pimd_init.LES.o out_evb.LES.o dbonds_anal2num.o bond_anal2num.o \
    evb_umb.LES.o evb_umb_primitive.LES.o

AMOEBAOBJ= \
    amoeba_valence.o amoeba_multipoles.o amoeba_recip.o amoeba_interface.o\
    amoeba_direct.o amoeba_mdin.o amoeba_adjust.o amoeba_self.o\
    amoeba_vdw.o amoeba_induced.o amoeba_runmd.o

POL_GAUSSOBJ= \
    pGM_multipoles.o pGM_recip.o pGM_interface.o\
    pGM_direct.o pGM_mdin.o pGM_self.o\
    pGM_induced.o

APBSOBJ= \
    file_io_dat.APBS.o apbs_vars.APBS.o apbs.APBS.o \
    constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    multisander.o sander.APBS.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.APBS.o runmin.o rdparm.o \
    mdread.APBS.o locmem.o runmd.o relaxmd.o quench.o getcor.o degcnt.o decnvh.o \
    fastwt.o parallel.o shake.o rattle.o ene.o mdwrit.o minrit.o \
    set.o dynlib.APBS.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_bspline.o\
    ew_fft.o nonbond_list.o short_ene.o pol_gauss.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o amopen.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.o thermo_int.o \
    matinv.o assert.o erfcfun.o sinr_t.o md_scheme.o\
    constantph.o constante.o external.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(EMILOBJ) $(AMOEBAOBJ) $(POL_GAUSSOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o sander_lib.o amd.o scaledMD.o barostats.o random.o resamplekin.o boost-gamma-distribution.o les.o music_module.o

QMOBJ= \
    abfqmmm_module.o qm_mm.o qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o \
    qm2_extern_nw_module.o qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_quick_module.o qm2_extern_reaxff_puremd_module.o \
    qm2_extern_genmpi_module.o quick_module.o tcpb_module.o \
    qmmm_adaptive_module.o prn_qmmm_dipole.o

ifdef LIOLIBS
QMOBJ += qm2_extern_lio_module.o
endif

ifdef FBLIBS
QMOBJ += qm2_extern_fb_module.o
endif

QMAPIOBJ = \
    abfqmmm_module.o qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o qm2_extern_nw_module.o \
    qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_quick_module.o qm2_extern_genmpi_module.o quick_module.o tcpb_module.o \
    qm2_extern_reaxff_puremd_module.o qmmm_adaptive_module.o prn_qmmm_dipole.o

ifdef LIOLIBS
QMAPIOBJ += qm2_extern_lio_module.o
endif

ifdef FBLIBS
QMAPIOBJ += qm2_extern_fb_module.o
endif

SEBOMDOBJ = \
    sebomd_module.o sebomd_arrays.o

SEBOMDLIB = ../sebomd/sebomd.a

MMOBJ=  file_io_dat.o constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    multisander.o sander.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o rdparm.o sander_lib.o \
    mdread.o locmem.o runmd.o relaxmd.o quench.o getcor.o degcnt.o decnvh.o \
    fastwt.o parallel.o shake.o rattle.o ene.o mdwrit.o minrit.o \
    set.o dynlib.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.o ew_box.o ew_bspline.o\
    ew_fft.o nonbond_list.o short_ene.o pol_gauss.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o amopen.o \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.o thermo_int.o \
    matinv.o assert.o erfcfun.o sinr_t.o md_scheme.o \
    constantph.o constante.o external.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o\
    $(EMILOBJ) $(AMOEBAOBJ) $(POL_GAUSSOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o amd.o scaledMD.o barostats.o random.o resamplekin.o boost-gamma-distribution.o   les.o music_module.o \
    CV_module.o multiCV_module.o splines_utilities.o string_utilities.o CV_utilities.o string_method.o

APIOBJ=  file_io_dat.o constants.o state.o memory_module.o stack.o \
    nose_hoover.o nose_hoover_init.o \
    trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
    icosasurf.o egb.o remd.o findmask.o \
    relax_mat.o nmr.o multitmd.o \
    trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o rdparm.API.o cpparm.o sander_lib.o \
    mdread.o locmem.o runmd.o relaxmd.o quench.o getcor.o degcnt.o decnvh.o \
    fastwt.API.o parallel.o shake.o rattle.o ene.o mdwrit.o minrit.o \
    set.API.o dynlib.API.o mdfil.o nmlsrc.o \
    ew_force.o ew_setup.API.o ew_box.API.o ew_bspline.o\
    ew_fft.o nonbond_list.API.o short_ene.o pol_gauss.o ew_recip.o pcshift.o align.o csa.o \
    rfree.o rgroup.o amopen.o \
    debug.API.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    mexit.o new_time.API.o extra_pts.API.o thermo_int.o les.o \
    matinv.o assert.o erfcfun.o sinr_t.o md_scheme.o \
    constantph.o constante.o external.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o\
    $(EMILOBJ) $(AMOEBAOBJ) $(POL_GAUSSOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.API.o crg_reloc.o amd.o scaledMD.o barostats.o interface.o random.o resamplekin.o boost-gamma-distribution.o qm_mm.API.o prmtop_type.o music_module.o

LESAPIOBJ= memory_module.o stack.o file_io_dat.o les.o \
    constants.o state.o nose_hoover.o nose_hoover_init_LES.o \
    trace.o rmsgrd.o rand2.o lmodC.o xminC.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o relax_mat.o nmr.LES.o multitmd.o \
    cshf.o nmrcal.o mtmdcall.o pearsn.o cpparm.LES.o \
    printe.o runmin.o nonbond_list.LESAPI.o force.LES.o rdparm.LESAPI.o \
    mdread.LES.o locmem.LES.o runmd.LES.o relaxmd.o quench.o degcnt.LES.o decnvh.o \
    fastwt.API.o parallel.LES.o shake.o rattle.o ene.LES.o \
    dynlib.LESAPI.o nmlsrc.o \
    ew_force.LES.o ew_setup.LESAPI.o ew_bspline.o\
    ew_fft.o short_ene.LES.o pol_gauss.o ew_recip.LES.o pcshift.o \
    align.o csa.o rfree.o rgroup.o \
    amopen.o sander_lib.o random.o resamplekin.o boost-gamma-distribution.o \
    debug.API.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    mexit.o new_time.API.o extra_pts.LESAPI.o thermo_int.o matinv.o assert.o \
    erfcfun.o sinr_t.o md_scheme.LES.o \
    constantph.o constante.o external.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(AMOEBAOBJ) $(POL_GAUSSOBJ) $(EMILOBJ) $(SANDER_RISM_INTERFACE) \
    $(SEBOMDOBJ) \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.API.o crg_reloc.o barostats.o prmtop_type.o \
    AmberNetcdf.o bintraj.LES.o binrestart.LES.o getcor.LES.o mdfil.LES.o mdwrit.LES.o \
    lmod.LES.o trajene.o ew_box.LESAPI.o minrit.LES.o set.LESAPI.o qm_mm.LESAPI.o \
    amd.o scaledMD.o interface.LESAPI.o music_module.o

FULLPIMDOBJ = pimd_vars.o pimd_force.o pimd_init.o cmd_vars.o cmd_matrix.o

PARTPIMDOBJ = pimd_vars.o pimd_force.LES.o pimd_init.LES.o cmd_vars.o cmd_matrix.LES.o

PUPILOBJ= file_io_dat.o constants.o state.o memory_module.o stack.o \
	nose_hoover.o nose_hoover_init.o \
	trace.o rmsgrd.o lmod.o rand2.o lmodC.o xminC.o decomp.o \
	pupildata.o icosasurf.o egb.o findmask.o \
	relax_mat.o nmr.o multitmd.o \
	multisander.o sander.PUPIL.o trajene.o cshf.o nmrcal.o mtmdcall.o pearsn.o \
	printe.PUPIL.o runmin.o force.PUPIL.o rdparm.o  \
	mdread.PUPIL.o locmem.o runmd.o relaxmd.o quench.o getcor.o degcnt.o decnvh.o \
	fastwt.o parallel.o shake.o rattle.o ene.o mdwrit.o minrit.o \
	set.o dynlib.PUPIL.o mdfil.PUPIL.o nmlsrc.o \
	ew_force.o ew_setup.o ew_box.o ew_bspline.o\
	ew_fft.o nonbond_list.o short_ene.o pol_gauss.o ew_recip.o pcshift.o align.o csa.o \
	rfree.o rgroup.o amopen.o sander_lib.o \
	debug.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
	mexit.PUPIL.o new_time.o extra_pts.o thermo_int.o \
	matinv.o assert.o erfcfun.o sinr_t.o md_scheme.o \
	constantph.o constante.o external.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
	$(EMILOBJ) $(AMOEBAOBJ) $(POL_GAUSSOBJ) $(SANDER_RISM_INTERFACE) AmberNetcdf.o bintraj.o binrestart.o \
	spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o \
	putvalues.o fixport.o  dssp.o charmm.o crg_reloc.o amd.o scaledMD.o \
	barostats.o random.o resamplekin.o boost-gamma-distribution.o  music_module.o

LSCIVROBJ = lscivr_vars.o lsc_init.o lsc_xp.o

LESOBJ= memory_module.o stack.o file_io_dat.o \
    constants.o state.o nose_hoover.o nose_hoover_init_LES.o \
    trace.o rmsgrd.o rand2.o lmodC.o xminC.o decomp.o icosasurf.o \
    egb.LES.o remd.LES.o findmask.o \
    relax_mat.o nmr.LES.o multitmd.o \
    multisander.LES.o sander.LES.o  cshf.o nmrcal.o mtmdcall.o pearsn.o \
    printe.o runmin.o nonbond_list.LES.o force.LES.o rdparm.LES.o \
    mdread.LES.o locmem.LES.o runmd.LES.o relaxmd.o quench.o degcnt.LES.o decnvh.o \
    fastwt.o parallel.LES.o shake.o rattle.o ene.LES.o \
    dynlib.LES.o nmlsrc.o \
    ew_force.LES.o ew_setup.LES.o ew_bspline.o\
    ew_fft.o short_ene.LES.o pol_gauss.o ew_recip.LES.o pcshift.o \
    align.o csa.o rfree.o rgroup.o \
    amopen.o sander_lib.o random.o resamplekin.o boost-gamma-distribution.o  \
    debug.o ew_recip_reg.o ew_handle_dips.o ew_bndpol.o ew_dipole_recip.o \
    mexit.o new_time.o extra_pts.LES.o thermo_int.o matinv.o assert.o \
    erfcfun.o sinr_t.o md_scheme.LES.o \
    constantph.o constante.o external.o prn_dipoles.o ips.o sglds.o iwrap2.o emap.o \
    $(AMOEBAOBJ) $(POL_GAUSSOBJ) $(EMILOBJ) $(SANDER_RISM_INTERFACE) \
    $(SEBOMDOBJ) les.o \
    spatial_recip.o spatial_fft.o parms.o softcore.o mbar.o linear_response.o molecule.o xref.o dssp.o \
    charmm.o crg_reloc.o barostats.o \
    qm_ewald.o qmmm_vsolv.o \
    qm2_extern_util_module.o qm2_extern_module.o \
    qm2_extern_adf_module.o qm2_read_adf_results.o KFReader.o \
    ArrayList.o qm2_extern_gms_module.o qm2_extern_tc_module.o \
    qm2_extern_gau_module.o qm2_extern_orc_module.o qm2_extern_nw_module.o \
    qm2_extern_qc_module.o qm2_extern_mrcc_module.o \
    qm2_extern_quick_module.o qm2_extern_reaxff_puremd_module.o \
    qm2_extern_genmpi_module.o quick_module.o tcpb_module.o \
    qmmm_adaptive_module.o abfqmmm_module.o prn_qmmm_dipole.o \
    AmberNetcdf.o bintraj.LES.o binrestart.LES.o getcor.LES.o mdfil.LES.o mdwrit.LES.o \
    lmod.LES.o trajene.o ew_box.LES.o minrit.LES.o set.LES.o qm_mm.LES.o amd.o scaledMD.o \
    music_module.o \
    CV_module.o multiCV_module.o splines_utilities.o string_utilities.o CV_utilities.o string_method.o

ifdef LIOLIBS
LESOBJ += qm2_extern_lio_module.o
endif

ifdef FBLIBS
LESOBJ += qm2_extern_fb_module.o
endif

XRAY_OBJS = \
   xray_fftpack.o xray_FFT.o xray_fourier.o xray_globals.o xray_interface.o \
   xray_real_space.o xray_reciprocal_space.o xray_utils.o

.PHONY: configured configured_serial configured_parallel serial parallel sander

install: $(INSTALLTYPE)

SERIALPROGS = $(BINDIR)/sander$(SFX) $(BINDIR)/sander.LES$(SFX) $(BINDIR)/ambmask$(SFX)

.NOTPARALLEL: $(SERIALPROGS)

serial: configured_serial $(SERIALPROGS) $(BUILD_SANDERAPI) make_test_links

serial_notest: $(SERIALPROGS) $(BUILD_SANDERAPI)

sander: configured_serial $(BINDIR)/sander$(SFX)

build_sanderapi: $(LIBDIR)/libsander$(SHARED_SUFFIX) \
                 $(LIBDIR)/libsanderles$(SHARED_SUFFIX) $(INCDIR)/sander.h

make_test_links:
	cd ../../../test/amd             && ln -s -f ../../bin/sander .
	cd ../../../test/chamber         && ln -s -f ../../bin/sander .
	cd ../../../test/nmropt          && ln -s -f ../../bin/sander .
	cd ../../../test/qmmm2           && ln -s -f ../../bin/sander .
	cd ../../../test/qmmm_DFTB       && ln -s -f ../../bin/sander .
	cd ../../../test/abfqmmm         && ln -s -f ../../bin/sander .
	cd ../../../test/rism3d          && ln -s -f ../../bin/sander .
	cd ../../../test/sander_pbsa_frc && ln -s -f ../../bin/sander .
	cd ../../../test/rism3d          && ln -s -f ../../bin/sander .

skip_sanderapi:
	@echo ""
	@echo "Skipping SANDER API build."
	@echo ""

CUDAPROGS = $(BINDIR)/sander.quick.cuda$(SFX) $(BINDIR)/sander.LES.quick.cuda$(SFX)

cuda_serial: configured_cuda $(CUDAPROGS)

CUDAPARALLELPROGS = $(BINDIR)/sander.quick.cuda.MPI$(SFX) $(BINDIR)/sander.LES.quick.cuda.MPI$(SFX)

cuda_parallel: configured_cuda configured_parallel $(CUDAPARALLELPROGS)

PARALLELPROGS = $(BINDIR)/sander.MPI$(SFX) $(BINDIR)/sander.LES.MPI$(SFX)

.NOTPARALLEL: $(PARALLELPROGS)

parallel: configured_parallel $(PARALLELPROGS)

all_serial_programs: $(SERIALPROGS)
	$(MAKE) $(BINDIR)/sander.PUPIL$(SFX)
	$(MAKE) $(BINDIR)/sander.APBS$(SFX)

#---------------------------------------------------------------------------
$(BINDIR)/sander$(SFX): libsqm $(MMOBJ) $(QMOBJ) \
              $(FULLPIMDOBJ) $(EMIL) $(SEBOMDOBJ) sebomd \
              $(LSCIVROBJ) force.o syslib ../lib/nxtsec.o netlib \
              $(NFE_OBJECTS) $(XRAY_OBJS) libpbsa libquick libtcpb librism $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(MMOBJ) $(QMOBJ) $(SEBOMDOBJ) \
	      $(FULLPIMDOBJ) $(LSCIVROBJ) force.o -L$(LIBDIR) -lsqm \
	      $(NFE_OBJECTS) $(XRAY_OBJS) \
	      -lFpbsa ../lib/nxtsec.o \
	      $(SEBOMDLIB) $(FBLIBS) \
	      ../lib/sys.a $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
	      $(LDFLAGS) $(AMBERLDFLAGS)  $(LOCALFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.quick.cuda$(SFX): configured_cuda libsqm $(MMOBJ) $(QMOBJ) \
              $(FULLPIMDOBJ) $(EMIL) $(SEBOMDOBJ) sebomd \
              $(LSCIVROBJ) force.o syslib ../lib/nxtsec.o netlib \
              $(NFE_OBJECTS) $(XRAY_OBJS) libpbsa libquick libtcpb librism $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(MMOBJ) $(QMOBJ) $(SEBOMDOBJ) \
	      $(FULLPIMDOBJ) $(LSCIVROBJ) force.o -L$(LIBDIR) -lsqm \
	      $(NFE_OBJECTS) $(XRAY_OBJS) \
	      -lFpbsa ../lib/nxtsec.o \
	      $(SEBOMDLIB) $(FBLIBS) \
	      ../lib/sys.a $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
	      $(LDFLAGS) $(AMBERLDFLAGS)  $(LOCALFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.MPI$(SFX): configured_parallel libsqm $(MMOBJ) $(QMOBJ) \
              $(FULLPIMDOBJ) $(EMIL) $(LSCIVROBJ) $(EVBOBJ) force.o \
              $(SEBOMDOBJ) sebomd $(XRAY_OBJS) \
              syslib ../lib/nxtsec.o netlib libpbsa libquick libtcpb librism_mpi $(NFE_OBJECTS) \
              $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(MMOBJ) $(QMOBJ) $(SEBOMDOBJ) \
              $(FULLPIMDOBJ) $(LSCIVROBJ) $(EVBOBJ) force.o \
	      -L$(LIBDIR) -lsqm -lFpbsa \
	      $(SEBOMDLIB) $(FBLIBS) $(XRAY_OBJS) \
	      ../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
	      $(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.quick.cuda.MPI$(SFX): configured_cuda configured_parallel libsqm $(MMOBJ) $(QMOBJ) \
              $(FULLPIMDOBJ) $(EMIL) $(LSCIVROBJ) $(EVBOBJ) force.o \
              $(SEBOMDOBJ) sebomd $(XRAY_OBJS) \
              syslib ../lib/nxtsec.o netlib libpbsa libquick libtcpb librism_mpi $(NFE_OBJECTS) \
              $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(MMOBJ) $(QMOBJ) $(SEBOMDOBJ) \
              $(FULLPIMDOBJ) $(LSCIVROBJ) $(EVBOBJ) force.o \
	      -L$(LIBDIR) -lsqm -lFpbsa \
	      $(SEBOMDLIB) $(FBLIBS) $(XRAY_OBJS) \
	      ../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
	      $(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.PUPIL$(SFX): configured_serial libsqm $(PUPILOBJ) $(QMOBJ) \
              $(FULLPIMDOBJ) $(EMIL) $(LSCIVROBJ) syslib $(XRAY_OBJS) \
              $(SEBOMDOBJ) sebomd \
              ../lib/nxtsec.o netlib libpbsa libquick libtcpb librism $(NFE_OBJECTS) \
              $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(PUPILOBJ) $(QMOBJ) $(SEBOMDOBJ) $(FULLPIMDOBJ) \
		$(LSCIVROBJ) -L$(LIBDIR) -lsqm -lFpbsa \
	        $(SEBOMDLIB) $(FBLIBS) $(XRAY_OBJS) \
		../lib/nxtsec.o $(EMILLIB) ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	    	$(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) \
		$(PUPILLIBS) $(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)


#---------------------------------------------------------------------------
$(BINDIR)/sander.LES$(SFX): libsqm $(LESOBJ) $(PARTPIMDOBJ) \
              syslib ../lib/nxtsec.o netlib \
              $(LSCIVROBJ) $(NFE_OBJECTS) $(XRAY_OBJS) \
              libpbsa libquick libtcpb librism $(EMIL) sebomd $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(LESOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		$(XRAY_OBJS) -L$(LIBDIR) -lsqm -lFpbsa \
	        $(SEBOMDLIB) $(FBLIBS) \
		../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	      	$(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
		$(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)


#---------------------------------------------------------------------------
$(BINDIR)/sander.LES.quick.cuda$(SFX): libsqm $(LESOBJ) $(PARTPIMDOBJ) \
              syslib ../lib/nxtsec.o netlib \
              $(LSCIVROBJ) $(NFE_OBJECTS) $(XRAY_OBJS) \
              libpbsa libquick libtcpb librism $(EMIL) sebomd $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(LESOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		$(XRAY_OBJS) -L$(LIBDIR) -lsqm -lFpbsa \
	        $(SEBOMDLIB) $(FBLIBS) \
		../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	      	$(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
		$(LDFLAGS) $(AMBERLDFLAGS)  $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.LES.MPI$(SFX): configured_parallel libsqm $(LESOBJ) $(EVBPIMD) \
              $(PARTPIMDOBJ) syslib ../lib/nxtsec.o \
              netlib libpbsa libquick libtcpb librism_mpi $(EMIL) $(XRAY_OBJS) \
              $(LSCIVROBJ) $(NFE_OBJECTS) sebomd $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(LESOBJ) $(EVBPIMD) \
              $(PARTPIMDOBJ) $(LSCIVROBJ) $(XRAY_OBJS) \
              -L$(LIBDIR) -lsqm -lFpbsa $(SEBOMDLIB) $(FBLIBS) \
              ../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
              $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
              $(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.LES.quick.cuda.MPI$(SFX): configured_cuda configured_parallel libsqm $(LESOBJ) $(EVBPIMD) \
              $(PARTPIMDOBJ) syslib ../lib/nxtsec.o \
              netlib libpbsa libquick libtcpb librism_mpi $(EMIL) $(XRAY_OBJS) \
              $(LSCIVROBJ) $(NFE_OBJECTS) sebomd $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $(LESOBJ) $(EVBPIMD) \
             $(PARTPIMDOBJ) $(LSCIVROBJ) $(XRAY_OBJS) \
             -L$(LIBDIR) -lsqm -lFpbsa $(SEBOMDLIB) $(FBLIBS) \
             ../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
             $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) $(EMILLIB) \
             $(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/sander.APBS$(SFX): configured_serial libsqm $(APBSOBJ) $(QMOBJ) \
              $(LSCIVROBJ) $(FULLPIMDOBJ) $(NFE_OBJECTS) \
              force.APBS.o syslib libquick libtcpb librism $(EMIL) \
              $(SEBOMDOBJ) sebomd $(XRAY_OBJS) \
              ../lib/nxtsec.o netlib $(PLUMED_DEPENDENCIES)
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(APBSOBJ) $(QMOBJ) $(SEBOMDOBJ) \
              $(FULLPIMDOBJ) $(LSCIVROBJ) force.APBS.o \
              -L$(APBS_LIBDIR) $(APBS_LIBS) \
              -L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	        $(SEBOMDLIB) $(FBLIBS) $(XRAY_OBJS) \
              ../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) $(FFTW3) $(SANDER_FLIBSF) \
              $(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS) $(PLUMED_LOAD)

#---------------------------------------------------------------------------
$(BINDIR)/ambmask$(SFX): 	ambmask.o findmask.o amopen.o parms.o \
              memory_module.o mexit.o ../lib/nxtsec.o
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ ambmask.o findmask.o \
		amopen.o mexit.o ../lib/nxtsec.o constants.o memory_module.o parms.o \
		$(LDFLAGS) $(AMBERLDFLAGS) $(LIOLIBS)


#-----------LIBS

$(LIBDIR)/libxray.a: $(XRAY_OBJS)
	@echo "[SANDER]  AR $@"
	$(AR) $@ $(XRAY_OBJS)
	$(RANLIB) $@

$(LIBDIR)/libsander$(SHARED_SUFFIX): libsqm $(APIOBJ) \
              $(QMAPIOBJ) $(FULLPIMDOBJ) $(EMIL) \
              $(SEBOMDOBJ) sebomd $(INCDIR)/sander_api.mod \
              $(LSCIVROBJ) force.o syslib $(PLUMED_DEPENDENCIES) \
              ../lib/nxtsec.o netlib \
              $(NFE_OBJECTS) $(XRAY_OBJS) libpbsa libquick libtcpb librism fftw3_dummy.o
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(MAKE_SHARED) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@  $(APIOBJ) \
		  $(QMAPIOBJ) $(SEBOMDOBJ) $(FULLPIMDOBJ) $(LSCIVROBJ) force.o \
		  -L$(LIBDIR) -lsqm \
	      $(NFE_OBJECTS) $(XRAY_OBJS) \
	      -lFpbsa ../lib/nxtsec.o $(EMILLIB) \
	      $(SEBOMDLIB) $(FBLIBS) \
	      ../lib/sys.a $(NETCDFLIBF) \
	      $(FLIBS_RISMSANDER) fftw3_dummy.o $(SANDER_FLIBSF) \
	      $(LIOLIBS) $(PLUMED_LOAD)

$(LIBDIR)/libsanderles$(SHARED_SUFFIX): libsqm $(LESAPIOBJ) \
		$(PARTPIMDOBJ) $(QMAPIOBJ) syslib \
		../lib/nxtsec.o netlib $(INCDIR)/sanderles_api.mod \
		$(LSCIVROBJ) $(NFE_OBJECTS) $(XRAY_OBJS) \
		libpbsa libquick libtcpb librism $(EMIL) sebomd $(PLUMED_DEPENDENCIES) fftw3_dummy.o
	@echo "[SANDER]  FC $@"
	$(VB)$(FC) $(MAKE_SHARED) $(SANDER_FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) \
       -o $@ $(LESAPIOBJ) $(QMAPIOBJ) $(PARTPIMDOBJ) $(LSCIVROBJ) \
		$(XRAY_OBJS) -L$(LIBDIR) -lsqm -lFpbsa $(EMILLIB) \
	        $(SEBOMDLIB) $(FBLIBS) $(PLUMED_LOAD) \
		../lib/nxtsec.o ../lib/sys.a $(NFE_OBJECTS) $(NETCDFLIBF) \
	      	$(FLIBS_RISMSANDER) fftw3_dummy.o $(SANDER_FLIBSF) $(LIOLIBS)

$(INCDIR)/sander_api.mod: interface.o
	$(VB)/bin/cp sander_api.mod $(INCDIR)

$(INCDIR)/sanderles_api.mod: interface.o
	$(VB)/bin/cp sanderles_api.mod $(INCDIR)

$(INCDIR)/sander.h: sander.h
	$(VB)/bin/cp sander.h $(INCDIR)

EMIL: memory_module.o
	$(MAKE) -C ../emil install

syslib:
	cd ../lib && $(MAKE) sys.a

netlib:
	cd ../blas && $(MAKE) $(BLAS)
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../arpack && $(MAKE) install

sebomd: $(SEBOMDOBJ)
	cd ../sebomd && $(MAKE) sebomd.a

../lib/nxtsec.o:  ../lib/nxtsec.F
	cd ../lib && $(MAKE) nxtsec.o

fftw3_dummy.o:  fftw3_dummy.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c fftw3_dummy.c

pmemd_clib:
	cd ../pmemd/src && $(CC) -c $(CPPFLAGS) $(CNOOPTFLAGS) \
		$(CFLAGS) $(AMBERCFLAGS) pmemd_clib.c

../lib/mexit.o:  ../lib/mexit.F
	cd ../lib && $(MAKE) mexit.o

ifndef PLUMED_INCLUDE_FILE
Plumed.o: Plumed.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
	        -o Plumed.o Plumed.c
endif

rand2.o:  ../sff/rand2.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o rand2.o ../sff/rand2.c

lmodC.o:  ../sff/lmodC.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o lmodC.o ../sff/lmodC.c

xminC.o:  ../sff/xminC.c
	@echo "[SANDER]  CC $<"
	$(VB)$(CC) -c -DSQM $(CNOOPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) \
		-o xminC.o ../sff/xminC.c

libsqm:
	cd ../sqm && $(MAKE) $@

libquick:
	if [ "$(QUICK)" = 'yes' ]; then \
		cd ../quick && $(MAKE) install ; \
	fi

libtcpb:
	if [ "$(TCPB)" = 'yes' ]; then \
		cd ../tcpb-cpp && $(MAKE) install ; \
	fi

librism librism_mpi:
	if [ -n "$(RISMSANDER)" ]; then \
		cd ../rism && $(MAKE) $@; \
	fi

libpbsa:
	cd ../pbsa && $(MAKE) -j1 libFpbsa.a
	@/bin/mv ../pbsa/libFpbsa.a $(LIBDIR)

clean:
	/bin/rm -f sander.PUPIL$(SFX)
	/bin/rm -f sander.APBS$(SFX)
	/bin/rm -f *.o
	/bin/rm -f *.mod
	/bin/rm -f *.d
	/bin/rm -f *nbflag
	/bin/rm -f *__genmod.f90  # intel compiler-generated interface modules

uninstall.serial:
	/bin/rm -f $(BINDIR)/mopac.sh
	-for file in $(SERIALPROGS) ; do \
	  /bin/rm -f $${file} ; \
	done
	/bin/rm -f $(LIBDIR)/libsander.so $(LIBDIR)/libsander.dylib
	/bin/rm -f $(LIBDIR)/libsanderles.so $(LIBDIR)/libsanderles.dylib
	/bin/rm -f $(INCDIR)/sander_api.mod $(INCDIR)/sander_api.MOD
	/bin/rm -f $(INCDIR)/sanderles_api.mod $(INCDIR)/sanderles_api.MOD
	/bin/rm -f $(INCDIR)/sander.h

uninstall.cuda:
	-for file in $(CUDAPROGS) ; do \
	  /bin/rm -f $${file} ; \
	done

uninstall.cuda_parallel:
	-for file in $(CUDAPARALLELPROGS) ; do \
	  /bin/rm -f $${file} ; \
	done

uninstall.parallel:
	-for file in $(PARALLELPROGS) ; do \
	  /bin/rm -f $${file} ; \
	done

uninstall: uninstall.serial uninstall.cuda uninstall.cuda_parallel uninstall.parallel

depend::
	./makedepend > depend

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) and do NOT specify -mpi." ;\
		exit 2 ;\
	fi ;\
	)

configured_cuda: configured
	@(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is not of type cuda !" ;\
			echo "  Rerun $(CONFIG_COMMAND) and specify -cuda." ;\
			exit 2 ;\
		} ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
			exit 2 ;\
		} ;\
	)
#----------HEADER DEPENDENCIES (constructed from ./makedepend > depend)
include depend
# DO NOT DELETE
