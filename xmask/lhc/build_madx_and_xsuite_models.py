import xtrack as xt
import xpart as xp

import xmask as xm

from .errors import install_correct_errors_and_synthesisize_knobs
from .knob_manipulations import rename_coupling_knobs_and_coefficients
from .knob_manipulations import define_octupole_current_knobs
from .knob_manipulations import add_correction_term_to_dipole_correctors


# Assumptions from the madx model and optics
# - Machine energy is stored in madx variable "nrj" (in GeV)
# - The variable "mylhcbeam" is set to 1 in the madx instance holding
#   b1 and b2 (both clockwise), and to 4 in the madx instance holding b4
#   (anti-clockwise).
# - Version of optics is stored in madx variable "ver_lhc_run" for LHC and
#   "ver_hllhc_optics" for HL-LHC
# - Macros are available to save/desable/load orbit bumps, which are called
#   "crossing_save", "crossing_disable", and "crossing_restore".


def build_xsuite_collider(
    sequence_b1, sequence_b2, sequence_b4, beam_config,
    enable_imperfections,
    install_apertures=False,
    enable_knob_synthesis=False,
    pars_for_imperfections={},
    ver_lhc_run=None,
    ver_hllhc_optics=None,
    call_after_last_use=None,):

    """
    Build xsuite collider from madx sequences and optics.

    Parameters
    ----------
    sequence_b1: cpymad.madx.Sequence
        Madx sequence for beam 1 (clockwise)
    sequence_b2: cpymad.madx.Sequence
        Madx sequence for beam 2 (clockwise)
    sequence_b4: cpymad.madx.Sequence
        Madx sequence for beam 4 (anti-clockwise)
    beam_config: dict
        Dictionary with beam configuration (see examples)
    enable_imperfections: bool
        If True, lattice imperfections are installed and corrected
    install_apertures: bool
        If True, apertures are installed
    enable_knob_synthesis: bool
        If True, knobs (linear coupling) are synthesized
    pars_for_imperfections: dict
        Dictionary with parameters for imperfections configuration (see examples)
    ver_lhc_run: str
        Version of LHC optics (None if HL-LHC)
    ver_hllhc_optics: str
        Version of HL-LHC optics (None if LHC)
    call_after_last_use: callable
        Function called after last madx use command


    Returns
    -------
    collider: xtrack.Multiline

    """

    if sequence_b4 is not None:
        assert sequence_b2 is not None, 'If b4 provided, b2 must be provided'

    for beam_n, sequence in zip((1, 2), [sequence_b1, sequence_b2]):
        if sequence is None: continue
        xm.attach_beam_to_sequence(sequence, beam_to_configure=beam_n,
                            beam_configuration=beam_config[sequence.name])
        sequence._madx.globals.nrj = beam_config[sequence.name]['beam_energy_tot']
        # Warm up (seems I need to twiss for mad to load everything)
        sequence._madx.use(sequence.name)
        sequence._madx.twiss()

    # Store energy in nrj
    if sequence_b4 is not None:
        sequence_b4._madx.globals.nrj = beam_config[sequence_b4.name]['beam_energy_tot']

    # Generate beam 4
    if sequence_b4 is not None:
        xm.configure_b4_from_b2(sequence_b4=sequence_b4, sequence_b2=sequence_b2)

    # Save lines for closed orbit reference
    lines_co_ref = xm.save_lines_for_closed_orbit_reference(
        sequence_clockwise=sequence_b1,
        sequence_anticlockwise=sequence_b4)

    lines_to_track = {}
    for sequence_to_track in [sequence_b1, sequence_b4]:

        if sequence_to_track is None:
            continue

        sequence_name = sequence_to_track.name
        mad_track = sequence_to_track._madx

        # Final use
        mad_track.use(sequence_name)

        if call_after_last_use is not None:
            call_after_last_use(mad_track)

        # We work exclusively on the flat machine
        mad_track.input('exec, crossing_disable;')
        mad_track.input('exec, crossing_save;') # In this way crossing_restore keeps the flat machine

        # Install and correct errors
        install_correct_errors_and_synthesisize_knobs(mad_track,
            enable_imperfections=enable_imperfections,
            enable_knob_synthesis=enable_knob_synthesis,
            pars_for_imperfections=pars_for_imperfections,
            ver_lhc_run=ver_lhc_run,
            ver_hllhc_optics=ver_hllhc_optics)

        # Prepare xsuite line
        line = xt.Line.from_madx_sequence(
            mad_track.sequence[sequence_name], apply_madx_errors=True,
            deferred_expressions=True,
            install_apertures=install_apertures,
            replace_in_expr={'bv_aux': 'bvaux_'+sequence_name})
        mad_beam = mad_track.sequence[sequence_name].beam
        line.particle_ref = xp.Particles(p0c = mad_beam.pc*1e9,
            q0 = mad_beam.charge, mass0 = mad_beam.mass*1e9)

        # Prepare coupling and octupole knobs
        rename_coupling_knobs_and_coefficients(line=line,
                                            beamn=int(sequence_name[-1]))
        define_octupole_current_knobs(line=line, beamn=int(sequence_name[-1]))
        lines_to_track[sequence_name] = line

    lines = {}
    lines.update(lines_to_track)
    lines.update(lines_co_ref)
    collider = xt.Multiline(lines=lines)

    if 'lhcb1' in lines_to_track:
        collider['lhcb1_co_ref'].particle_ref = collider['lhcb1'].particle_ref.copy()
    if 'lhcb2' in lines_to_track:
        collider['lhcb2_co_ref'].particle_ref = collider['lhcb2'].particle_ref.copy()

    add_correction_term_to_dipole_correctors(collider)

    return collider