import xtrack as xt
import xpart as xp

import pymaskmx as pm

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
    enable_knob_synthesis,
    pars_for_imperfections,
    ver_lhc_run,
    ver_hllhc_optics):

    pm.attach_beam_to_sequence(sequence_b1, beam_to_configure=1,
                                beam_configuration=beam_config['lhcb1'])
    pm.attach_beam_to_sequence(sequence_b2, beam_to_configure=2,
                                beam_configuration=beam_config['lhcb2'])

    # Warm up (seems I need to twiss for mad to load everything)
    for seq in [sequence_b1, sequence_b2]:
        mm = seq._madx
        mm.use(seq.name)
        mm.twiss()

    # Generate beam 4
    pm.configure_b4_from_b2(
        sequence_b4=sequence_b4,
        sequence_b2=sequence_b2)

    # Save lines for closed orbit reference
    lines_co_ref = pm.save_lines_for_closed_orbit_reference(
        sequence_clockwise=sequence_b1,
        sequence_anticlockwise=sequence_b4)

    lines_to_track = {}
    for sequence_to_track in [sequence_b1, sequence_b4]:

        sequence_name = sequence_to_track.name
        mad_track = sequence_to_track._madx

        # Final use
        mad_track.use(sequence_name)

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
            replace_in_expr={'bv_aux': 'bvaux_'+sequence_name})
        mad_beam = mad_track.sequence[sequence_name].beam
        line.particle_ref = xp.Particles(p0c = mad_beam.pc*1e9,
            q0 = mad_beam.charge, mass0 = mad_beam.mass*1e9)

        # Prepare coupling and octupole knobs
        rename_coupling_knobs_and_coefficients(line=line,
                                            beamn=int(sequence_name[-1]))
        define_octupole_current_knobs(line=line, beamn=int(sequence_name[-1]))
        lines_to_track[sequence_name] = line


    collider = xt.Multiline(
        lines={
            'lhcb1': lines_to_track['lhcb1'],
            'lhcb2': lines_to_track['lhcb2'],
            'lhcb1_co_ref': lines_co_ref['lhcb1_co_ref'],
            'lhcb2_co_ref': lines_co_ref['lhcb2_co_ref'],
        })

    collider['lhcb1_co_ref'].particle_ref = collider['lhcb1'].particle_ref.copy()
    collider['lhcb2_co_ref'].particle_ref = collider['lhcb2'].particle_ref.copy()

    add_correction_term_to_dipole_correctors(collider)

    return collider