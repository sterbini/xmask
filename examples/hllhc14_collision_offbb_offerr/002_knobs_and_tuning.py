import xtrack as xt
import xmask as xm

# Load collider anf build trackers
collider = xt.Multiline.from_json('collider_00_from_mad.json')
collider.build_trackers()

# Read knobs and tuning settings from config file
with open('config.yaml','r') as fid:
    config = xm.yaml.load(fid)
conf_knobs_and_tuning = config['config_knobs_and_tuning']

# Set all knobs (crossing angles, dispersion correction, rf, crab cavities,
# experimental magnets, etc.)
for kk, vv in conf_knobs_and_tuning['knob_settings'].items():
    collider.vars[kk] = vv

# Tunings
for line_name in ['lhcb1', 'lhcb2']:

    knob_names = conf_knobs_and_tuning['knob_names'][line_name]

    targets = {
        'qx': conf_knobs_and_tuning['qx'][line_name],
        'qy': conf_knobs_and_tuning['qy'][line_name],
        'dqx': conf_knobs_and_tuning['dqx'][line_name],
        'dqy': conf_knobs_and_tuning['dqy'][line_name],
    }

    xm.machine_tuning(line=collider[line_name],
        enable_closed_orbit_correction=True,
        enable_linear_coupling_correction=True,
        enable_tune_correction=True,
        enable_chromaticity_correction=True,
        knob_names=knob_names,
        targets=targets,
        line_co_ref=collider[line_name+'_co_ref'],
        co_corr_config=conf_knobs_and_tuning['closed_orbit_correction'][line_name],
        tune_knob_step=conf_knobs_and_tuning['tune_knob_step'],
        chromaticity_knob_step=conf_knobs_and_tuning['chromaticity_knob_step'],
        tune_tol=conf_knobs_and_tuning['tune_tol'],
        chromaticity_tol=conf_knobs_and_tuning['chromaticity_tol'],)

collider.to_json('collider_02_tuned_bb_off.json')

#!end-doc-part

# Checks

import numpy as np
for line_name in ['lhcb1', 'lhcb2']:

    assert collider[line_name].particle_ref.q0 == 1
    assert np.isclose(collider[line_name].particle_ref.p0c, 7e12,
                      atol=0, rtol=1e-5)
    assert np.isclose(collider[line_name].particle_ref.mass0, 0.9382720813e9,
                        atol=0, rtol=1e-5)

    tw = collider[line_name].twiss()

    assert np.isclose(tw.qx, 62.31, atol=1e-4, rtol=0)
    assert np.isclose(tw.qy, 60.32, atol=1e-4, rtol=0)
    assert np.isclose(tw.qs, 0.00212, atol=1e-4, rtol=0) # Checks that RF is well set

    assert np.isclose(tw.dqx, 5, atol=0.1, rtol=0)
    assert np.isclose(tw.dqy, 5, atol=0.1, rtol=0)

    assert np.isclose(tw.c_minus, 0, atol=2.5e-4, rtol=0)
    assert np.allclose(tw.zeta, 0, rtol=0, atol=1e-4) # Check RF phase

    # Check separations
    assert np.isclose(tw['ip1', 'x'], 0, rtol=0, atol=7e-8) # sigma is 4e-6
    assert np.isclose(tw['ip1', 'y'], 0, rtol=0, atol=5e-8) # sigma is 4e-6
    assert np.isclose(tw['ip5', 'x'], 0, rtol=0, atol=6e-8) # sigma is 4e-6
    assert np.isclose(tw['ip5', 'y'], 0, rtol=0, atol=5e-8) # sigma is 4e-6

    assert np.isclose(tw['ip2', 'x'],
            -0.14465034135119917e-3 * {'lhcb1': 1, 'lhcb2': 1}[line_name], # set separation
            rtol=0, atol=4e-6)
    assert np.isclose(tw['ip2', 'y'], 0, rtol=0, atol=12e-8)

    assert np.isclose(tw['ip8', 'x'], 0, rtol=0, atol=13e-8)
    assert np.isclose(tw['ip8', 'y'],
            -0.03428196013423046e-3 * {'lhcb1': 1, 'lhcb2': -1}[line_name], # set separation
            rtol=0, atol=5e-8)

    # Check crossing angles
    assert np.isclose(tw['ip1', 'px'],
            250e-6 * {'lhcb1': 1, 'lhcb2': -1}[line_name], rtol=0, atol=0.7e-6)
    assert np.isclose(tw['ip1', 'py'], 0, rtol=0, atol=0.5e-6)
    assert np.isclose(tw['ip5', 'px'], 0, rtol=0, atol=0.5e-6)
    assert np.isclose(tw['ip5', 'py'], 250e-6, rtol=0, atol=0.5e-6)

    assert np.isclose(tw['ip2', 'px'], 0, rtol=0, atol=0.5e-6)
    assert np.isclose(tw['ip2', 'py'], -100e-6 , rtol=0, atol=0.5e-6) # accounts for spectrometer

    assert np.isclose(tw['ip8', 'px'],
            -115e-6* {'lhcb1': 1, 'lhcb2': -1}[line_name], rtol=0, atol=0.5e-6) # accounts for spectrometer
    assert np.isclose(tw['ip8', 'py'], 2e-6, rtol=0, atol=0.5e-6) # small effect from spectrometer (titled)

    assert np.isclose(tw['ip1', 'betx'], 15e-2, rtol=3e-2, atol=0) # beta beating coming from on_disp
    assert np.isclose(tw['ip1', 'bety'], 15e-2, rtol=4e-2, atol=0)
    assert np.isclose(tw['ip5', 'betx'], 15e-2, rtol=2e-2, atol=0)
    assert np.isclose(tw['ip5', 'bety'], 15e-2, rtol=2e-2, atol=0)

    assert np.isclose(tw['ip2', 'betx'], 10., rtol=4e-2, atol=0)
    assert np.isclose(tw['ip2', 'bety'], 10., rtol=4e-2, atol=0)

    assert np.isclose(tw['ip8', 'betx'], 1.5, rtol=3e-2, atol=0)
    assert np.isclose(tw['ip8', 'bety'], 1.5, rtol=2e-2, atol=0)

    # Check crab cavities
    z_crab_test = 1e-2
    phi_crab_1 = ((
        collider[line_name].twiss(method='4d', zeta0=z_crab_test)['ip1', 'x']
      - collider[line_name].twiss(method='4d', zeta0=-z_crab_test)['ip1', 'x'])
      / 2 / z_crab_test)

    phi_crab_5 = ((
        collider[line_name].twiss(method='4d', zeta0=z_crab_test)['ip5', 'y']
      - collider[line_name].twiss(method='4d', zeta0=-z_crab_test)['ip5', 'y'])
      / 2 / z_crab_test)

    assert np.isclose(phi_crab_1, -190e-6 * {'lhcb1': 1, 'lhcb2': -1}[line_name],
                      rtol=1e-2, atol=0)
    assert np.isclose(phi_crab_5, -190e-6, rtol=1e-2, atol=0)

# 
