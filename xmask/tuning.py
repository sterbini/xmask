import json
from pathlib import Path

import xtrack as xt

def machine_tuning(line,
        enable_closed_orbit_correction=True,
        enable_linear_coupling_correction=True,
        enable_tune_correction=True,
        enable_chromaticity_correction=True,
        knob_names=None,
        targets=None,
        line_co_ref=None, co_corr_config=None,
        tune_knob_step=1e-8,
        chromaticity_knob_step=1e-4,
        tune_tol=1e-4,
        chromaticity_tol=0.05,
        verbose=False):

    # Correct closed orbit
    if enable_closed_orbit_correction:
        print(f'Correcting closed orbit')
        assert line_co_ref is not None
        assert co_corr_config is not None
        if isinstance(co_corr_config, (str, Path)):
            with open(co_corr_config, 'r') as fid:
                co_corr_config = json.load(fid)

        line.correct_closed_orbit(
                                reference=line_co_ref,
                                correction_config=co_corr_config)

    if enable_linear_coupling_correction:
        assert knob_names is not None
        assert 'c_minus_knob_1' in knob_names
        assert 'c_minus_knob_2' in knob_names
        # Match coupling
        print(f'Matching linear coupling')
        line.match(verbose=verbose,
            vary=[
                xt.Vary(name=knob_names['c_minus_knob_1'],
                        limits=[-0.5e-2, 0.5e-2], step=1e-5),
                xt.Vary(name=knob_names['c_minus_knob_2'],
                        limits=[-0.5e-2, 0.5e-2], step=1e-5)],
            targets=[xt.Target('c_minus', 0, tol=1e-4)])

    # Match tune and chromaticity
    if enable_tune_correction or enable_chromaticity_correction:

        vary = []
        match_targets = []

        if enable_tune_correction:
            assert knob_names is not None
            assert 'q_knob_1' in knob_names
            assert 'q_knob_2' in knob_names
            assert targets is not None
            assert 'qx' in targets
            assert 'qy' in targets

            vary.append(xt.Vary(knob_names['q_knob_1'], step=tune_knob_step))
            vary.append(xt.Vary(knob_names['q_knob_2'], step=tune_knob_step))
            match_targets.append(xt.Target('qx', targets['qx'], tol=tune_tol))
            match_targets.append(xt.Target('qy', targets['qy'], tol=tune_tol))


        if enable_chromaticity_correction:
            assert knob_names is not None
            assert 'dq_knob_1' in knob_names
            assert 'dq_knob_2' in knob_names
            assert targets is not None
            assert 'dqx' in targets
            assert 'dqy' in targets

            vary.append(xt.Vary(knob_names['dq_knob_1'], step=chromaticity_knob_step))
            vary.append(xt.Vary(knob_names['dq_knob_2'], step=chromaticity_knob_step))
            match_targets.append(xt.Target('dqx', targets['dqx'], tol=chromaticity_tol))
            match_targets.append(xt.Target('dqy', targets['dqy'], tol=chromaticity_tol))

        print(f'Matching tune and chromaticity')
        line.match(verbose=verbose, vary=vary, targets=match_targets)
