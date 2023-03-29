import time
#import shutil
#import pickle
import json
import numpy as np

import xtrack as xt
import xfields as xf

def norm(x):
    return np.sqrt(np.sum(np.array(x) ** 2))

atol = 1e-12
rtol = 5.4e-4
strict = False
beta_beat_tol = 5e-2

with open('collider_03_tuned_bb_on.json', 'r') as fid:
    dct = json.load(fid)
    collider = xt.Multiline.from_dict(dct)


for sequence in ['lhcb1', 'lhcb2']:
    print(f'Checking sequence {sequence}')
    
    ltest = collider[sequence]

    if sequence == 'lhcb1':
        lref = xt.Line.from_json('/home/sterbini/2023_03_06_xsuite_test/lhcmask/python_examples/'
                            'hl_lhc_collisions_python_onbb_onerrors/xsuite_lines/line_bb_for_tracking.json')
    elif sequence == 'lhcb2':
        lref = xt.Line.from_json('/home/sterbini/2023_03_06_xsuite_test/lhcmask/python_examples/'
                                'hl_lhc_collisions_python_onbb_onerrors_b2/xsuite_lines/line_bb_for_tracking.json')


    # we swithch off the dipolar elements very weakly powered
    for ll in (lref, ltest):
        for ee in ll.elements:
            if isinstance(ee, xt.Multipole):
                if np.abs(ee.knl[0])<1e-9:
                    ee.knl[0] = 0
                if np.abs(ee.ksl[0])<1e-9:
                    ee.ksl[0] = 0

    # we copy the mss ksl[2]
    # Driving term computed with orbit in lref and with flat machine in ltest

    ks2_max = np.max([np.abs(lref[nn].ksl[2]) for nn in lref.element_names if nn.startswith('mss.')])
            
    for nn in ltest.element_names:
        if nn.startswith('mss.'):
            assert np.isclose(ltest[nn].ksl[2], lref[nn].ksl[2], atol=ks2_max*0.13, rtol=0)
            ltest[nn].ksl[2] = lref[nn].ksl[2]

    # It seems that MADX does not like the coupling of the matching
    if sequence == 'lhcb1':
        lref.vars['cmrskew'] = ltest.vars['c_minus_re_b1']
        lref.vars['cmiskew'] = ltest.vars['c_minus_im_b1']
    else:
        lref.vars['cmrskew'] = ltest.vars['c_minus_re_b2']
        lref.vars['cmiskew'] = ltest.vars['c_minus_im_b2']

    lref.particle_ref = ltest.particle_ref
    lref.build_tracker()
    ltest.build_tracker()

    tw_lref = lref.twiss(matrix_stability_tol=0.05)
    tw_ltest = ltest.twiss()

    lref.unfreeze()
    ltest.unfreeze()


    original_length = ltest.get_length()
    assert (lref.get_length() - original_length) < 1e-6

    # Simplify the two machines
    for ll in (ltest, lref):
        ll._var_management = None
        ll.remove_markers(inplace=True)
        ll.remove_inactive_multipoles(inplace=True)
        ll.remove_zero_length_drifts(inplace=True)
        ll.merge_consecutive_drifts(inplace=True)
        ll.merge_consecutive_multipoles(inplace=True)

        # Remove inactive RFMultipoles and normalize phase
        for ii, ee in enumerate(ll.elements):
            if ee.__class__.__name__ == 'RFMultipole':
                if ( np.max(np.abs(ee.knl)) < 1e-20 and
                        np.max(np.abs(ee.ksl)) < 1e-20 and
                        abs(ee.voltage) < 1e-20):
                    ll.element_names[ii] = None
                # # Normalize phase
                # for kkll, pp in [[ee.knl, ee.pn],
                #                  [ee.ksl, ee.ps]]:
                #     for ii, vv in enumerate(kkll):
                #         if vv < 0:
                #             kkll[ii] = -vv
                #             pp[ii] += 180
                #         if pp[ii]>180:
                #             pp[ii] -= 360

        while None in ll.element_names:
            ll.element_names.remove(None)

    # Check that the two machines are identical
    assert len(ltest) == len(lref)

    assert (ltest.get_length() - original_length) < 1e-6
    assert (lref.get_length() - original_length) < 1e-6

    for ii, (ee_test, ee_ref, nn_test, nn_ref) in enumerate(
        zip(ltest.elements, lref.elements, ltest.element_names, lref.element_names)
    ):
        print(f'Checking element {ii} {nn_test} {nn_ref}         ', end='\r', flush=True)
        assert type(ee_test) == type(ee_ref)

        dtest = ee_test.to_dict()
        dref = ee_ref.to_dict()

        for kk in dtest.keys():

            # Check if they are identical
            if np.isscalar(dref[kk]) and dtest[kk] == dref[kk]:
                continue

            if isinstance(dref[kk], dict):
                if kk=='fieldmap':
                    continue
                if kk=='boost_parameters':
                    continue
                if kk=='Sigmas_0_star':
                    continue

            # Check if the relative error is small
            val_test = dtest[kk]
            val_ref = dref[kk]
            if not np.isscalar(val_ref):
                if len(val_ref) != len(val_test):
                    diff_rel = 100
                else:
                    for iiii, (vvr, vvt) in enumerate(list(zip(val_ref, val_test))):
                        if not np.isclose(vvr, vvt, atol=atol, rtol=rtol):
                            print(f'Issue found on `{kk}[{iiii}]`')
                            diff_rel = 1000
                        else:
                            diff_rel = 0
            else:
                if val_ref > 0:
                    diff_rel = np.abs((val_test - val_ref)/val_ref)
                else:
                    diff_rel = 100

            if diff_rel < rtol:
                continue

            # Check if absolute error is small

            if not np.isscalar(val_ref) and len(val_ref) != len(val_test):
                diff_abs = 1000
            else:
                diff_abs = norm(np.array(val_test) - np.array(val_ref))
            if diff_abs > 0:
                print(f"{nn_test}[{kk}] - test:{dtest[kk]} six:{dref[kk]}")
            if diff_abs < atol:
                continue

            # Exception: drift length (100 um tolerance)
            if not(strict) and isinstance(ee_test, xt.Drift):
                if kk == "length":
                    if diff_abs < 1e-4:
                        continue

            # Exception: multipole lrad is not passed to sixtraxk
            if isinstance(ee_test, xt.Multipole):
                if kk == "length":
                    if np.abs(ee_test.hxl) + np.abs(ee_test.hyl) == 0.0:
                        continue
                if kk == "order" or kk == "inv_factorial_order":
                    # Checked through bal
                    continue
                if kk == 'knl' or kk == 'ksl' or kk == 'bal':
                    if len(val_ref) != len(val_test):
                        lmin = min(len(val_ref), len(val_test))
                        for vv in [val_ref,val_test]:
                            if len(vv)> lmin:
                                for oo in range(lmin, len(vv)):
                                    # we do not care about errors above 10
                                    if vv[oo] != 0 and oo < {'knl':10,
                                                            'ksl':10, 'bal':20}[kk]:
                                        raise ValueError(
                                            'Missing significant multipole strength')

                        val_ref = val_ref[:lmin]
                        val_test = val_test[:lmin]

                    if len(val_ref) == 0 and len(val_test) == 0:
                        continue
                    else:
                        for iiii, (vvr, vvt) in enumerate(list(zip(val_ref, val_test))):                            
                            if not np.isclose(vvr, vvt, atol=atol, rtol=rtol):
                                if iiii==0 and np.isclose(vvr, vvt, atol=1e-9, rtol=0):
                                    diff_rel = 0
                                    pass
                                else:
                                    print(f'Issue found on `{kk}[{iiii}]`')
                                    diff_rel = 1000
                                    break
                            else:
                                diff_rel = 0
                        if diff_rel < rtol:
                            continue

            # Exception: correctors involved in lumi leveling
            passed_corr = False
            for nn_corr in [
                'mcbcv.5l8.b1', 'mcbyv.a4l8.b1', 'mcbxv.3l8',
                'mcbyv.4r8.b1', 'mcbyv.b5r8.b1',
                'mcbyh.b5l2.b1', 'mcbyh.4l2.b1', 'mcbxh.3l2', 'mcbyh.a4r2.b1',
                'mcbch.5r2.b1',
                'mcbcv.5l8.b2', 'mcbyv.a4l8.b2', 'mcbxv.3l8',
                'mcbyv.4r8.b2', 'mcbyv.b5r8.b2',
                'mcbyh.b5l2.b2', 'mcbyh.4l2.b2', 'mcbxh.3l2', 'mcbyh.a4r2.b2',
                'mcbch.5r2.b2', 'mcbch.a5r2.b2', 'mcbyh.4r2.b2', 'mcbxh.3r2',
                'mcbyh.a4l2.b2', 'mcbyh.5l2.b2', 'mcbyv.5r8.b2', 'mcbyv.a4r8.b2',
                'mcbxv.3r8', 'mcbyv.4l8.b2', 'mcbcv.b5l8.b2']:
                if nn_corr in nn_test and kk in ['knl','ksl','bal']:
                    assert len(val_ref)<=2
                    assert len(val_test)<=2
                    diff_rel = norm(val_test-val_ref)/norm(val_ref)
                    passed_corr = (diff_rel < 1e-2)
                    if passed_corr:
                        break
            if not(strict) and  passed_corr:
                continue

            # Exceptions BB4D (separations are recalculated)
            if not(strict) and isinstance(ee_test, xf.BeamBeamBiGaussian2D):
                if kk == "other_beam_shift_x" or kk == "ref_shift_x":
                    diff_abs = ee_test.other_beam_shift_x + ee_test.ref_shift_x - (ee_ref.other_beam_shift_x + ee_ref.ref_shift_x)
                    if diff_abs / np.sqrt(dtest["other_beam_Sigma_11"]) < 0.033: # This is neede to accommodate different leveling routines (1% difference)
                        continue
                if kk == "other_beam_shift_y" or kk == "ref_shift_y":
                    diff_abs = ee_test.other_beam_shift_y + ee_test.ref_shift_y - (ee_ref.other_beam_shift_y + ee_ref.ref_shift_y)
                    if diff_abs / np.sqrt(dtest["other_beam_Sigma_33"]) < 0.03:
                        continue
                
                if kk == "other_beam_Sigma_11":
                    # beta-beating from errors not includede in the ref lattice
                    if diff_rel < beta_beat_tol:
                        continue
                if kk == "other_beam_Sigma_33":
                    if diff_rel < beta_beat_tol:
                        continue
                if kk.startswith('post_subtract'):
                    continue
            if isinstance(ee_test, xf.BeamBeamBiGaussian2D):
                if kk == 'other_beam_q0' or kk == 'other_beam_num_particles':
                    # ambiguity due to old interface
                    if np.abs(ee_test.other_beam_num_particles*ee_test.other_beam_q0 -
                            ee_ref.other_beam_num_particles*ee_ref.other_beam_q0 ) < 1.: # charges
                        continue

            # Exceptions BB6D (angles and separations are recalculated)
            if not(strict) and isinstance(ee_test, xf.BeamBeamBiGaussian3D):
                if kk == "_cos_alpha":
                    if diff_abs < 5e-4:
                        continue
                if kk == "_sin_alpha":
                    if diff_abs < 5e-4:
                        continue
                if kk in  ['_sin_phi','_cos_phi','_tan_phi']:
                    if diff_abs < 1e-6:
                        continue
                if kk == "ref_shift_x" or kk == "other_beam_shift_x":
                    if diff_abs / np.sqrt(dtest["slices_other_beam_Sigma_11_star"][0]) < 0.015:
                        continue
                if kk == "ref_shift_y" or kk == "other_beam_shift_y":
                    if diff_abs / np.sqrt(dtest["slices_other_beam_Sigma_33_star"][0]) < 0.015:
                        continue
                if kk == "ref_shift_zeta":
                    if diff_abs <1e-5:
                        continue
                if kk == "ref_shift_pzeta":
                    if diff_abs <1e-5:
                        continue
                if kk == "ref_shift_px" or kk == 'ref_shift_py':
                    if diff_abs <30e-9:
                        continue
                if kk.startswith('slices_other_beam_Sigma_'):
                    print(kk)
                    ok_sigma = False
                    for nn in '11 12 13 14 22 23 24 33 34 44'.split(' '):
                        # import pdb; pdb.set_trace()
                        if nn in ['12', '34'] and diff_abs < 5e-11:
                            ok_sigma = True
                            break
                        
                        if kk == 'slices_other_beam_Sigma_'+nn+'_star':
                            diff_rel = np.abs((val_test - val_ref)/val_ref)
                            assert len(diff_rel)==1
                            diff_rel = diff_rel[0]
                            if diff_rel < beta_beat_tol:
                                ok_sigma = True
                                break
                    if ok_sigma:
                        continue
                if kk == "other_beam_shift_x" or kk == "ref_shift_x":
                    diff_abs = ee_test.other_beam_shift_x + ee_test.ref_shift_x - (ee_ref.other_beam_shift_x + ee_ref.ref_shift_x)
                    if diff_abs / np.sqrt(dtest["other_beam_Sigma_11"]) < 0.02: # This is neede to accommodate different leveling routines (1% difference)
                        continue
                if kk == "other_beam_shift_y" or kk == "ref_shift_y":
                    diff_abs = ee_test.other_beam_shift_y + ee_test.ref_shift_y - (ee_ref.other_beam_shift_y + ee_ref.ref_shift_y)
                    if diff_abs / np.sqrt(dtest["other_beam_Sigma_33"]) < 0.02:
                        continue
    # bb_ho.l5b1_05
                if kk.startswith('post_subtract'):
                    continue
            if isinstance(ee_test, xt.XYShift):
                if kk in ['dx','dy']:
                    if diff_abs <1e-9:
                        continue
            if isinstance(ee_test, xt.SRotation):
                if kk in ['sin_z', 'cos_z', 'angle']:
                    if diff_abs <1e-9:
                        continue
            # If it got here it means that no condition above is met
            raise ValueError("Too large discrepancy!")
    print(
        f"""
    *********************************************************************************

    The line {sequence} from test sequence and the line from reference are identical!
    
    *********************************************************************************
    """
    )

    time.sleep(1)