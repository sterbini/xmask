import numpy as np
import xmask.lhc as xmlhc

def check_optics_orbit_etc(collider, line_names):

    for line_name in line_names:

        assert collider[line_name].particle_ref.q0 == 1
        assert np.isclose(collider[line_name].particle_ref.p0c, 7e12,
                        atol=0, rtol=1e-5)
        assert np.isclose(collider[line_name].particle_ref.mass0, 0.9382720813e9,
                            atol=0, rtol=1e-5)

        tw = collider[line_name].twiss()

        if line_name == 'lhcb1':
            assert np.isclose(tw.qx, 62.31, atol=1e-4, rtol=0)
            assert np.isclose(tw.qy, 60.32, atol=1e-4, rtol=0)
            assert np.isclose(tw.dqx, 5, atol=0.1, rtol=0)
            assert np.isclose(tw.dqy, 7, atol=0.1, rtol=0)
        elif line_name == 'lhcb2':
            assert np.isclose(tw.qx, 62.315, atol=1e-4, rtol=0)
            assert np.isclose(tw.qy, 60.325, atol=1e-4, rtol=0)
            assert np.isclose(tw.dqx, 6, atol=0.1, rtol=0)
            assert np.isclose(tw.dqy, 8, atol=0.1, rtol=0)
        else:
            raise ValueError(f'Unknown line name {line_name}')

        assert np.isclose(tw.qs, 0.00212, atol=1e-4, rtol=0) # Checks that RF is well set

        assert np.isclose(tw.c_minus, 0, atol=1e-4, rtol=0)
        assert np.allclose(tw.zeta, 0, rtol=0, atol=1e-4) # Check RF phase

        # Check separations
        assert np.isclose(tw['x', 'ip1'], 0, rtol=0, atol=5e-8) # sigma is 4e-6
        assert np.isclose(tw['y', 'ip1'], 0, rtol=0, atol=5e-8) # sigma is 4e-6
        assert np.isclose(tw['x', 'ip5'], 0, rtol=0, atol=5e-8) # sigma is 4e-6
        assert np.isclose(tw['y', 'ip5'], 0, rtol=0, atol=5e-8) # sigma is 4e-6

        assert np.isclose(tw['x', 'ip2'],
                -0.138e-3 * {'lhcb1': 1, 'lhcb2': 1}[line_name], # set separation
                rtol=0, atol=4e-6)
        assert np.isclose(tw['y', 'ip2'], 0, rtol=0, atol=5e-8)

        assert np.isclose(tw['x', 'ip8'], 0, rtol=0, atol=5e-8)
        assert np.isclose(tw['y', 'ip8'],
                -0.043e-3 * {'lhcb1': 1, 'lhcb2': -1}[line_name], # set separation
                rtol=0, atol=5e-8)

        # Check crossing angles
        assert np.isclose(tw['px', 'ip1'],
                250e-6 * {'lhcb1': 1, 'lhcb2': -1}[line_name], rtol=0, atol=0.5e-6)
        assert np.isclose(tw['py', 'ip1'], 0, rtol=0, atol=0.5e-6)
        assert np.isclose(tw['px', 'ip5'], 0, rtol=0, atol=0.5e-6)
        assert np.isclose(tw['py', 'ip5'], 250e-6, rtol=0, atol=0.5e-6)

        assert np.isclose(tw['px', 'ip2'], 0, rtol=0, atol=0.5e-6)
        assert np.isclose(tw['py', 'ip2'], -100e-6 , rtol=0, atol=0.5e-6) # accounts for spectrometer

        assert np.isclose(tw['px', 'ip8'],
                -115e-6* {'lhcb1': 1, 'lhcb2': -1}[line_name], rtol=0, atol=0.5e-6) # accounts for spectrometer
        assert np.isclose(tw['py', 'ip8'], 2e-6, rtol=0, atol=0.5e-6) # small effect from spectrometer (titled)

        assert np.isclose(tw['betx', 'ip1'], 15e-2, rtol=2e-2, atol=0) # beta beating coming from on_disp
        assert np.isclose(tw['bety', 'ip1'], 15e-2, rtol=3e-2, atol=0)
        assert np.isclose(tw['betx', 'ip5'], 15e-2, rtol=2e-2, atol=0)
        assert np.isclose(tw['bety', 'ip5'], 15e-2, rtol=2e-2, atol=0)

        assert np.isclose(tw['betx', 'ip2'], 10., rtol=4e-2, atol=0)
        assert np.isclose(tw['bety', 'ip2'], 10., rtol=3e-2, atol=0)

        assert np.isclose(tw['betx', 'ip8'], 1.5, rtol=3e-2, atol=0)
        assert np.isclose(tw['bety', 'ip8'], 1.5, rtol=2e-2, atol=0)

        # Check crab cavities
        z_crab_test = 1e-2
        phi_crab_1 = ((
            collider[line_name].twiss(method='4d', zeta0=z_crab_test)['x', 'ip1']
        - collider[line_name].twiss(method='4d', zeta0=-z_crab_test)['x', 'ip1'])
        / 2 / z_crab_test)

        phi_crab_5 = ((
            collider[line_name].twiss(method='4d', zeta0=z_crab_test)['y', 'ip5']
        - collider[line_name].twiss(method='4d', zeta0=-z_crab_test)['y', 'ip5'])
        / 2 / z_crab_test)

        assert np.isclose(phi_crab_1, -190e-6 * {'lhcb1': 1, 'lhcb2': -1}[line_name],
                        rtol=1e-2, atol=0)
        assert np.isclose(phi_crab_5, -170e-6, rtol=1e-2, atol=0)

        # Check one octupole strength
        if line_name == 'lhcb1':
            assert np.isclose(collider['lhcb1']['mo.33l4.b1'].knl[3], -2.2169*200/235,
                          rtol=1e-3, atol=0)
        elif line_name == 'lhcb2':
            assert np.isclose(collider['lhcb2']['mo.33r4.b2'].knl[3], -2.2169,
                            rtol=1e-3, atol=0)



def build_sequence(mad, mylhcbeam, **kwargs):

    # Select beam
    mad.input(f'mylhcbeam = {mylhcbeam}')

    mad.input(

    f'''
    ! Get the toolkit
    call,file=
        "acc-models-lhc/toolkit/macro.madx";
    '''
    '''
    ! Build sequence
    option, -echo,-warn,-info;
    if (mylhcbeam==4){
        call,file="acc-models-lhc/../runIII/lhcb4.seq";
    } else {
        call,file="acc-models-lhc/../runIII/lhc.seq";
    };
    option, -echo, warn,-info;
    '''
    f'''
    !Install HL-LHC
    call, file=
        "acc-models-lhc/hllhc_sequence.madx";
    '''
    '''
    ! Slice nominal sequence
    exec, myslice;
    ''')

    xmlhc.install_errors_placeholders_hllhc(mad)

    mad.input(
    '''
    !Cycling w.r.t. to IP3 (mandatory to find closed orbit in collision in the presence of errors)
    if (mylhcbeam<3){
        seqedit, sequence=lhcb1; flatten; cycle, start=IP3; flatten; endedit;
    };
    seqedit, sequence=lhcb2; flatten; cycle, start=IP3; flatten; endedit;

    ! Install crab cavities (they are off)
    call, file='acc-models-lhc/toolkit/enable_crabcavities.madx';
    on_crab1 = 0;
    on_crab5 = 0;

    ! Set twiss formats for MAD-X parts (macro from opt. toolkit)
    exec, twiss_opt;

    '''
    )


def apply_optics(mad, optics_file):
    mad.call(optics_file)
    # A knob redefinition
    mad.input('on_alice := on_alice_normalized * 7000./nrj;')
    mad.input('on_lhcb := on_lhcb_normalized * 7000./nrj;')

orbit_correction_config = {}
orbit_correction_config['lhcb1'] = {
    'IR1 left': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.r8.b1',
        end='e.ds.l1.b1',
        vary=(
            'corr_co_acbh14.l1b1',
            'corr_co_acbh12.l1b1',
            'corr_co_acbv15.l1b1',
            'corr_co_acbv13.l1b1',
            ),
        targets=('e.ds.l1.b1',),
    ),
    'IR1 right': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.r1.b1',
        end='s.ds.l2.b1',
        vary=(
            'corr_co_acbh13.r1b1',
            'corr_co_acbh15.r1b1',
            'corr_co_acbv12.r1b1',
            'corr_co_acbv14.r1b1',
            ),
        targets=('s.ds.l2.b1',),
    ),
    'IR5 left': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.r4.b1',
        end='e.ds.l5.b1',
        vary=(
            'corr_co_acbh14.l5b1',
            'corr_co_acbh12.l5b1',
            'corr_co_acbv15.l5b1',
            'corr_co_acbv13.l5b1',
            ),
        targets=('e.ds.l5.b1',),
    ),
    'IR5 right': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.r5.b1',
        end='s.ds.l6.b1',
        vary=(
            'corr_co_acbh13.r5b1',
            'corr_co_acbh15.r5b1',
            'corr_co_acbv12.r5b1',
            'corr_co_acbv14.r5b1',
            ),
        targets=('s.ds.l6.b1',),
    ),
    'IP1': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.l1.b1',
        end='s.ds.r1.b1',
        vary=(
            'corr_co_acbch6.l1b1',
            'corr_co_acbcv5.l1b1',
            'corr_co_acbch5.r1b1',
            'corr_co_acbcv6.r1b1',
            'corr_co_acbyhs4.l1b1',
            'corr_co_acbyhs4.r1b1',
            'corr_co_acbyvs4.l1b1',
            'corr_co_acbyvs4.r1b1',
        ),
        targets=('ip1', 's.ds.r1.b1'),
    ),
    'IP2': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.l2.b1',
        end='s.ds.r2.b1',
        vary=(
            'corr_co_acbyhs5.l2b1',
            'corr_co_acbchs5.r2b1',
            'corr_co_acbyvs5.l2b1',
            'corr_co_acbcvs5.r2b1',
            'corr_co_acbyhs4.l2b1',
            'corr_co_acbyhs4.r2b1',
            'corr_co_acbyvs4.l2b1',
            'corr_co_acbyvs4.r2b1',
        ),
        targets=('ip2', 's.ds.r2.b1'),
    ),
    'IP5': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.l5.b1',
        end='s.ds.r5.b1',
        vary=(
            'corr_co_acbch6.l5b1',
            'corr_co_acbcv5.l5b1',
            'corr_co_acbch5.r5b1',
            'corr_co_acbcv6.r5b1',
            'corr_co_acbyhs4.l5b1',
            'corr_co_acbyhs4.r5b1',
            'corr_co_acbyvs4.l5b1',
            'corr_co_acbyvs4.r5b1',
        ),
        targets=('ip5', 's.ds.r5.b1'),
    ),
    'IP8': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.l8.b1',
        end='s.ds.r8.b1',
        vary=(
            'corr_co_acbch5.l8b1',
            'corr_co_acbyhs4.l8b1',
            'corr_co_acbyhs4.r8b1',
            'corr_co_acbyhs5.r8b1',
            'corr_co_acbcvs5.l8b1',
            'corr_co_acbyvs4.l8b1',
            'corr_co_acbyvs4.r8b1',
            'corr_co_acbyvs5.r8b1',
        ),
        targets=('ip8', 's.ds.r8.b1'),
    ),
}

orbit_correction_config['lhcb2'] = {
    'IR1 left': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.l1.b2',
        end='e.ds.r8.b2',
        vary=(
            'corr_co_acbh13.l1b2',
            'corr_co_acbh15.l1b2',
            'corr_co_acbv12.l1b2',
            'corr_co_acbv14.l1b2',
            ),
        targets=('e.ds.r8.b2',),
    ),
    'IR1 right': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.l2.b2',
        end='s.ds.r1.b2',
        vary=(
            'corr_co_acbh12.r1b2',
            'corr_co_acbh14.r1b2',
            'corr_co_acbv13.r1b2',
            'corr_co_acbv15.r1b2',
            ),
        targets=('s.ds.r1.b2',),
    ),
    'IR5 left': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='e.ds.l5.b2',
        end='e.ds.r4.b2',
        vary=(
            'corr_co_acbh13.l5b2',
            'corr_co_acbh15.l5b2',
            'corr_co_acbv12.l5b2',
            'corr_co_acbv14.l5b2',
            ),
        targets=('e.ds.r4.b2',),
    ),
    'IR5 right': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.l6.b2',
        end='s.ds.r5.b2',
        vary=(
            'corr_co_acbh12.r5b2',
            'corr_co_acbh14.r5b2',
            'corr_co_acbv13.r5b2',
            'corr_co_acbv15.r5b2',
            ),
        targets=('s.ds.r5.b2',),
    ),
    'IP1': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.r1.b2',
        end='e.ds.l1.b2',
        vary=(
            'corr_co_acbch6.r1b2',
            'corr_co_acbcv5.r1b2',
            'corr_co_acbch5.l1b2',
            'corr_co_acbcv6.l1b2',
            'corr_co_acbyhs4.l1b2',
            'corr_co_acbyhs4.r1b2',
            'corr_co_acbyvs4.l1b2',
            'corr_co_acbyvs4.r1b2',
        ),
        targets=('ip1', 'e.ds.l1.b2',),
    ),
    'IP2': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.r2.b2',
        end='e.ds.l2.b2',
        vary=(
            'corr_co_acbyhs5.l2b2',
            'corr_co_acbchs5.r2b2',
            'corr_co_acbyvs5.l2b2',
            'corr_co_acbcvs5.r2b2',
            'corr_co_acbyhs4.l2b2',
            'corr_co_acbyhs4.r2b2',
            'corr_co_acbyvs4.l2b2',
            'corr_co_acbyvs4.r2b2',
        ),
        targets=('ip2', 'e.ds.l2.b2'),
    ),
    'IP5': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.r5.b2',
        end='e.ds.l5.b2',
        vary=(
            'corr_co_acbch6.r5b2',
            'corr_co_acbcv5.r5b2',
            'corr_co_acbch5.l5b2',
            'corr_co_acbcv6.l5b2',
            'corr_co_acbyhs4.l5b2',
            'corr_co_acbyhs4.r5b2',
            'corr_co_acbyvs4.l5b2',
            'corr_co_acbyvs4.r5b2',
        ),
        targets=('ip5', 'e.ds.l5.b2',),
    ),
    'IP8': dict(
        ref_with_knobs={'on_corr_co': 0, 'on_disp': 0},
        start='s.ds.r8.b2',
        end='e.ds.l8.b2',
        vary=(
            'corr_co_acbchs5.l8b2',
            'corr_co_acbyhs5.r8b2',
            'corr_co_acbcvs5.l8b2',
            'corr_co_acbyvs5.r8b2',
            'corr_co_acbyhs4.l8b2',
            'corr_co_acbyhs4.r8b2',
            'corr_co_acbyvs4.l8b2',
            'corr_co_acbyvs4.r8b2',
        ),
        targets=('ip8', 'e.ds.l8.b2',),
    ),
}

knob_settings_yaml_str = """
  # Orbit knobs
  on_x1: 250            # [urad]
  on_sep1: 0            # [mm]
  on_x2: -170           # [urad]
  on_sep2: 0.138        # [mm]
  on_x5: 250            # [urad]
  on_sep5: 0            # [mm]
  on_x8: -250           # [urad]
  on_sep8: -0.043       # [mm]
  on_a1: 0              # [urad]
  on_o1: 0              # [mm]
  on_a2: 0              # [urad]
  on_o2: 0              # [mm]
  on_a5: 0              # [urad]
  on_o5: 0              # [mm]
  on_a8: 0              # [urad]
  on_o8: 0              # [mm]
  on_disp: 1            # Value to choose could be optics-dependent

  # Crab cavities
  on_crab1: -190        # [urad]
  on_crab5: -170        # [urad]

  # Magnets of the experiments
  on_alice_normalized: 1
  on_lhcb_normalized: 1
  on_sol_atlas: 0
  on_sol_cms: 0
  on_sol_alice: 0

  # RF voltage and phases
  vrf400:       16.0            # [MV]
  lagrf400.b1:   0.5            # [rad]
  lagrf400.b2:   0.             # [rad]

  # Octupoles
  i_oct_b1:     -200            # [A]
  i_oct_b2:     -235            # [A]
"""

knob_names_yaml_str = """
lhcb1:
    q_knob_1: kqtf.b1
    q_knob_2: kqtd.b1
    dq_knob_1: ksf.b1
    dq_knob_2: ksd.b1
    c_minus_knob_1: c_minus_re_b1
    c_minus_knob_2: c_minus_im_b1

lhcb2:
    q_knob_1: kqtf.b2
    q_knob_2: kqtd.b2
    dq_knob_1: ksf.b2
    dq_knob_2: ksd.b2
    c_minus_knob_1: c_minus_re_b2
    c_minus_knob_2: c_minus_im_b2
"""
tune_chroma_yaml_str = """
qx:
  lhcb1: 62.31
  lhcb2: 62.315
qy:
  lhcb1: 60.32
  lhcb2: 60.325
dqx:
  lhcb1: 5
  lhcb2: 6
dqy:
  lhcb1: 7
  lhcb2: 8

"""

def _get_z_centroids(ho_slices, sigmaz):
    from scipy.stats import norm
    z_cuts = norm.ppf(
        np.linspace(0, 1, ho_slices + 1)[1:int((ho_slices + 1) / 2)]) * sigmaz
    z_centroids = []
    z_centroids.append(-sigmaz / np.sqrt(2*np.pi)
        * np.exp(-z_cuts[0]**2 / (2 * sigmaz * sigmaz)) * float(ho_slices))
    for ii,jj in zip(z_cuts[0:-1],z_cuts[1:]):
        z_centroids.append(-sigmaz / np.sqrt(2*np.pi)
            * (np.exp(-jj**2 / (2 * sigmaz * sigmaz))
               - np.exp(-ii**2 / (2 * sigmaz * sigmaz))
            ) * ho_slices)
    return np.array(z_centroids + [0] + [-ii for ii in z_centroids[-1::-1]])