import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c as clight
from scipy.constants import e as qe
from scipy.constants import epsilon_0

import xtrack as xt

num_particles = 1e11
nemitt_x = 2.5e-6
nemitt_y = 2.5e-6

collider = xt.Multiline.from_json('collider_00_from_mad.json')
collider.build_trackers()

# Switch off the crab cavities
collider.vars['on_crab1'] = 0
collider.vars['on_crab5'] = 0


# Check that orbit is flat
tw = collider.twiss(method='4d')
assert np.max(np.abs(tw.lhcb1.x))< 1e-7

# Install head-on only
collider.discard_trackers()

collider.install_beambeam_interactions(
    clockwise_line='lhcb1',
    anticlockwise_line='lhcb2',
    ip_names=['ip1'],
    num_long_range_encounters_per_side=[0],
    num_slices_head_on=1,
    harmonic_number=35640,
    bunch_spacing_buckets=10,
    sigmaz=1e-6
)

collider.build_trackers()
# Switch on RF (assumes 6d)
collider.vars['vrf400'] = 16

collider.configure_beambeam_interactions(
    num_particles=num_particles,
    nemitt_x=nemitt_x, nemitt_y=nemitt_y,
    crab_strong_beam=False
)

collider.lhcb1.matrix_stability_tol = 1e-2
collider.lhcb2.matrix_stability_tol = 1e-2

tw_bb_on = collider.lhcb1.twiss()
collider.vars['beambeam_scale'] = 0
tw_bb_off = collider.lhcb1.twiss()

tune_shift_x = tw_bb_on.qx - tw_bb_off.qx
tune_shift_y = tw_bb_on.qy - tw_bb_off.qy

# Analytical tune shift

q0 = collider.lhcb1.particle_ref.q0
mass0 = collider.lhcb1.particle_ref.mass0 # eV

# classical particle radius
r0 = 1 / (4 * np.pi * epsilon_0) * q0**2 * qe / mass0

plt.close('all')
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
collider.vars['beambeam_scale'] = 1
fp_b1_bb_on = collider.lhcb1.get_footprint(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
fp_b1_bb_on.plot(ax=ax1, color='k')
plt.plot(np.mod(tw_bb_on.qx, 1), np.mod(tw_bb_on.qy, 1), 'ko', markersize=10)

collider.vars['beambeam_scale'] = 0
fp_b2_bb_off = collider.lhcb1.get_footprint(nemitt_x=nemitt_x, nemitt_y=nemitt_y)
fp_b2_bb_off.plot(ax=ax1, color='g')
plt.plot(np.mod(tw_bb_off.qx, 1), np.mod(tw_bb_off.qy, 1), 'go', markersize=10)


plt.show()


