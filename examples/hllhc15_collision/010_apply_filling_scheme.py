import numpy as np

import xtrack as xt

collider = xt.Multiline.from_json('./collider_04_tuned_and_leveled_bb_on.json')
collider.build_trackers()

filling_pattern_b1 = np.zeros(3564, dtype=int)
filling_pattern_b2 = np.zeros(3564, dtype=int)

# Fill 50 bunches around bunch 500
filling_pattern_b1[500-25:500+25] = 1
filling_pattern_b2[500-25:500+25] = 1


harmonic_number = 35640
bunch_spacing_buckets = 10

ip_names = collider._bb_config['ip_names'] # is ['ip1', 'ip2', 'ip5', 'ip8']
delay_at_ips_slots = [0, 891, 0, 2670] # Defined as anticlockwise bunch id that
                                       # meets bunch 0 of the clockwise beam

ring_length_in_slots = harmonic_number / bunch_spacing_buckets

for orientation in ['clockwise', 'anticlockwise']:

    if orientation == 'clockwise':
        delay_at_ips_dict = {iipp: dd
                             for iipp, dd in zip(ip_names, delay_at_ips_slots)}
    elif orientation == 'anticlockwise':
        delay_at_ips_dict = {iipp: np.mod(ring_length_in_slots - dd, ring_length_in_slots)
                             for iipp, dd in zip(ip_names, delay_at_ips_slots)}
    else:
        raise ValueError('?!')

    bbdf = collider._bb_config['dataframes'][orientation]

    delay_in_slots = []

    for nn in bbdf.index.values:
        ip_name = bbdf.loc[nn, 'ip_name']
        this_delay = delay_at_ips_dict[ip_name]

        if nn.startswith('bb_lr.'):
            if orientation == 'clockwise':
                this_delay += bbdf.loc[nn, 'identifier']
            elif orientation == 'anticlockwise':
                this_delay -= bbdf.loc[nn, 'identifier']
            else:
                raise ValueError('?!')

        delay_in_slots.append(this_delay)

    bbdf['delay_in_slots'] = delay_in_slots


# Some checks
dframes = collider._bb_config['dataframes']
assert (dframes['clockwise'].loc['bb_ho.c1b1_00', 'delay_in_slots'] == 0)
assert (dframes['clockwise'].loc['bb_ho.c5b1_00', 'delay_in_slots'] == 0)
assert (dframes['clockwise'].loc['bb_ho.c2b1_00', 'delay_in_slots'] == 891)
assert (dframes['clockwise'].loc['bb_ho.c8b1_00', 'delay_in_slots'] == 2670)

assert (dframes['anticlockwise'].loc['bb_ho.c1b2_00', 'delay_in_slots'] == 0)
assert (dframes['anticlockwise'].loc['bb_ho.c5b2_00', 'delay_in_slots'] == 0)
assert (dframes['anticlockwise'].loc['bb_ho.c2b2_00', 'delay_in_slots'] == 3564 - 891)
assert (dframes['anticlockwise'].loc['bb_ho.c8b2_00', 'delay_in_slots'] == 3564 - 2670)
