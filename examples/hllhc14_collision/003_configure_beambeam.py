import xmask as xm
import xtrack as xt

# Read beam-beam config from config file
with open('config.yaml','r') as fid:
    config = xm.yaml.load(fid)
config_bb = config['config_beambeam']

# Load collider and build trackers
collider = xt.Multiline.from_json('collider_02_tuned_bb_off.json')
collider.build_trackers()

# Configure beam-beam lenses
print('Configuring beam-beam lenses...')
collider.configure_beambeam_interactions(
    num_particles=config_bb['num_particles_per_bunch'],
    nemitt_x=config_bb['nemitt_x'],
    nemitt_y=config_bb['nemitt_y'])

collider.to_json('collider_03_tuned_bb_on.json')
