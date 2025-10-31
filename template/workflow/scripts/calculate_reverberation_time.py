# %%

import numpy as np
import sofar as sf
import pyfar as pf
import pyrato as pr
import os


try: 
    snakemake_exists = hasattr(snakemake, 'input')
except NameError:

    snakemake_exists = False

if not snakemake_exists:
    input_file = os.path.join('..', '..', 'resources', 'ShoeboxRoom_receiver_3_2_1.rir.sofa')
    output_file = os.path.join('..', '..', 'results', 'ShoeboxRoom_receiver_3_2_1.rt30.csv')
    rt = '30'

else:
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    rt = snakemake.params.rt

# %%
# Create output folder
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# %%
# Load RIR and convert from SOFA format
sofa = sf.read_sofa(input_file)

rir, source_positions, receiver_positions = pf.io.convert_sofa(sofa)

# %%
# Octave band filtering
frequency_range = (100, 1e3)
num_octave_fractions = 1
center_frequencies = pf.dsp.filter.fractional_octave_frequencies(
    num_fractions=num_octave_fractions,
    frequency_range=frequency_range,
    )[0]

rir_octave = pf.dsp.filter.fractional_octave_bands(
    rir,
    num_fractions=num_octave_fractions,
    frequency_range=frequency_range,
)[:, 0, 0]

# %%
# Calculate energy decay curves

edc = pr.energy_decay_curve_chu(rir_octave, channel_independent=True)

# %%
# Calculate the reverberation time

reverberation_time = pr.reverberation_time_linear_regression(edc, T=f'T{rt}')

# %%
# write reverberation time to disk

data_write = np.vstack((center_frequencies, reverberation_time))

np.savetxt(output_file, data_write)

# %%
