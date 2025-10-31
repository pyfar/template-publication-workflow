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
    input_file = os.path.join('..', '..', 'resources', 'shoebox_room_receiver_3_2.5_1.5.rir.sofa')
    output_file = os.path.join('..', '..', 'results', 'shoebox_room_receiver_3_2.5_1.5.rt30.sofa')
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
# edc = pf.dsp.normalize(edc, reference_method='max')

# %%
# Calculate the reverberation time

reverberation_time = pr.reverberation_time_linear_regression(edc, T=f'T{rt}')
print(reverberation_time)
reverberation_time = pf.FrequencyData(reverberation_time, center_frequencies)


# %%
# Create SOFA object
sofa_out = sf.Sofa('GeneralTF')

# %%
# A first impression of the data can be obtained with

sofa_out.inspect()

# %%
# Fill SOFA Data
# check out the documentation of the sofa convention to first fill the data
# https://sofar.readthedocs.io/en/stable/resources/conventions.html#generaltf-1-0

sofa_out.Data_Real = np.real(reverberation_time.freq[np.newaxis, :])
sofa_out.Data_Imag = np.imag(reverberation_time.freq[np.newaxis, :])
sofa_out.N = center_frequencies

# %%
# fill metadata
# next check out the required and optional metadata for the GeneralTF convention
# and fill them accordingly


# %%
# check if sofa data are consistent and meet the standard
sofa_out.verify()


# %%
# write sofa file
sf.write_sofa(output_file, sofa_out)

# %%
