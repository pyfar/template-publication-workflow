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
    input_file = os.path.join('results', 'example.rir.sofa')
    output_file = os.path.join('results', 'example.rt.sofa')

else:
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

# %%
# Create output folder
input_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', input_file))
output_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', output_file))
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# %%
# Load example RIR
sofa = sf.read_sofa(input_path)

rir, source_positions, receiver_positions = pf.io.convert_sofa(sofa)

# %%
# Third-octave band filtering

rir_third_octave = pf.dsp.filter.fractional_octave_bands(
    rir,
    num_fractions=3
)

# %%
# Calculate energy decay curves

edc = pr.schroeder_integration(rir_third_octave)

# %%
# A first impression of the data can be obtained with

sofa.inspect()

# %%
# Fill SOFA Data
# check out the documentation of the sofa convention to first fill the data
# https://sofar.readthedocs.io/en/stable/resources/conventions.html#generalfir-1-0

sofa.Data_IR = rir.time[np.newaxis, :, :]

# %%
# fill metadata
# next check out the required and optional metadata for the GeneralFIR convention
# and fill them accordingly


# %%
# check if sofa data are consistence and meet the standard
sofa.verify()


# %%
# write sofa file
sf.write_sofa(output_file, sofa)

# %%
