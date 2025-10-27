# %%

import numpy as np
import sofar as sf
import pyfar as pf
import os


try: 
    snakemake_exists = hasattr(snakemake, 'input')
except NameError:

    snakemake_exists = False

if not snakemake_exists:
    output_file = os.path.join('results', 'example.rir.sofa')

else:
    output_file = snakemake.output[0]

# %%
# Create output folder
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# %%
# Load example RIR
rir = pf.signals.files.room_impulse_response()

# %%
# SOFA is a standard format to store spatial acoustic data, such as
# room impulse responses (RIRs). Here you can find an overview of
# the SOFA standards https://www.sofaconventions.org/

# we use sofar (https://sofar.readthedocs.io/) to create and write SOFA files
# and create a GeneralFIR SOFA object to store the RIR and required metadata.
sofa = sf.Sofa('GeneralFIR')

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
