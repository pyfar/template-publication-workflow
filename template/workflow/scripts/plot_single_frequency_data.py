# %%

import sofar as sf
import pyfar as pf
import os
import matplotlib.pyplot as plt


try: 
    snakemake_exists = hasattr(snakemake, 'input')
except NameError:

    snakemake_exists = False

if not snakemake_exists:
    input_file = os.path.join('results', 'example.rt.sofa')
    output_file = os.path.join('results', 'example.rt.png')

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

rts, source_positions, receiver_positions = pf.io.convert_sofa(sofa)

# %%
fig = plt.figure()
pf.plot.freq(rts, dB=False)

fig.savefig(output_path)

# %%