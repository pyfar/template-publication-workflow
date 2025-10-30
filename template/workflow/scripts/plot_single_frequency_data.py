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
    input_file = os.path.join('..', '..', 'results', 'example.rt30.sofa')
    output_file = os.path.join('..', '..', 'results', 'example.rt30.png')
else:
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]

# %%
# Create output folder
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# %%
# Load example RIR
sofa = sf.read_sofa(input_file)

rts, source_positions, receiver_positions = pf.io.convert_sofa(sofa)

# %%
fig = plt.figure()
ax = pf.plot.freq(rts, dB=False)
ax.set_ylabel('Reverberation Time (s)')
ax.set_ylim((0.2, 0.4))

fig.savefig(output_file)

# %%