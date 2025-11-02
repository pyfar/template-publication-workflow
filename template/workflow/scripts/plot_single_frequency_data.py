# %%

import numpy as np
import pyfar as pf
import os
import matplotlib.pyplot as plt

input_file = os.path.join('..', '..', 'results', 'ShoeboxRoom_receiver_3_2_1.rt30.csv')
output_file = os.path.join('..', '..', 'results', 'ShoeboxRoom_receiver_3_2_1.rt30.png')

# %%
# Set the plot style for a manuscript (sets font sizes, latex rendering, etc..)
plt.style.use('../../config/manuscript.mplstyle')
# %%
# Create output folder
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# %%
# Load example RIR
data_raw = np.loadtxt(input_file)
center_frequencies = data_raw[0, :]
reverberation_times = data_raw[1, :]
reverberation_times = pf.FrequencyData(reverberation_times, center_frequencies) 
# %%
fig = plt.figure()
ax = pf.plot.freq(reverberation_times, dB=False)
ax.set_ylabel('Reverberation Time (s)')
ax.set_ylim((0.2, 0.4))

fig.savefig(output_file)

# %%