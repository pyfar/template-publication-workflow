# %%

import numpy as np
import sofar as sf
import pyfar as pf
import os
from pyrato.analytic import rectangular_room_rigid_walls


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

room_dimensions = [7, 4, 3]  # in meters

source_position = [6, 1.5, 1.8]  # in meters
receiver_position = [3, 2.5, 1.5]  # in meters

reverberation_time_target = 0.3  # in seconds

temperature = 20  # in degree Celsius
humidity = 50  # in percent

speed_of_sound = 331.3 + 0.606 * temperature + 0.0018 * humidity  # in m/s

sampling_rate = 8e3  # in Hz

# simulate a room impulse response
room_impulse_response = rectangular_room_rigid_walls(
    dimensions=room_dimensions,
    source=source_position,
    receiver=receiver_position,
    reverberation_time=reverberation_time_target,
    max_freq=3e3,
    speed_of_sound=speed_of_sound,
    samplingrate=sampling_rate,
    n_samples=2**12,
)[0]

noise_power = np.max(np.abs(room_impulse_response.time), axis=-1) * 10**(-60/20)
awgn = pf.signals.noise(
    room_impulse_response.n_samples, spectrum='white', rms=noise_power,
    sampling_rate=room_impulse_response.sampling_rate)

noisy_room_impulse_response = room_impulse_response + awgn

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

sofa.Data_IR = room_impulse_response.time[np.newaxis, :, :]

# %%
# Fill metadata
# next check out the required and optional metadata for the GeneralFIR convention
# and fill them accordingly
# You can also add custom metadata if needed. 
# For reference, have a look at the SOFA conventions SingleRoomSRIR which is
# intended for storing room impulse responses with a large number of receivers.


# %%
# check if sofa data are consistence and meet the standard
sofa.verify()


# %%
# write sofa file
sf.write_sofa(output_file, sofa)

# %%
