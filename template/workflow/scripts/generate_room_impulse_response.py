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
    output_file = os.path.join('..', '..', 'resources', 'ShoeboxRoom_receiver_3_2_1.rir.sofa')
    r_x, r_y, r_z = 3, 2, 1  # in meters
else:
    output_file = snakemake.output[0]
    r_x = float(snakemake.params.receiver_position_x)
    r_y = float(snakemake.params.receiver_position_y)
    r_z = float(snakemake.params.receiver_position_z)


# %%
# Create output folder
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# %%
# Define room parameters and generate room impulse response

room_dimensions = [7, 4, 3]  # in meters
source_position = [6, 1.5, 1.8]  # in meters
receiver_position = [r_x, r_y, r_z]  # in meters

assert all(0 < v1 < v2 for v1, v2 in zip(receiver_position, room_dimensions)), \
    "Receiver position must be inside the room."

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
    n_samples=2**14,
)[0]

noise_power = np.max(np.abs(room_impulse_response.time), axis=-1) * 10**(-60/20)
awgn = pf.signals.noise(
    room_impulse_response.n_samples, spectrum='white', rms=noise_power,
    sampling_rate=room_impulse_response.sampling_rate)

noisy_room_impulse_response = room_impulse_response + awgn

# %%
sofa = sf.Sofa('GeneralFIR')

# %%
# A first impression of the data can be obtained with

sofa.inspect()

# %%
# Fill SOFA Data
# check out the documentation of the sofa convention to first fill the data
# https://sofar.readthedocs.io/en/stable/resources/conventions.html#generalfir-1-0

sofa.Data_IR = noisy_room_impulse_response.time[np.newaxis, :, :]
sofa.Data_SamplingRate = noisy_room_impulse_response.sampling_rate

# %%
# Fill metadata
# next check out the required and optional metadata for the GeneralFIR convention
# and fill them accordingly
# You can also add custom metadata if needed. 
# For reference, have a look at the SOFA conventions SingleRoomSRIR which is
# intended for storing room impulse responses with a large number of receivers.


# %%
# check if sofa data are consistent and meet the standard
sofa.verify()


# %%
# write sofa file
sf.write_sofa(output_file, sofa)

# %%
