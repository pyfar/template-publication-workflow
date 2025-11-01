

rule generate_room_impulse_response_shoebox_room:
    output:
        "resources/ShoeboxRoom_receiver_{r_x}_{r_y}_{r_z}.rir.sofa"
    params:
        receiver_position_x = "{r_x}",
        receiver_position_y = "{r_y}",
        receiver_position_z = "{r_z}",
    log:
        "logs/generate_room_impulse_response.{r_x}_{r_y}_{r_z}.log"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/generate_room_impulse_response.py"
