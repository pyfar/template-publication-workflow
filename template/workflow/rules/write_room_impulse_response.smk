

rule write_room_impulse_response:
    output:
        "results/example.rir.sofa"
    log:
        "logs/write_room_impulse_response.log"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/write_room_impulse_response.py"
