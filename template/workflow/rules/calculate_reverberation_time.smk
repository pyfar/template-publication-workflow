

rule calculate_reverberation_time:
    input:
        "resources/{filename}.rir.sofa"
    output:
        "results/{filename}.rt{rt}.sofa"
    params:
        rt="{rt}"  # which decay time to calculate (e.g., 60 for RT60)
    log:
        "logs/calculate_reverberation_time.{filename}.rt{rt}.log"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/calculate_reverberation_time.py"
