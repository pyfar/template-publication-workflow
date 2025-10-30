

rule plot_single_frequency_data:
    input:
        "results/{filename}.rt{rt}.sofa"
    output:
        "figures/{filename}.rt{rt}.{postfix}"
    log:
        "logs/plot_single_frequency_data.{filename}.rt{rt}.{postfix}.log"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_single_frequency_data.py"
