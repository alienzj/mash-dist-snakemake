#!/usr/bin/env snakemake

import os
import pandas as pd
from pprint import pprint


GENOMES = pd.read_csv(config["input"]["genomes"], sep="\t").set_index("genome_dir")
#GENOMES = pd.read_csv(config["input"]["genomes"], sep="\t").set_index("genome_id")


rule extract_2btag:
    input:
    	lambda wildcards: str(GENOMES.loc[wildcards.genome]["genome_path"])
    output:
        directory(os.path.join(config["output"]["2bRAD_M"], "database/{genome}"))
    wildcard_constraints:
        genome = "GC[A|F]/\d{3}/\d{3}/\d{3}"
    log:
        os.path.join(config["output"]["2bRAD_M"], "logs/{genome}", "extract_2btag.log")
    benchmark:
        os.path.join(config["output"]["2bRAD_M"], "benchmark/{genome}", "extract_2btag.txt")
    params:
        IIBRADM = config["params"]["2bRAD_M"]["path"],
	log_dir = os.path.join(config["output"]["2bRAD_M"], "logs/{genome}"),
	bench_dir = os.path.join(config["output"]["2bRAD_M"], "benchmark/{genome}"),
        out_prefix = lambda wildcards: GENOMES.loc[wildcards.genome, "genome_id"]
    threads: 1
    shell:
        '''
	mkdir -p {params.log_dir}
	mkdir -p {params.bench_dir}
        mkdir -p {output}

        for s in {{1..16}}
        do
            perl {params.IIBRADM}/scripts/2bRADExtraction.pl \
            -i {input} \
            -t 1 \
            -s $s \
            -od {output} \
            -op {params.out_prefix} \
	    >> {log} 2>&1
        done

        cat {output}/{params.out_prefix}.*.fa.gz > {output}/{params.out_prefix}.2bRAD.fa.gz
        '''


rule check_finished:
    input:
        expand(os.path.join(config["output"]["2bRAD_M"],
                            "database/{genome}"),
               genome=GENOMES.index.unique())
    output:
        os.path.join(config["output"]["2bRAD_M"], "finished")
    shell:
        '''
	touch {output}
	'''

rule all:
    input:
        os.path.join(config["output"]["2bRAD_M"], "finished")


localrules:
    check_finished,
    all
