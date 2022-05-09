#!/usr/bin/env snakemake

import os
import pandas as pd
from pprint import pprint
from glob import glob


GENOMES = pd.read_csv(config["input"]["genomes"], sep="\t").set_index("genome_dir")
#GENOMES = pd.read_csv(config["input"]["genomes"], sep="\t").set_index("genome_id")

ENZYMES = ["Full", "2bRAD",
           "AlfI", "AloI", "BaeI", "BcgI",
           "BplI", "BsaXI", "BslFI", "Bsp24I",
           "CjeI", "CjePI", "CspCI", "FalI",
           "HaeIV", "Hin4I", "PpiI", "PsrI"]

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
	ln -s {input} {output}/{params.out_prefix}.Full.fa.gz
        '''


checkpoint mash_sketch_prepare:
    input:
        genomes = expand(os.path.join(config["output"]["2bRAD_M"],
                                      "database/{genome}"),
                         genome=GENOMES.index.unique())
    output:
        list_dir = directory(os.path.join(config["output"]["2bRAD_M"],
                                          "genome_list/{enzyme}"))
    params:
        batch_num = config["params"]["mash"]["batch_num"],
        enzyme = "{enzyme}"
    run:
        import os
        fa_list = []
        for genome_dir in input.genomes:
            fa_list_ = glob(f'''{genome_dir}/*.{params.enzyme}.fa.gz''')
            for fa in fa_list_:
                if os.path.getsize(fa) > 20:
                    fa_list.append(fa)

        os.makedirs(output.list_dir, exist_ok=True)

        for batch in range(0, len(fa_list), params.batch_num):
            batch_file = os.path.join(output.list_dir, f"batch_{batch}.tsv")
            with open(batch_file, "w") as oh:
                oh.write("\n".join(fa_list[batch:batch+params.batch_num]))


rule mash_sketch:
    input:
        os.path.join(config["output"]["2bRAD_M"], "genome_list/{enzyme}/batch_{batch}.tsv")
    output:
        os.path.join(config["output"]["2bRAD_M"], "mash_sketch/{enzyme}/batch_{batch}.msh")
    wildcard_constraints:
        batch = "\d+"
    benchmark:
        os.path.join(config["output"]["2bRAD_M"],
                     "benchmark/mash_sketch/{enzyme}/batch_{batch}.mash_sketch.benchmark.txt")
    log:
        os.path.join(config["output"]["2bRAD_M"],
                     "logs/mash_sketch/{enzyme}/batch_{batch}.mash_sketch.log")
    threads:
        config["params"]["mash"]["threads"]
    params:
        out_prefix = os.path.join(config["output"]["2bRAD_M"], "mash_sketch/{enzyme}/batch_{batch}"),
        kmer_size = config["params"]["mash"]["kmer_size"],
        sketch_size = config["params"]["mash"]["sketch_size"],
        seed = config["params"]["mash"]["seed"],
        probability_threshold = config["params"]["mash"]["probability_threshold"]
    shell:
        '''
        mash sketch \
        -p {threads} \
        -l {input} \
        -o {params.out_prefix} \
        -k {params.kmer_size} \
        -s {params.sketch_size} \
        -S {params.seed} \
        -w {params.probability_threshold} \
        >{log} 2>&1
        '''
   
rule mash_dist:
    input: 
        reference = os.path.join(config["output"]["2bRAD_M"],
                                 "mash_sketch/{enzyme}/batch_{reference}.msh"),
        query = os.path.join(config["output"]["2bRAD_M"],
                             "mash_sketch/{enzyme}/batch_{query}.msh")
    output:
        dist = temp(os.path.join(config["output"]["2bRAD_M"],
                            "mash_dist/{enzyme}/{reference}_vs_{query}.mash_dist.txt.gz"))
    benchmark:
        os.path.join(config["output"]["2bRAD_M"],
                     "benchmark/mash_dist/{enzyme}/{reference}_vs_{query}.mash_dist.benchmark.txt")
    log:
        os.path.join(config["output"]["2bRAD_M"],
                     "logs/mash_dist/{enzyme}/{reference}_vs_{query}.mash_dist.log")
    threads:
        config["params"]["mash"]["threads"]
    params:
        max_pvalue = config["params"]["mash"]["max_pvalue"],
        max_distance = config["params"]["mash"]["max_distance"]
    shell:
        '''
        mash dist \
        -p {threads} \
        -v {params.max_pvalue} \
        -d {params.max_distance} \
        {input.reference} {input.query} \
        > {output.dist}.tmp 2> {log}
        
        cat {output.dist}.tmp | \
        sed 's#.fa.gz##g' | \
        awk -F'\t' '{{split($1,a,"/"); split($2,b,"/"); print a[length(a)] "\t" b[length(b)] "\t" $3 "\t" $4 "\t" $5}}' | \
        gzip -c > {output.dist}

        rm -rf {output.dist}.tmp 
        '''


def aggregate_mash_dist_output(wildcards):
    checkpoint_output = checkpoints.mash_sketch_prepare.get(**wildcards).output.list_dir
    batch=sorted(list(set([i for i in glob_wildcards(os.path.join(checkpoint_output, "batch_{batch}.tsv")).batch])))
    #pprint(batch)
    #return expand(os.path.join(config["output"]["2bRAD_M"],
    #                           "mash_dist/{enzyme}/{reference}_vs_{query}.mash_dist.txt"),
    #              enzyme=wildcards.enzyme,
    #              reference=sorted(list(set([i for i in glob_wildcards(os.path.join(checkpoint_output, "batch_{batch}.tsv")).batch]))),
    #              query=sorted(list(set([i for i in glob_wildcards(os.path.join(checkpoint_output, "batch_{batch}.tsv")).batch]))))
   
    # reduce the number of pair-wise comparison 
    reference_list = []
    query_list = []
    enzyme_list = []
    for i in range(0, len(batch)):
        for j in range(i, len(batch)):
            reference_list.append(batch[i])
            query_list.append(batch[j])
	    enzyme_list.append(wildcards.enzyme)

    return expand(os.path.join(config["output"]["2bRAD_M"],
                               "mash_dist/{enzyme}/{reference}_vs_{query}.mash_dist.txt.gz"),
                  zip,
                  enzyme=enzyme_list,
                  reference=reference_list,
                  query=query_list)


rule mash_dist_merge:
    input:
        aggregate_mash_dist_output
    output:
        os.path.join(config["output"]["2bRAD_M"], "mash_dist_merge/{enzyme}.mash_dist.txt.gz")
    shell:
        '''
        cat {input} > {output}
        '''


rule check_finished:
    input:
        expand(os.path.join(config["output"]["2bRAD_M"], "mash_dist_merge/{enzyme}.mash_dist.txt.gz"),
               enzyme=ENZYMES)
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
    mash_sketch_prepare,
    mash_dist_merge,
    check_finished,
    all
