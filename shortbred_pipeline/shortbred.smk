
from os.path import join
from glob import glob

# Load configurations
configfile: "config.yaml"

# Auto-generate samples list
SAMPLES, = glob_wildcards(join(config["sample_dir"],"{sample}_cat.fastq.gz"))

rule all:
    input:
        expand("results/renamed_results/{sample}.tsv", sample=SAMPLES)

rule build_marker:
    input:
        shortbred_identify=config["shortbred_identify_path"],
        target_fa=config["target_fa_path"],
        uniref_ref=config["uniref_path"],
        usearch_path=config["usearch_path"]
    output:
        marker="marker/colibactin_marker.fa",
        tmp_dir="results/tmp_marker"
    conda:
        "shortbred.yaml"
    threads: 10
    shell:
        '''
        python {input.shortbred_identify} --goi {input.target_fa} \
            --ref {input.uniref_ref} --markers {output.marker} --usearch {input.usearch_path} \
            --tmp {output.tmp_dir}
        '''

rule quantify_reads:
    input:
        shortbred_quantify=config["shortbred_quantify_path"],
        reads=join(config["sample_dir"], "{sample}_cat.fastq.gz"),
        usearch_path=config["usearch_path"],
        clb_marker="marker/colibactin_marker.fa",
    output:
        results="results/{sample}.tsv",
        tmpdir="results/tmp_quant/{sample}"
    conda:
        "shortbred.yaml"
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 28000
    shell:
        '''
        python {input.shortbred_quantify} --markers {input.clb_marker} \
            --wgs {input.reads} --results {output.results} --tmp {output.tmpdir} \
            --usearch {input.usearch_path} --threads {threads}
        '''

rule rename_files:
    input:
        tsv="results/{sample}.tsv"
    output:
        renamed="results/renamed_results/{sample}.tsv"
    shell:
        '''
        awk 'BEGIN{{ FS = OFS = "\t" }} {{ sample_name = "{wildcards.sample}"; sub(/\.tsv$/, "", sample_name); print (NR==1? "Sample" : sample_name), $0 }}' {input.tsv} > {output.renamed}
        '''


