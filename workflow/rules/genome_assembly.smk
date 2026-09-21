"""
snakemake -c 12 \
--use-singularity --singularity-args "--bind $CDATA" \
--use-conda \
--cluster "sbatch --export ALL -e logs/{name}.{jobid}.err -o logs/{name}.{jobid}.out --job-name {name}.{jobid}" \
--jobs 60 --latency-wait 30 geno_lobcar -n
"""

SPEC = ["lobcar"]

rule geno_lobcar:
    input: 

rule assembly:
    input:
      coverage_i = "../results/crabeater_QC/illumina_100kb_coverage.tsv",
      coverage_p = "../results/crabeater_QC/pacbio_100kb_coverage.tsv",
      busco = expand("../results/busco/{spec}", spec = SPEC),
      blob_db = expand("../results/{spec}_assembly", spec = SPEC),
      svg = expand("../results/snail/{spec}.snail.svg", spec = SPEC),
      pdf = expand("../results/snail/{spec}.snail.pdf", spec = SPEC)

rule index_assembly_fasta:
    input:
        fasta = "../results/crabeater_QC/lobcar_assembly.fa"
    output:
        fai = "../results/crabeater_QC/lobcar_assembly.fa.fai"
    container: c_popgen
    shell:
        """
        samtools faidx {input.fasta}
        """

rule map_illumina_to_assembly:
    input:
        reference = "../results/crabeater_QC/lobcar_assembly.fa",
        r1 = "../crabeaters/genome/BMK240604-CB590-ZX01-0101/BMK_DATA_20250925095607_1/Data/rawdata_Illumina_assembly/Unknown_CB590-004R0001_1.fq.gz",
        r2 = "../crabeaters/genome/BMK240604-CB590-ZX01-0101/BMK_DATA_20250925095607_1/Data/rawdata_Illumina_assembly/Unknown_CB590-004R0001_2.fq.gz"
    output:
        bam = "../results/crabeater_QC/illumina.sorted.bam",
        bai = "../results/crabeater_QC/illumina.sorted.bam.bai"
    threads: 16
    container: c_geno
    shell:
        """
        bwa mem \
            -t {threads} \
            {input.reference} \
            {input.r1} \
            {input.r2} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam}

        samtools index \
            -@ {threads} \
            {output.bam}
        """

rule map_pacbio_to_assembly:
    input:
        reference = "../results/crabeater_QC/lobcar_assembly.fa",
        reads = "../crabeaters/BMK240604-CB590-ZX01-0201/BMK_DATA_20250925103401_1/CB590-005P0001/cell/CB590-005P0001.ccs.fastq.gz"
    output:
        bam = "../results/crabeater_QC/pacbio.sorted.bam",
        bai = "../results/crabeater_QC/pacbio.sorted.bam.bai"
    threads: 24
    conda: "minimap2"
    shell:
        """
        minimap2 \
            -ax map-hifi \
            -t {threads} \
            {input.reference} \
            {input.reads} \
        | samtools sort \
            -@ {threads} \
            -o {output.bam}

        samtools index \
            -@ {threads} \
            {output.bam}
        """

rule make_coverage_windows:
    input:
        fai = "../results/crabeater_QC/lobcar_assembly.fa.fai"
    output:
        windows = "../results/crabeater_QC/genome_100kb_windows.bed"
    container: c_popgen
    shell:
        """
        cut -f1,2 {input.fai} \
        | bedtools makewindows -g - -w 100000 \
        > {output.windows}
        """

rule coverage_illumina:
    input:
        bam = "../results/crabeater_QC/illumina.sorted.bam",
        bai = "../results/crabeater_QC/illumina.sorted.bam.bai",
        bed = "../results/crabeater_QC/genome_100kb_windows.bed"
    output:
        coverage_i = "../results/crabeater_QC/illumina_100kb_coverage.tsv"
    container: c_popgen
    shell:
        """
        samtools bedcov \
            {input.bed} \
            {input.bam} \
        | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$NF/($3-$2)}}' \
        > {output}
        """

rule coverage_pacbio:
    input:
        bam = "../results/crabeater_QC/pacbio.sorted.bam",
        bai = "../results/crabeater_QC/pacbio.sorted.bam.bai",
        bed = "../results/crabeater_QC/genome_100kb_windows.bed"
    output:
        coverage_p = "../results/crabeater_QC/pacbio_100kb_coverage.tsv"
    container:
        c_popgen
    shell:
        """
        samtools bedcov \
            {input.bed} \
            {input.bam} \
        | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$NF/($3-$2)}}' \
        > {output.coverage_p}
        """


rule busco:
  input: 
    spec = "../{spec}/{spec}_filtered.fna.gz",
    busco_db = "../data/carnivora_odb10"
  output: "../results/busco/{spec}"
  log: "../results/busco/logs/{spec}.log"
  container: "$CDATA/apptainer_local/busco_v5.6.1_cv1.sif"
  shell:
    """
    busco -i {input.spec} -l {input.busco_db} --offline -o {output} -m genome -c 16 &> {log}
    """

rule blob_cr:
  input:
    fa = "../{spec}/{spec}_filtered.fa",
    busco = "../results/busco/{spec}/run_carnivora_odb10/full_table.tsv"
  output:
    blob_db = directory("../results/{spec}_assembly/")
  params: 
    dir = "../results/{spec}_assembly"
  container: "$CDATA/apptainer_local/blobtoolkit_latest.sif"
  shell:
    """
    blobtools add \
        --fasta {input.fa} \
        --busco {input.busco} \
        --create \
        {output.blob_db}
    """

rule snail_plot:
  input:
    fa = "../{spec}/{spec}_filtered.fa"
  output:
    svg = "../results/snail/{spec}.snail.svg",
    pdf = "../results/snail/{spec}.snail.pdf"
  params:
    prefix = "../results/snail/{spec}"
  container:
    "$CDATA/apptainer_local/r-snailplot.sif"
  shell:
    """
    Rscript /opt/snail_plot.R {input.fa} {params.prefix}
    """
