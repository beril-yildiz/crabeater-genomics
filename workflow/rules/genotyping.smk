
rule geno:
    input: 
      #fq_fw = expand( "../results/{ref}/trimmed_fq/{sample}_1.fq.gz",ref=REF_GENOME,sample = SAMPLES),
      #fq_rv = expand( "../results/{ref}/trimmed_fq/{sample}_2.fq.gz",ref=REF_GENOME,sample = SAMPLES),
      #bam = expand( "../results/{ref}/mapped/raw/{sample}.raw.bam", ref=REF_GENOME,sample = SAMPLES),
      #raw_stats= expand( "../results/{ref}/mapped/raw/qc/{sample}.stats.txt", ref=REF_GENOME,sample = SAMPLES),
      #dedup = expand("../results/{ref}/mapped/sorted/{sample}.dedup.bam", ref=REF_GENOME,sample = SAMPLES),
      #bai = expand( "../results/{ref}/mapped/sorted/{sample}.dedup.bam.bai", ref=REF_GENOME,sample = SAMPLES),
      #filt_stats= expand( "../results/{ref}/mapped/sorted/qc/{sample}_dedup.stats.txt", ref=REF_GENOME,sample = SAMPLES),
      #vcf_filt1 = expand( "../results/{ref}/genotypes/snps/variants.snps_filtered_{ref}.vcf.gz", ref=REF_GENOME),
      #vcf_filt2 = expand( "../results/{ref}/genotypes/snps/variants.snps_masked_{ref}.vcf.gz", ref=REF_GENOME),
      #vcf_filt3 = expand( "../results/{ref}/genotypes/snps/variants.snps_miss_filter_{ref}.vcf.gz", ref=REF_GENOME),
      #vcf = expand( "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz", ref=REF_GENOME),
      #tsv = expand( "../results/{ref}/genotypes/snps/missingness_per_sample_{ref}.tsv", ref=REF_GENOME),
      #metrics = expand( "../results/{ref}/genotypes/snps/metrics/metrics.snps_raw_dp_{ref}.tsv.gz", ref=REF_GENOME),
      #depth= expand( "../results/{ref}/mapped/sorted/{sample}.dedup_depth.tsv", ref=REF_GENOME,sample = SAMPLES),
      #sum_tsv= expand( "../results/{ref}/mapped/sorted/{sample}.depth_summary.tsv", ref=REF_GENOME,sample = SAMPLES),
      "../results/synk/lobcar_karyotype.txt",
      "../results/synk/zalcal_karyotype.txt",
      "../results/synk/lobcar_zalcal"

# create a file to fix the sample order within bcf files
rule sample_order:
  output: "../data/sample_order.txt"
  params:
    sample_order = "\n".join(sorted(SAMPLES))
  shell:
    """
    echo -e "{params.sample_order}" > {output}
    """

# create a fai-index for the refence genome 
rule faidx_index:
    input:
      fa = "../data/genomes/{ref}_filtered.fna.gz",
    output:
      fai = "../data/genomes/{ref}_filtered.fna.gz.fai"
    resources:
      mem_mb=8192
    container: c_geno
    shell:
      """
      samtools faidx {input.fa}
      """

# index the refence genome for mapping with bwa
rule bwa_index:
    input:
      fa = "../data/genomes/{ref}_filtered.fna.gz",
      fai = "../data/genomes/{ref}_filtered.fna.gz.fai"
    output: "../data/genomes/{ref}_filtered.fna.gz.bwt"
    log: "logs/bwa_idx_{ref}_filtered.log"
    resources:
      mem_mb=8192
    container: c_geno
    shell:
      """
      bwa index {input.fa} &> {log}
      """

# quality-trim reads 
rule trim_raw_data:
    input:
      fq_fw = lambda wc: "../data/{ref}/" + gather_seq_center_id(wc) + "_good_1.fq.gz",
      fq_rv = lambda wc: "../data/{ref}/" + gather_seq_center_id(wc) + "_good_2.fq.gz"
    output:
      fq_fw = "../results/{ref}/trimmed_fq/{sample}_1.fq.gz",
      fq_rv = "../results/{ref}/trimmed_fq/{sample}_2.fq.gz",
      html_rep = "../results/{ref}/trimmed_fq/{sample}_report.html",
    conda: "qc_tools"
    shell:
      """
      fastp \
        -w 8 \
        -i {input.fq_fw} \
        -I {input.fq_rv} \
        -o {output.fq_fw} \
        -O {output.fq_rv} \
        -q 20 \
        -h {output.html_rep}
      """

# align reads to reference genome
rule bwa_alignment:
    input:
      ref = "../data/genomes/{ref}_filtered.fna.gz",
      ref_index = "../data/genomes/{ref}_filtered.fna.gz.bwt",
      fq_fw = "../results/{ref}/trimmed_fq/{sample}_1.fq.gz",
      fq_rv = "../results/{ref}/trimmed_fq/{sample}_2.fq.gz"
    output:
      bam = "../results/{ref}/mapped/raw/{sample}.raw.bam" 
    threads: 14
    container: c_geno
    shell:
      """
      bwa mem \
        -t {threads} \
        -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' \
        {input.ref} \
        {input.fq_fw} \
        {input.fq_rv} | \
        samtools view \
        -b \
        --threads {threads} \
        -o {output.bam} \
        -
      """

# check the statisctics of the raw bam files      
rule samtools_rawstats:
    input:
      raw_bam = "../results/{ref}/mapped/raw/{sample}.raw.bam"
    output:
      raw_stats="../results/{ref}/mapped/raw/qc/{sample}.stats.txt"
    conda: "popgen_basics"
    shell:
      """
      samtools stats {input.raw_bam} > {output.raw_stats}
      """

# sort bam file
rule samtools_rsort:
    input:
      bam = "../results/{ref}/mapped/raw/{sample}.raw.bam"
    output:
      bam = temp( "../results/{ref}/mapped/sorted/{sample}.rsorted.bam" )
    log: "logs/samtools_rsort/{ref}/{sample}.log"
    threads: 10
    conda: "popgen_basics"
    shell:
      """
      samtools sort \
        -n \
        --threads {threads} \
        -o {output.bam} \
        {input.bam}
      """

# double-check alignment with respect to the paired-end info of the reads
rule samtools_fixmate:
    input:
      bam = "../results/{ref}/mapped/sorted/{sample}.rsorted.bam"
    output:
      bam = temp( "../results/{ref}/mapped/mate_fix/{sample}.fixmate.bam" )
    log: "logs/samtools_fixmate/{ref}/{sample}.log"
    threads: 10
    conda: "popgen_basics"
    shell:
      """
      samtools fixmate \
        --threads {threads} \
        -m \
        {input.bam} \
        {output.bam}
      """

# sort mapped reads
rule samtools_sort:
    input:
      bam = "../results/{ref}/mapped/mate_fix/{sample}.fixmate.bam"
    output:
      bam = temp( "../results/{ref}/mapped/sorted/{sample}.sorted.bam" )
    log: "logs/samtools_sort/{ref}/{sample}.log"
    threads: 10
    conda: "popgen_basics"
    shell:
      """
      samtools sort \
        --threads {threads} \
        -o {output.bam} \
        {input.bam}
      """

# mark duplicates within the mapped reads
rule samtools_markdup:
    input:
      bam = "../results/{ref}/mapped/sorted/{sample}.sorted.bam"
    output:
      dedup = "../results/{ref}/mapped/sorted/{sample}.dedup.bam"
    log:
      "logs/samtools_markdup/{ref}/{sample}.log"
    threads: 10
    conda: "popgen_basics"
    shell:
      """
      samtools markdup \
        --threads {threads} \
        -r \
        {input.bam} \
        {output.dedup}
      """

# filter out the marked duplicates
rule samtools_index_rmdup:
    input:
      bam = "../results/{ref}/mapped/sorted/{sample}.dedup.bam"
    output:
      bai = "../results/{ref}/mapped/sorted/{sample}.dedup.bam.bai"
    log:
      "logs/samtools_index/{ref}/{sample}.log"
    conda: "popgen_basics"
    shell:
      """
      samtools \
        index \
        {input.bam}
      """

# check the statisctics of the filtered bam files      
rule samtools_filtstats:
    input:
      filt_bam = "../results/{ref}/mapped/sorted/{sample}.dedup.bam"
    output:
      filt_stats= "../results/{ref}/mapped/sorted/qc/{sample}_dedup.stats.txt"
    conda: "popgen_basics"
    shell:
      """
      samtools stats {input.filt_bam} > {output.filt_stats}
      """

# check depth for each mapped individual
rule individual_depth:
    input:
        bam= "../results/{ref}/mapped/sorted/{sample}.dedup.bam",
        bai= "../results/{ref}/mapped/sorted/{sample}.dedup.bam.bai"
    output:
        depth= "../results/{ref}/mapped/sorted/{sample}.dedup_depth.tsv"
    shell:
        """
        samtools depth -a {input.bam} > {output.depth}
        """

# start first step of genotyping
rule bcftools_mpileup:
    input:
        ref = "../data/genomes/{ref}_filtered.fna.gz",
        bams = expand( "../results/{ref}/mapped/sorted/{sample}.dedup.bam", ref=REF_GENOME, sample = SAMPLES ),
        bais = expand( "../results/{ref}/mapped/sorted/{sample}.dedup.bam.bai", ref=REF_GENOME, sample = SAMPLES )
    output:
        bcf = "../results/{ref}/genotypes/raw/variants.pileup_{ref}.bcf"
    params:
      pileup_par = "--min-MQ 20 --min-BQ 20 -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD -O b"
    container: c_popgen
    shell:
      """
      bcftools mpileup \
        --threads 14 \
        {params.pileup_par} \
        -f {input.ref} \
        {input.bams} \
        -o {output.bcf}
      """

# second step of genotyping
rule bcftools_call_snps:
    input:
      bcf_pileup = "../results/{ref}/genotypes/raw/variants.pileup_{ref}.bcf"
    output:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_raw_{ref}.bcf"
    params:
      call_par = "-m -a GQ -v -O b",
      all_samples = ",".join(map(str,SAMPLES))
    container: c_popgen
    shell:
      """
      bcftools call \
        {params.call_par} \
        --threads 14 \
        -o {output.bcf} \
        --samples {params.all_samples} \
        {input.bcf_pileup}
      """

# force the same sample order
rule bcftools_sample_order:
    input:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_raw_{ref}.bcf",
      txt = "../data/sample_order_{ref}.txt"
    output:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_sample_order_{ref}.bcf"
    params:
      order_par = "-O b",
    container: c_popgen
    threads: 1
    shell:
      """
      bcftools view \
        {params.order_par} \
        --threads {threads} \
        -o {output.bcf} \
        -S {input.txt} \
        {input.bcf}
      """

# for filtering: average read-depth over all samples, so it is added to the bcf 
rule bcftools_avg_dp:
    input:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_raw_{ref}.bcf"
    output:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_raw_dp_{ref}.bcf"
    container: c_popgen
    shell:
      """
      bcftools +fill-tags \
        -O b \
        {input.bcf} -o {output.bcf} -- \
        -t 'MEAN_DP:1=float(avg(FORMAT/DP))'
      """

# export the distribution of metrics we want to filter on, to visualize them for an informed threshold choice
rule raw_metrics:
    input:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_raw_dp_{ref}.bcf"
    output:
      metrics = "../results/{ref}/genotypes/snps/metrics/metrics.snps_raw_dp_{ref}.tsv.gz"
    params:
      pre = "../results/{ref}/genotypes/snps/metrics.snps_raw_dp_{ref}.tsv"
    container: c_popgen
    shell:
      """
      echo -e "chrom\tpos\tqual\tavg_depth" > {params.pre}
      bcftools query -f '%CHROM\t%POS\t%QUAL\t%MEAN_DP\n' {input.bcf} >> {params.pre}

      bgzip {params.pre}

      tabix \
        -S 1 -s 1 -b 2 -e 2 \
        {output.metrics}
      """

FILTER_PAR = "-i 'INFO/MEAN_DP > 5 & INFO/MEAN_DP < 40 & QUAL > 100 & TYPE=\"snp\"' -O z" 

# Filter the called SNPs based on quality thresholds over all samples
rule bcftools_filter_snps:
    input:
      bcf = "../results/{ref}/genotypes/snps/variants.snps_raw_dp_{ref}.bcf"
    output:
      vcf_filt1 = "../results/{ref}/genotypes/snps/variants.snps_filtered_{ref}.vcf.gz"
    container: c_popgen
    threads: 5
    params:
      filter_par = FILTER_PAR
    shell:
      """
      bcftools filter \
        {params.filter_par} \
        --threads {threads} \
        -o {output.vcf_filt1} \
        {input.bcf}
      
      tabix -p vcf {output.vcf_filt1}
      """

# define masking parameters for turning certain SNPs into missing data ./.
MASK_PAR = "-i 'FORMAT/DP<8 | FORMAT/DP>80'"

# mask individual genotypes that have made it through & get missingness info at each SNP
rule bcftools_mask_weird_depth_calls_individually:
    input:
      vcf = "../results/{ref}/genotypes/snps/variants.snps_filtered_{ref}.vcf.gz"
    output:
      vcf_filt2 = "../results/{ref}/genotypes/snps/variants.snps_masked_{ref}.vcf.gz"
    container: c_popgen
    threads: 1
    params:
      mask_par = MASK_PAR
    shell:
      """
      bcftools +setGT \
        {input.vcf} -- -t q \
        {params.mask_par} \
        -n . | \
        bcftools \
        +fill-tags /dev/stdin \
         -- -t 'FMISS=F_MISSING,NMISS=N_MISSING' | \
        bgzip > {output.vcf_filt2}
      
      tabix -p vcf {output.vcf_filt2}
      """

# if too many samples are masked at a particular SNP, filter that position.
MISSING_FILTER_PAR = "-i 'INFO/FMISS<0.33' -O z"

# filter the called SNPs based on missingness threshold
rule missingness_filter:
    input:
      vcf = "../results/{ref}/genotypes/snps/variants.snps_masked_{ref}.vcf.gz"
    output:
      vcf_filt3 = "../results/{ref}/genotypes/snps/variants.snps_miss_filter_{ref}.vcf.gz"
    container: c_popgen
    threads: 5
    params:
      filter_par = MISSING_FILTER_PAR
    shell:
      """
      bcftools filter \
        {params.filter_par} \
        --threads {threads} \
        -o {output.vcf_filt3} \
        {input.vcf}
      
      tabix -p vcf {output.vcf_filt3}
      """

# create a bi-allelic subset
rule bi_allelic_snp_filter:
    input:
      vcf = "../results/{ref}/genotypes/snps/variants.snps_miss_filter_{ref}.vcf.gz"
    output:
      vcf = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz"
    benchmark:
      "benchmark/genotyping/{ref}/snp_filtering_vcftools_{ref}.tsv"
    resources:
      mem_mb=40960
    container: c_popgen
    shell:
      """
      bcftools view \
        -m2 -M2 \
        -v snps \
        -O z \
        {input.vcf} > {output.vcf}
    
      tabix -p vcf {output.vcf}
      """

# check missingness per sample
rule missingness_per_sample:
    input:
      vcf = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz"
    output:
      tsv = "../results/{ref}/genotypes/snps/missingness_per_sample_{ref}.tsv"
    resources:
      mem_mb=20480
    container: c_popgen
    shell:
      """
      paste \
        <(bcftools query -f '[%SAMPLE\t]\n' {input.vcf} | head -1 | tr '\t' '\n') \
        <(bcftools query -f '[%GT\t]\n' {input.vcf} | awk -v OFS="\t" '{{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 }} END {{for (i in sum) print i, sum[i] / NR }}' | sort -k1,1n | cut -f 2) > {output.tsv}
      """
