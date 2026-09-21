
s_bind_paths="$CDATA"
BOOTSTRAPS = range(1,21)
CHROM_NUM = range(1,25) # depends on the species

rule het:
  input: 
    pi= expand( "../results/{ref}/genotypes/heterozygosity/window-pi-10kb_{ref}.windowed.pi", ref=REF_GENOME),

# calculate pi for each 10kb window of the genome
rule pi:
  input:
    vcf= "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz" 
  output:
    pi= "../results/{ref}/genotypes/heterozygosity/window-pi-10kb_{ref}.windowed.pi"
  container: c_popgen
  shell:
    """
    vcftools \
      --gzvcf {input.vcf} \
      --window-pi 10000 \
      --out ../results/{wildcards.ref}/genotypes/heterozygosity/window-pi-10kb_{wildcards.ref}
    """

# detect sex scaffolds when they are not known
rule sex_scf_detection: 
    input:
        hal = "../results/{ref}/cactus/13_pinniped.hal" 
    output:
        chain = "../results/{ref}/cactus/chain/{ref}_zalcal.chain.gz",
        bc = "../results/{ref}/cactus/chain/{ref}_zalcal.bigChain.bb",
        bc_link = "../results/{ref}/cactus/chain/{ref}_zalcal.bigChain.link.bb"
    log:"logs/{ref}_sexchr_det_cactus_chain.log"
    params:
        sif = c_cactus,
        js = "../results/{ref}/cactus/scratch/sexchr/",
        local_js = "js_sexchr_{ref}",
        tmp = "../results/{ref}/cactus/scratch/sexchr/tmp",
        target = "{ref}"
    resources:
      mem_mb=20480
    shell:
        """
        readonly CACTUS_IMAGE={params.sif} 
        readonly CACTUS_SCRATCH={params.local_js} 

        # create local jobStore and tmp directories
        mkdir -p ${{CACTUS_SCRATCH}}/tmp
        mkdir -p {params.js}
        mkdir -p {params.tmp}
        mkdir -p ../results/{params.target}/cactus/chain

        apptainer exec --cleanenv \
            --overlay ${{CACTUS_SCRATCH}} \
            --bind ${{CACTUS_SCRATCH}}/tmp:/tmp,{params.js}:/run,$(pwd),{s_bind_paths} \
            --env PYTHONNOUSERSITE=1 \
            {params.sif} \
            cactus-hal2chains \
                --defaultDisk 150G \
                {params.tmp}/{params.local_js} \
                {input.hal} \
                ../results/{params.target}/cactus/chain \
                --targetGenomes {params.target} \
                --queryGenomes zalcal \
                --bigChain \
            2> {log}
        """

# remove the identified sex scaffolds
rule remove_sex:
  input:
    vcf = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz",
    scaffolds = "../data/sex_scaf_lepwed.txt"
  output:
    vcf_nosex = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}_nosexchr.vcf.gz"
  container: c_popgen
  shell:
    """
    # Build a string of multiple --not-chr arguments
    EXCLUDE=$(awk '{{printf "--not-chr %s ", $1}}' {input.scaffolds})
    vcftools --gzvcf {input.vcf} \
      $EXCLUDE \
      --recode --stdout | bgzip -c > {output.vcf_nosex}
    """

# index the resulting vcf file without sex scaffolds
rule idx:
  input:
    vcf_nosex = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}_nosexchr.vcf.gz"
  output:
    idx = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}_nosexchr.vcf.gz.tbi"
  container: c_popgen
  shell:
    """
    tabix -p vcf {input.vcf_nosex} 
    """ 

# get allele frequencies
rule vcftools_freq:
    input:
        vcf = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}_nosexchr.vcf.gz"
    output:
        freq = "../results/{ref}/freq/{ref}.frq"
    container: c_popgen
    shell: 
        """
        vcftools \
            --gzvcf {input.vcf} \
            --freq \
            --out ../results/{wildcards.ref}/freq/{wildcards.ref}
        """

# create a reduced SNP density vcf for ROH calling
rule snp_subset: 
  input: 
    vcf= "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}_nosexchr.vcf.gz"
  output: 
    vcf = "../results/{ref}/genotypes/snps/subset/variants.bi_allelic_{ref}_nosexchr_subset3kb.vcf.gz"
  container: c_popgen
  shell:
    """
    prefix=$(echo {output.vcf} | sed 's/.vcf.gz$//')
    vcftools \
      --gzvcf {input.vcf} \
      --thin 3000 \
      --recode \
      --recode-INFO-all \
      --out $prefix

    bgzip -c $prefix.recode.vcf > {output.vcf}
    rm $prefix.recode.vcf
    """

# create bed files necessary for plink
rule bed:
  input:
    vcf = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz"
  output:
    bed = "../results/{ref}/genotypes/variants.bi_allelic_{ref}"
  threads: 8
  container: c_popgen
  shell:
    """
    plink --vcf {input.vcf} --maf 0.01 --allow-extra-chr --make-bed --out {output.bed}
    """

# generate PCA with plink
rule PCA:
  input:
    bed = "../results/{ref}/genotypes/variants.bi_allelic_{ref}"
  output:
    pca = "../results/{ref}/genotypes/variants.bi_allelic_{ref}"
  threads: 8
  container: c_popgen
  shell:
    """
    plink --bfile {input.bed} --pca --allow-extra-chr --out {output.pca}
    """

# ROH analysis with bcftools
rule roh_bcf:
  input:
    vcf = "../results/{ref}/genotypes/snps/subset/variants.bi_allelic_{ref}_nosexchr_subset3kb.vcf.gz"
  output:
    roh = "../results/{ref}/ROH/bcftools/{ref}_subset3kb_roh.tsv.gz",
    roh_snps  = "../results/{ref}/ROH/bcftools/{ref}_subset3kb_roh_snps.tsv.gz"
  params:
    samples = expand( "{smp}", smp = SAMPLES )
  log:
    "logs/roh/bcftools_{ref}_subset3kb_ROH.log"
  container: c_popgen
  shell:
    """
    SAMPLES=$(echo {params.samples} | sed 's/ /,/g')

    bcftools \
      roh {input.vcf} \
      -e - \
      -s $SAMPLES \
      -O rz > {output.roh} 2> {log}

    echo -e "-----------------" >> {log} 
    bcftools \
      roh {input.vcf} \
      -e - \
      -s $SAMPLES \
      -O sz > {output.roh_snps} 2>> {log}
    """

# transform vcf into smc file
rule vcf_to_smc:
    input:
        vcf = "../results/{ref}/genotypes/snps/variants.bi_allelic_{ref}.vcf.gz"
    output:
        smc="../results/{ref}/smcpp/smc/{chrom}.smc.gz"
    params:  
        samples_smcpp=SAMPLES_smcpp
    container: c_smc
    shell:
        """
        smc++ vcf2smc \
          {input.vcf} \
          {output.smc} \
          {wildcards.chrom} \
          {wildcards.ref}:{params.samples_smcpp}
        """

# smc++ analysis for demographic reconstruction
rule smcpp_estimate:
    input:
        smc=expand("../results/{ref}/smcpp/smc/{chrom}.smc.gz",
                   ref=REF_GENOME,
                   chrom=CHROMS)
    output:
        model="../results/{ref}/smcpp/estimate/" 
    params:
        mu=0.7e-8 
    container: c_smc
    shell:
        """
        smc++ estimate \
          --cores 12 \
          -o {output.model} \
          {params.mu} \
          {input.smc}
        """

# plot smc++ results
rule smcpp_plot:
    input:
        model="../results/{ref}/smcpp/estimate/model.final.json/model.final.json" #
    output:
        png = "../results/{ref}/smcpp/plots/{ref}_ne.png"
    params:
        g=14.9 # depends on the species
    container: c_smc
    shell:
        """
        smc++ plot \
          -g {params.g} \
          -c \
          {output.png} \
          {input.model}
        """

### generate bootstraps to calculate CI
rule smcpp_bootstrap:    
    input: 
      smc=expand("../results/{ref}/smcpp/smc/{chrom}.smc.gz", ref=REF_GENOME, chrom=CHROMS) 
    output: 
      boot_chr="../results/{ref}/smcpp/bootstrap_{boot}/bootstrap_chr{chrom_num}.gz"
    params:
      nboot=20,
      chunks_per_chr=14, # depends on the species
      chunk=5000000, # depends on the species
      fake_chr=25 # depends on the species
    conda: "biopython"
    shell:
        """
        smc++ estimate \
          --cores 12 \
          -o ../results/{wildcards.ref}/smcpp/bootstrap_estimates \
          {params.mu} \
          {input} 
        """

# run smc++ for each bootstrap
rule smcpp_estimate_bootstrap:
    input:
        expand("../results/{ref}/smcpp/bootstrap_{boot}/bootstrap_chr{chrom_num}.gz",
                   ref=REF_GENOME,
                   boot=BOOTSTRAPS,
                   chrom_num=CHROM_NUM)
    output:
       "../results/{ref}/smcpp/bootstrap_estimates/{boot}/model.final.json"
    params:
        mu=0.7e-8
    container: c_smc
    shell:
        """
        smc++ estimate \
          --cores 12 \
          -o ../results/{wildcards.ref}/smcpp_36_knots6/bootstrap_estimates/{wildcards.boot} \
          {params.mu} \
          {input}
        """
# plot smc++ results      
rule smcpp_plot_bootstrap:
    input:
        "../results/{ref}/smcpp/bootstrap_estimates/{boot}/model.final.json" 
    output:
        png_boot = "../results/{ref}/smcpp/bootstrap_estimates/{boot}/{ref}_bootstraps_ne.png" 
    container: c_smc
    shell:
        """
        smc++ plot \
          -g 14.9 \
          -c \
          {output.png_boot} \
          {input}
        """