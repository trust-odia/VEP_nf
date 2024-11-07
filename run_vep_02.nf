#! /appl/nextflow-23.04.1.5866/bin/nextflow
nextflow.enable.dsl=2

/*
 * pipeline input parameters
 */
params.vep_prefix = '/project/drivas_shared/projects/VEP/PMBB-Release-2020-2.0_genetic_exome_chr'
//params.chromosome_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
params.chromosome_list = ["21","22"]
params.skip = false
//params.a = "a"
//params.b = "b"
params.gnomadv4 = "/project/drivas_shared/projects/VEP_nf/gnomadv4/" 
//params.out = '/project/drivas_shared/projects/VEP_nf/out/'
params.step1 = '/project/drivas_shared/projects/VEP_nf/out/'
params.cadd_indels = "/project/drivas_shared/projects/CADD_indels/cadd/cadd_out_all/"
params.cadd_files = "/project/drivas_shared/projects/VEP_nf/CADD_files/"

process run_vep {
   publishDir "${launchDir}/out", mode:'copy'

    input:
    val chromosome
    val vep_prefix
    //path x


    output:
    path "finalAnnot.chr${chromosome}.txt"
    //path "finalAnnot.chr21.txt"
    //path "finalAnnot.chr22.txt"
    

    script:
    """
    module unload perl/5.20.2
    module load perl/5.20.2
    module load bcftools/1.17
    module load htslib/1.9
    #module load ensemblvep/110.1
    
    
    vep \
    -i ${vep_prefix}${chromosome}_NF.vcf.split.gz \
    --tab \
    --plugin CADD,snv=/project/drivas_shared/projects/VEP/caad_files_1/whole_genome_SNVs_inclAnno.tsv.gz,indels=/project/drivas_shared/projects/VEP/caad_files_1/whole_genome_SNVs_inclAnno.tsv.gz \
    --plugin LoF,loftee_path:/home/todia/loftee,human_ancestor_fa:/project/verma_shared/projects/PMBB_VEP_Annotations/human_ancestor.fa.rz,conservation_file:/project/verma_shared/projects/PMBB_VEP_Annotations/phylocsf_gerp.sql \
    --plugin dbNSFP,/project/verma_shared/projects/PMBB_VEP_Annotations/dbNSFP/dbNSFP4.0a.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred,REVEL_score \
    --plugin SpliceAI,snv=/project/verma_shared/projects/PMBB_VEP_Annotations/SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/project/verma_shared/projects/PMBB_VEP_Annotations/SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
    -everything \
    --buffer_size 1000 \
    --offline \
    --fork 4 \
    --no_escape \
    --fasta /project/drivas_shared/projects/VEP/Homo_sapiens_assembly38_nochr.fasta \
    --hgvs \
    --dir_plugins /project/verma_shared/projects/PMBB_VEP_Annotations \
    --dir_cache /project/drivas_shared/projects/  \
    -o finalAnnot.chr${chromosome}.txt

    
    """

    stub:
    """
    
    
    
    
    """
}

process add_alphamiss {
    publishDir "${launchDir}", mode:'copy'
    
    input:
    val chromosome
    //val finalAnnot
    val run_vep_ch
    
    output:
    //path add_alphamiss_ch
    path "AlphaMissense_VEP_chr${chromosome}.txt"

    script:
    """
    #! ${params["my_R"]}

    
    #load packages

    library(tidyverse)
    library(purrr)
    library(fs)


    #load in PMBB NF annotated files
    chr_by_var_cp = read_delim("/project/drivas_shared/projects/VEP_nf/out/finalAnnot.chr${chromosome}.txt", skip = 119, col_names = T)
    chr_by_var_cp = as.data.frame(chr_by_var_cp)
    chr_by_var_cp = chr_by_var_cp %>% separate(col = ${params["VEP_cols"]["uploaded_variation_col"]}, into = c("chromosome", "POS_REF_ALT"), sep = "_", remove = FALSE)

    #load in Alphamissense prediction data
    AlphaMissense_hg38_united = read_delim("/project/drivas_shared/projects/AlphaMissense/splitted/AlphaMissense_hg38_pmbb_split_chr${chromosome}.txt", col_names = T) %>% select(-${params["VEP_cols"]["chromosome_col"]})

    #join PMBB_data with Alphamissense predictions
    AlphaMissense_VEP = left_join(chr_by_var_cp, AlphaMissense_hg38_united, by = c("#Uploaded_variation"="variant_id"))
    write.table(AlphaMissense_VEP, "/project/drivas_shared/projects/VEP_nf/AlphaMissense_hg38_pmbb_splitted/AlphaMissense_VEP_chr${chromosome}.txt", col.names=T,quote = FALSE, row.names = F, sep = "\t")
    system("tee AlphaMissense_VEP_chr${chromosome}.txt")
    """
}


process merge_alphamiss_cadd{
    publishDir "${launchDir}", mode:'copy'
    //memory '32 GB'
    
    input:
    val chromosome
    //val cadd
    val  add_alphamiss_ch
    
    
    
    

    output:
    path "VEP_AM_CADD_chr${chromosome}.txt"
    //tuple  val(chromosome), path("VEP_AM_CADD_chr${chromosome}.txt")

    script:
    """

    #! ${params["my_R"]}

    
    #load packages

    library(tidyverse)
    library(purrr)
    library(fs)


    #load CADD indel annotations
    cadd_files_mgd_united = read_delim("/project/drivas_shared/projects/VEP_nf/CADD_files/CADD_chr${chromosome}.txt") %>% as.data.frame()

    #load previous PMBB VEP annotation data
    AlphaMissense_VEP = read_delim("/project/drivas_shared/projects/VEP_nf/AlphaMissense_hg38_pmbb_splitted/AlphaMissense_VEP_chr${chromosome}.txt")

    #join the previous VEP data with CADD indel scores
    PMBB_CADD_indels_full = left_join(AlphaMissense_VEP, cadd_files_mgd_united %>% 
    select(${params["CADD_cols"]["uploaded_variation_col"]},${params["CADD_cols"]["RawScore_col"]},${params["CADD_cols"]["PHRED_col"]}), by = c("#Uploaded_variation"))

    #change datatype to character (RawScore)
    PMBB_CADD_indels_full = PMBB_CADD_indels_full %>% mutate(RawScore = as.character(${params["CADD_cols"]["RawScore_col"]}))

    #fill in the missing values with cadd and PHRED scores
    PMBB_CADD_indels_full_filled = PMBB_CADD_indels_full %>% mutate(across(c(${params["VEP_cols"]["CADD_RAW_col"]}), na_if, "-")) %>% 
    mutate(CADD_RAW_RawScore = coalesce(${params["VEP_cols"]["CADD_RAW_col"]},${params["CADD_cols"]["RawScore_col"]}))

    #change datatype to character (PHRED)
    PMBB_CADD_indels_full_filled = PMBB_CADD_indels_full_filled %>% mutate(PHRED = as.character(${params["CADD_cols"]["PHRED_col"]}))

    PMBB_CADD_indels_full_filled = PMBB_CADD_indels_full_filled %>% mutate(across(c(${params["VEP_cols"]["CADD_PHRED_col"]}), na_if, "-")) %>% 
    mutate(CADD_PHRED_merge = coalesce(${params["VEP_cols"]["CADD_PHRED_col"]},${params["CADD_cols"]["PHRED_col"]}))

    #save result
    write.table(PMBB_CADD_indels_full_filled, "/project/drivas_shared/projects/VEP_nf/VEP_AM_CADD/VEP_AM_CADD_chr${chromosome}.txt", col.names=T,quote = FALSE, row.names = F, sep = "\t")
    system("tee VEP_AM_CADD_chr${chromosome}.txt")
    """

    stub:
    """
    
    """
}


process run_vep_gnomadv4{
    publishDir "${launchDir}/gnomadv4", mode:'copy'
    //memory '32 GB'
    input:
    val vep_prefix 
    val chromosome
    val merge_alphamiss_cadd_ch
    


    output:
    path "finalAnnot_gmv4_chr${chromosome}.txt"
    //tuple  val(chromosome), path("finalAnnot_gmv4_chr${chromosome}.txt")


    script:
    """
    module unload perl/5.20.2
    module load perl/5.20.2
    module load bcftools/1.17
    module load htslib/1.9
    #module load ensemblvep/110.1
    
    vep \
    -i ${vep_prefix}${chromosome}_NF.vcf.split.gz \
    --tab \
    --offline \
    --dir_cache /project/drivas_shared/projects/ \
    --custom file=/project/drivas_shared/projects/gnomad4/gnomad4.0_GRCh38_combined_af.vcf.bgz,short_name=gnomad4,format=vcf,type=exact,fields=AF%AF_afr%AF_ami%AF_amr%AF_asj%AF_eas%AF_fin%AF_nfe%AF_mid%AF_sas%AF_remaining \
    -o finalAnnot_gmv4_chr${chromosome}.txt
    
    
    """

    stub:
    """
    
    
    """


}
process gnomadv4_freq{
    publishDir "${launchDir}", mode:'copy'

    input:
    val chromosome
    
    val  run_vep_gnomadv4_ch
    
    output:
    path "VEP_AM_CADD_GM4_chr${chromosome}.txt"

    script:
    """

    #! ${params["my_R"]}

    #load packages
    library(tidyverse)
    library(purrr)
    library(fs)


    gnomad4_chr = read_delim("/project/drivas_shared/projects/VEP_nf/gnomadv4/finalAnnot_gmv4_chr${chromosome}.txt", skip = 53, col_names = T)
    gnomad4_chr = as.data.frame(gnomad4_chr)
    gnomad4_chr = gnomad4_chr %>% separate(col = ${params["VEP_cols"]["uploaded_variation_col"]}, into = c("chromosome", "POS_REF_ALT"), sep = "_", remove = FALSE)
    
    PMBB_CADD_indels_full_filled = read_delim("/project/drivas_shared/projects/VEP_nf/VEP_AM_CADD/VEP_AM_CADD_chr${chromosome}.txt")
    PMBB_CADD_indels_full_filled_gmd4 = inner_join(PMBB_CADD_indels_full_filled, gnomad4_chr %>% 
    select(${params["VEP_cols"]["uploaded_variation_col"]}, ${params["VEP_cols"]["gnomad4_col"]}, ${params["VEP_cols"]["gnomad4_AF_col"]}, ${params["VEP_cols"]["gnomad4_AF_afr_col"]}, ${params["VEP_cols"]["gnomad4_AF_ami_col"]}, ${params["VEP_cols"]["gnomad4_AF_amr_col"]}, ${params["VEP_cols"]["gnomad4_AF_asj_col"]}, ${params["VEP_cols"]["gnomad4_AF_eas_col"]}, ${params["VEP_cols"]["gnomad4_AF_fin_col"]}, ${params["VEP_cols"]["gnomad4_AF_nfe_col"]}, ${params["VEP_cols"]["gnomad4_AF_mid_col"]}, ${params["VEP_cols"]["gnomad4_AF_sas_col"]}, ${params["VEP_cols"]["gnomad4_AF_remaining_col"]}) %>% distinct())
    write.table(PMBB_CADD_indels_full_filled_gmd4, "/project/drivas_shared/projects/VEP_nf/VEP_AM_CADD_GM4/VEP_AM_CADD_GM4_chr${chromosome}.txt", col.names=T,quote = FALSE, row.names = F, sep = "\t")
    system("tee VEP_AM_CADD_GM4_chr${chromosome}.txt")
    """

    stub:
    """
    
    """
}




workflow {
    chromosome = Channel.fromList(params.chromosome_list)
    vep_prefix = params.vep_prefix
    finalAnnot = params.step1
    cadd_indels = params.cadd_indels
    gnomadv4 = params.gnomadv4
    gnomadv4_a = params.gnomadv4
    cadd = params.cadd_files  
    
    
    
    
    run_vep_ch = run_vep(chromosome, vep_prefix)
    add_alphamiss_ch = add_alphamiss(chromosome,run_vep_ch)
    merge_alphamiss_cadd_ch = merge_alphamiss_cadd(chromosome, add_alphamiss_ch)
    run_vep_gnomadv4_ch = run_vep_gnomadv4(vep_prefix, chromosome ,merge_alphamiss_cadd_ch)
    gnomadv4_freq_ch = gnomadv4_freq(chromosome, run_vep_gnomadv4_ch)
    
    
}

