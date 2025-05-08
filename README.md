# NanoMethSeg

NanoMethSeg is a software tool developed to aid in the analysis of methylation data assertained by Nanopore long-read sequencing. The aim of the software is to segment the methylome into regions of like methylation status, which can reduce multiple testing and identify regions of interest in several different scenarios. Specific uses could be identifying potential transcription start sites (TSSs), imprinting control regions (ICRs) or loci undergoing X-inactivation. 

# The algorithm

To segment the genome into regions of like methylation status a Hidden Markov Model is used. A HMM is a statistical model that combines the concepts and rules of a Markov chain and Bayes rule to analyze and model a system, where the underlying states are not directly observable. In this case our observable states are the methylation data for each of the reference CpG sites and the hidden states are four types of methylation pattern, which are methylated regions (~100% of reads are methylated), partially methylated (~50% of reads are methylated), unmethylated (~0% of reads are methylated) and differentially methylated (~100% of reads are methylated on one haplotype (maternal/paternal) and ~0% of reads are methylated on the other haplotype).

... figure showing the four methylation pattern classes

The algorithm uses Baum-Welch parameter optimization to improve the HMM parameters. Through cycles of parameter estimation and maximisation (EM-cycles) the parameters tangent on the set of HMM parameter values the maximises the likelihood of the observed data.


# Tutorial

The following steps are performed on all samples, that will be included in the HMM model training, prior to running the NanoMethSeg software.

## Installation 

'''bash
conda create -n NMS -c bioconda -c conda-forge -c HCC dorado clair3
conda activate NMS
'''

Tips and tricks, before you start...

{text} will be used to indicate a path for a file, directory or model. I will use the same names throughout to ease the understanding of the flow of intermediary outputs. 

In the following decription all input files necessary for the NMS algorithm is placed in the {NMS_input_dir}, but you could chose any file organization structure you please. The only rule which needs to be followed, to easily be able to run NMS, is that the paths to the phased modification data and haploblock bed files have the same path name, with the exception of the sample name. For example, if you have two samples; sample1 and sample2, then the paths to the haploblock bed file could be; sample1/sample1_haploblocks.bed and sample2/sample2_haploblocks.bed. The only distinguising part of the two paths are the sample name, or the file path can be reduced to; sample_name/sample_name_haploblocks.bed. This is necessary when using more than one sample for HMM model training in NMS.

## Basecalling and modification calling

For this you need;

* A dorado basecalling model
* A remora modification calling model (it has to include "5mCG" or "5mCG and 5hmCG")
* A directory which contain all pod5 files for that sample

The basecalling and modification models should match the specs used for nanopore long-read sequencing. For more information on how to choose the correct models go to {link to dorado}. It is recommended to use GPUs for this process and it often takes days to finish.

'''bash

dorado basecaller {dorado_basecalling_model} {pod5_dir} --modified-bases-models {remora_modification_model} > {basecall_modbam}

'''

This should result in a modbam file containing read sequences and modification probabilities for CpG sites.


## Aligning reads

The reads are aligned using Oxford Nanopore Technologies own basecalling software dorado, which uses minimap2. It is possible to use minimap2 directly as well, but then you in the flags have to make sure no modification data is lost, which is not the standard. Therefore, using Dorado is much more straight forward.

For this you need;
* A reference genome
* The modbam file produced in the previous step

The recommended amount of threads is 32 or 64. Depending on the size of the modbam this may take a few hours.

The alignments are sorted and indexed for vizualization (IGV is recommended). The alignments are also filtered to only include uniquely mapped reads. In this case I filter for a minimum query length of 500 and minimum mapping quality of 10. These thresholds will depend on the data provided. Use quality control to determine proper thresholds. 

'''bash

dorado aligner -t {threads} {reference} {basecall_modbam} | samtools sort > {alignment}

samtools index -@ {threads} {alignment}

samtools view \
        --threads {threads} \
        -b -h \
        -m 500 \
        -G 1284 \
        -q 10 \
        -o {filtered_alignment} \
        {alignment}

samtools index -@ {threads} {filtered_alignment}

'''

## Variant calling and phasing

In the next section, we are creating the phased alignment file and haploblock bed file, which is necessary for the NMS algorithm. 

For this you need;
* A reference genome (same as the one used for alignment)
* The filtered alignment file produced in the previous step
* A clair3 model, which matches the specs of the basecalling and modification calling models

Variants are first called using Clair3 (reference), whereafter whatshap (reference) is used to phase and haplotag reads.

'''bash

run_clair3.sh \
        --bam_fn={filtered_alignment} \
        --ref_fn={reference} \
        --output={clair3_output_dir} \
        --threads={threads} \
        --platform="ont" \
        --model_path={clair3_model} \
        --include_all_ctgs

whatshap phase \
        --ignore-read-groups \
        --distrust-genotypes \
        --reference {reference} \
        --output {intermediary_output_dir}/{phased_vcf_prefix}.gz \
        {clair3_output_dir}/merged_output.vcf.gz {filtered_alignment}

tabix -p vcf {intermediary_output_dir}/{phased_vcf_prefix}.vcf.gz

whatshap stats --gtf {intermediary_output_dir}/{haploblock_prefix}.gft {intermediary_output_dir}/{phased_vcf_prefix}.vcf.gz

awk -F'\t' '{"{print $1, $4, $5}"}' {intermediary_output_dir}/{haploblock_prefix}.gft > {NMS_input_dir}/{haploblock_prefix}.bed

whatshap haplotag \
        --out-threads {threads} \
        --reference {reference} \
        --ignore-read-groups \
        --skip-missing-contigs \
        --tag-supplementary \
        {intermediary_output_dir}/{phased_vcf_prefix}.vcf.gz {filtered_alignment} \
    | samtools view \
        --threads {threads} \
        --bam --with-header \
    > {phased_filtered_alignment}

'''

## Extracting modification data

The last step to be performed on each sample individually before running NMS is extraction of the phased methylation data. 

'''bash

modkit pileup \
        {phased_filtered_alignment} \
        {NMS_input_dir} \
        --ref {reference} \
        --partition-tag HP \
        --prefix {modification_count_prefix} \
        --cpg \
        --only-tabs \
        --threads {threads} \
        --combine-strands \
        --combine-mods

'''

## Segmenting the genome

Once you have phased modification count data and haploblock bed files for each sample you want to include in your model, you can start running the NMS algorithm to train a HMM model and use this model to segment the genome into regions of like methylation. 

There are two main methods to indicate which samples should be included in the HMM model training; The first is to make a comma seperated list of the sample names you want to include, like sample1,sample2. As mentioned before, it is important that the file path to the phased modification count and the haploblock bed files follow the same format.

'''bash

python run_segmentation_HMM.py \
    --infile_format {phased_modification_count_format} \
    --haploblock_format {haploblock_bed_format} \
    --samples {comma-seperated list of samples}

'''

The second and recommended way to name the samples you want to include is by creating a sample meta file. This method is preferable due to its ability to carry other relevant sample information, into post analysis, like group-wise comparisons of methylation (not implemented yet). 

The sample meta file should be a tab-seperated list with the following format;

column 1: sample_ID (required, but can just be the same as sample_name)
column 2: sample_name (required)
column 3 to inf: any other covariat (optional)

There should only be one row per sample.

The sample_ID can be the same as the sample name, but can also be another identifying name for the sample. The sample_name is the name used to distinguish between outfiles.  

An example of a sample meta file could be as the following table.

| sample_ID   | sample_name | bmi         | sex         | age         |
|    :----:   |    :----:   |    :----:   |    :----:   |    :----:   |
| sample1     | sample1     | 19.5        | male        | 48          |
| sample2     | sample2     | 20.6        | female      | 32          |

If you want to use all samples in the sample meta file the following command is sufficient.

'''bash

python run_segmentation_HMM.py \
    --infile_format {phased_modification_count_format} \
    --haploblock_format {haploblock_bed_format} \
    --sample_meta {sample_meta_file}

'''

If you want to limit your training to a certain set of samples or chromosomes, you can use the flags;

* --subsample {comma-seperated list of samples found in the sample meta file}
* --chromosomes {comma-seperated list of chromosomes matching the contig names of the reference}

If you want to limit to a certain group of samples you can use specific columns of the sample meta file. This is only possible for discrete variables. 

* --subgroup {column name} {comma-seperated list of factor levels}

In the example sample meta file the flag "--group sex male" would limit the included samples to sample1