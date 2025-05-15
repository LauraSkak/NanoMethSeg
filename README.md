# NanoMethSeg

NanoMethSeg is a software tool developed to analyse methylation data obtained by the Oxford Nanopore Technology long-read DNA sequencing. The aim of this software is to segment the methylome into regions of similar methylation status, which can reduce multiple testing and identify regions of interest in several different scenarios. Specific uses could be identification of transcription start sites (TSSs), imprinting control regions (ICRs) or loci undergoing X chromosome inactivation. 

### The algorithm

To segment the genome into regions of like methylation status a Hidden Markov Model (HMM) was used. A HMM is a statistical model that combines the concepts and rules of a Markov chain and Bayes rule to analyze and model a system, where the underlying states are not directly observable. In this case our observable states are the methylation data for each of the reference CpG sites and the hidden states are four types of methylation pattern;

1. Methylated regions (~100% of reads are methylated)
2. Partially methylated (~50% of reads are methylated)
3. Unmethylated (~0% of reads are methylated) 
4. Differentially methylated (~100% of reads are methylated on one haplotype (maternal/paternal) and ~0% of reads are methylated on the other haplotype).

... figure showing the four methylation pattern classes

*Complicated option*
NanoMethSeg uses the Baum-Welch algorithm for parameter optimization. The Baum-Welch algorithm is a special case of the Expectation-Maximization algorithm (EM), and is used to find optimal parameters for the HMM model. Through cycles of calculating the **e**xpected log-likelihood using the current parameter estimates, $indsæt ligning$, and updating this estimate, $indsæt ligning$ through **m**aximization of this function (EM-cycles) the algorithm converges to a set of parameter values, which maximizes the likelihood of the observed data. 

*Less complicated option*
NanoMethSeg uses the Baum-Welch algorithm for parameter optimization to improve the HMM parameters. Through cycles of parameter estimation and maximisation (EM-cycles) the parameters tangent on the set of HMM parameter values that maximises the likelihood of the observed data. 


## Tutorial

The following steps are performed on all samples, that will be included in the HMM model training, prior to running the NanoMethSeg software.

### Installation 

!!! NOT ACCURATE YET !!!

```bash
conda create -n NMS -c bioconda -c conda-forge -c HCC dorado clair3
conda activate NMS
```

Tips and tricks, before you start...

{text} will be used to indicate a path for a file, directory or model. I will use the same names throughout to ease the understanding of the flow of intermediary outputs. 

In the following decription all input files necessary for the NanoMethSeg algorithm is placed in the {NMS_input_dir}, but you could chose any file organization structure you please. The only rule which needs to be followed, to easily be able to run NanoMethSeg, is that the paths to the phased modification data and haploblock bed files have the same path name, with the exception of the sample name. For example, if you have two samples; sample1 and sample2, then the paths to the haploblock bed file could be; sample1/sample1_haploblocks.bed and sample2/sample2_haploblocks.bed. The only distinguising part of the two paths are the sample name, or the file path can be reduced to; sample_name/sample_name_haploblocks.bed. This is necessary when using more than one sample for HMM model training in NanoMethSeg.

** Input figure of the overall workflow **

### Basecalling and modification calling

For this you need;

* A dorado basecalling model
* A remora modification calling model (it has to include "5mCG" or "5mCG and 5hmCG")
* A directory which contain all pod5 files for that sample

The basecalling and modification models should match the specifications used for Oxford Nanopore long-read sequencing. For more information on how to choose the correct models go to {link to dorado}. It is recommended to use GPUs for this process, as it often takes days to finish. # Add example of how long it can take

```bash

dorado basecaller {dorado_basecalling_model} {pod5_dir} --modified-bases-models {remora_modification_model} > {basecall_modbam}

```

This should result in a modbam file containing read sequences and modification probabilities for CpG sites.


### Aligning reads

The reads are aligned using Oxford Nanopore Technology's own basecalling software dorado, which uses minimap2. It is possible to use minimap2 seperately as well, but then you in the flags have to make sure no modification data is lost, which is not the standard. Therefore, using Dorado is much more straight forward.

For this you need;
* A reference genome
* The modbam file produced in the previous step

The recommended amount of threads is 32 or 64. Depending on the size of the modbam this may take a few hours.

The alignments are sorted and indexed for vizualization (IGV is recommended). The alignments are also filtered to only include uniquely mapped reads. In this case I filter for a minimum query length of 500 and minimum mapping quality of 10. These thresholds will depend on the data provided. Use quality control to determine proper thresholds. 

```bash

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

```

### Variant calling and phasing

In the next section, we are creating the phased alignment file and haploblock bed file, which is necessary for the NMS algorithm. 

For this you need;
* A reference genome (same as the one used for alignment)
* The filtered alignment file produced in the previous step
* A clair3 model, which matches the specs of the basecalling and modification calling models

Variants are first called using Clair3 (reference), whereafter whatshap (reference) is used to phase and haplotag reads.

```bash

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

```

### Extracting modification data

The last step to be performed on each sample individually before running NMS is extraction of the phased methylation data. 

```bash

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

```

### Segmenting the genome

Once you have phased modification count data and haploblock bed files for each sample you want to include in your model, you can start running the NMS algorithm to train a HMM model and use this model to segment the genome into regions of like methylation. 

There are two main methods to indicate which samples should be included in the HMM model training; The first is to make a comma seperated list of the sample names you want to include, like sample1,sample2. As mentioned before, it is important that the file path to the phased modification count and the haploblock bed files follow the same format.

```bash

python run_segmentation_HMM.py \
    --infile_format {phased_modification_count_format} \
    --haploblock_format {haploblock_bed_format} \
    --samples {comma-seperated list of samples}

```

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

```bash

python run_segmentation_HMM.py \
    --infile_format {phased_modification_count_format} \
    --haploblock_format {haploblock_bed_format} \
    --sample_meta {sample_meta_file}

```

#### Choosing your training dataset

Tips and tricks: It is generally adviced to build seperate models for autosomes and the X chromosome. As the Y chromosome is generally monoploid, any phasing of the Y chromosome would be due to alignment artifacts. 

There are many situations in which you want to restrict the data included in the HMM model training, based on the goal of your research. If you want to research X chromosome inactivation, for example, it would be more fitting to train your model only on samples that would experience X-inactivation and limit to only the X chomosome. 

If you want to limit your training to a certain set of samples or chromosomes, you can use the flags;

* --subsample {comma-seperated list of samples found in the sample meta file}
* --chromosomes {comma-seperated list of chromosomes matching the contig names of the reference}

If you want to limit to a certain group of samples you can use specific columns of the sample meta file. This is only possible for discrete variables. 

* --subgroup {column name} {comma-seperated list of factor levels}

In the example sample meta file the flag "--group sex male" would limit the included samples to sample1.

If you are not using a sample meta file you would use the --samples instead of --subsample. In this case you would not be able to use the --subgroup flag

### Optimizing input parameters

Optional parameters include;

* --threads {thread count}
* --flank {count of flanking CpGs that should be included in calculations}
* --min_iterations {The minimum required iterations}
* --max_iterations {The maximal iteration count}
* --min_sample_coverage {The minimum number of samples to have data covering a CpG site for the CpG to be included in the observable sequence}
* --min_read_coverage {The minimum number of reads covering af CpG for the CpG to be included in the observable sequence}

The recommended specifications for running NanoMethSeg depends on the size of the input data. A large training data set covering a large part of the genome and with multiple samples way require more than 500 GB of RAM. It is also recommended to use a thread count of 16 or 32. More threads are not likely to improve running time, and can in some cases cause the running time to increase.

To limit the number of EM iterations that is run, you can set a max iteration count. In most cases the changes HMM parameter values will have plataeued after 20 to 40 iterations. You can follow the parameter value progress in the log file.

... insert picture of plot showing the progress of parameter optimization ... (Not implemented yet)

If the NanoMethSeg terminates prematurely, you can reinitiate the training from the last saved parameter values by running the same command again. If you want to re-run your data ignoring previous training you can add the flag "--reinitiate". Keep in mind that this will perminently remove any existing training data. 

The default --flank count is 0, but this can be increased to smooth out methylation data to prevent many short regions range over only one or few CpGs. It is not recommended to choose a value higher than 3 # not implemented yet. the "3" suggestion should be investigated more

If you want to limit your HMM model training and genome segmentation to regions which have a certain minimum amount of data you can set a minimum read and/or sample count for a CpG site to be included. This is generally not necessary since a lower coverage also leads to lower posterior probabilities, which can then be taken into consideration in post analysis. 


## Suggestions for uses

The data used here is the Genome in a bottle consortium and a 

!! Show an describe how the output can be used for research in genomic imprinting, X-inactivation and group-wize differential methylation !!


## FAQ

Q: How can I contribute to this software?
A: Since this software is new I would appreciate any input on how to improve the usability of NanoMethSeg. Please create a GitHub issue or write me an e-mail.


## The author (As of May 2025) 

I, Laura Skak, am a masters graduate of bioinformatics from Aarhus University and now work as a PhD student at the Department of Molecular Medicine (MOMA) in Aarhus, Denmark working on exploring the epigenetics of individuals with sex chromosome aneuploidies.