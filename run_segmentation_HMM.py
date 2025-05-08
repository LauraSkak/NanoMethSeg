
description = '''
MISSING
''' # FIXME

'''

conda activate pythonHMM

python SCRIPTS/run_segmentation_HMM_v4.py --infile_format RESULTS/MODIFICATION_COUNTS/sample_ID/sample_ID_T2T_masked_Y_phased_filtered_alignment_strands_and_mods_combined_modkit/sample_ID_T2T_masked_Y_phased_filtered_alignment_strands_and_mods_combined_modkit_haplotype.bed --haploblock_format RESULTS/VARIANT_CALLS/sample_ID/sample_ID_T2T_masked_Y_filtered_alignment_variant_calls/haplo_blocks.gft --outdir RESULTS/HMM_SEGMENTATION --sample_meta DATA/sample_meta.tsv --parameters parameter_values --min_iterations 2 --max_iterations 10 --stop_conditions 1

python SCRIPTS/run_segmentation_HMM_v4.py --infile_format TEST/RESULTS/MODIFICATION_COUNTS/sample_ID/sample_ID_T2T_masked_Y_phased_filtered_alignment_strands_and_mods_combined_modkit/sample_ID_T2T_masked_Y_phased_filtered_alignment_strands_and_mods_combined_modkit_haplotype.bed --haploblock_format RESULTS/VARIANT_CALLS/sample_ID/sample_ID_T2T_masked_Y_filtered_alignment_variant_calls/haplo_blocks.gft --outdir TEST/RESULTS/HMM_SEGMENTATION --sample_meta DATA/sample_meta.tsv --parameters parameter_values --min_iterations 2 --max_iterations 10 --stop_conditions 1

'''

####################################################################################################
# LOAD PACKAGES                                                                                    #
####################################################################################################

import sys, os
import argparse
import time
import numpy as np
import threading
import psutil
from concurrent.futures import thread, process
from scipy.stats import binom
import functools
from mergedeep import merge

np.seterr(divide = 'ignore') 

####################################################################################################
# ARGPARSER ARGUMENTS                                                                              #
####################################################################################################

parser = argparse.ArgumentParser(#prog = "Merge region bedMethyl data",
                                 description = description,
                                 epilog = "Good luck with using this software! Best wishes, Laura Skak")

parser.add_argument("-i","--infile_format", 
                    type = str, 
                    help = "Should be the path(s) to the phased alignment file(s). Can either be a comma seperated list of all files to be included or one file path, where the sample name is replaced by 'SAMPLE_ID'. --sample_meta is required for the latter case.") # FIXME: list option not implemented

parser.add_argument("-b", "--haploblock_bed_format", 
                    type = str, 
                    help = "Should be the path(s) to the haploblock.gtf produced by the Clair3 in the variant calling step of preprocessing. Can either be provided as a list of also be format for the paths to haploblock gtf.") # FIXME: the text is wrong

parser.add_argument("--samples", # FIXME: not implemented yet
                    type = str, 
                    default = "all",
                    help = "Comma seperated list of sample names to be included in the HMM model training set.")

parser.add_argument("-m","--sample_meta", 
                    type = str, 
                    default = "all",
                    help = "Should be the path to the sample meta file")

parser.add_argument("-s","--subsample", 
                    type = str, 
                    default = "all",
                    help = "Comma seperated list of sample_name to be included")

parser.add_argument("--chromosomes", 
                    type = str, 
                    default = False,
                    help = "A comma seperated list of chromosomes to include.")

parser.add_argument("-g","--subgroup", 
                    type = str, 
                    default = "all",
                    help = "")

parser.add_argument("-o","--outdir", 
                    type = str, 
                    help = "Should be the path to the desired outfile directory.")

parser.add_argument("-f", "--flank", 
                    type = int, 
                    default=0,
                    help = "The integer of flanking CpGs that should be included in the calculations for each CpG position. The higher the number the more smooth the state transitions are.") # FIXME: more precise description please

parser.add_argument("--min_sample_coverage", 
                    type = int, 
                    default = 1,
                    help = "Is the minimum sample count for a CpG to be added as an observation in the HMM.")

parser.add_argument("--min_read_coverage", 
                    type = int,  
                    default = 1,
                    help = "Is the minimum read count for a CpG to be added as an observation in the HMM observable sequence.")

parser.add_argument("--max_iterations", 
                    type = int, 
                    help = "Should be an upper limit to the amount of iterations")

parser.add_argument("--min_iterations", 
                    type = int, 
                    help = "Should be a lower limit to the amount of iterations")

parser.add_argument("--stop_conditions", 
                    type = float, 
                    help = "Should be a lower limit to the amount of iterations")

parser.add_argument("-t", "--threads", 
                    type = int, 
                    help = "Should be a lower limit to the amount of iterations")

parser.add_argument("-p","--parameters", 
                    type = str, 
                    default= "parameter_values",
                    help = "Should be the desired name of a parameter file or the file containing the parameter values produced in a previous run. A new one is created if it is the first run or none exist.")

parser.add_argument("--only_autosomes", # FIXME: should maybe not be a thing
                    action='store_true',
                    help = "Will exclude sex chromosomes and mitochondrial chromosome.")

parser.add_argument("--reinitiate", # FIXME: not implemented yet
                    action='store_true',
                    help = "This will remove any existing training data and restart model training from the initial parameter values.")

args = parser.parse_args()

haploblock_bed = f'{args.outdir}/HMM_haploblocks.bed'

CpG_site_lookup_dict_file_tsv = f'{args.outdir}/CpG_site_lookup'

thread_count = 32


####################################################################################################
# FUNCTIONS                                                                                        #
####################################################################################################

def create_sample_dict(sample_file):
    
    sample_dict = {}
    sample_index_dict = {}
    index_count = 0

    with open(sample_file, 'r') as file:
        lines = file.readlines()

        for line in lines[1:]:
            row = line.strip().split("\t")

            sampleid = row[1]
            
            values = [sampleid]
            
            if len(row)> 2:
                for i in range(2,len(row)):
                    samplegroup = row[i]
                    
                    values.append(samplegroup)
            
            sample_dict[sampleid] = values

            sample_index_dict[sampleid] = index_count
            index_count += 1
            
    return sample_dict, sample_index_dict


def create_sample_dict_with_subsample(sample_file, subsample_list):
    
    sample_dict = {}
    sample_index_dict = {}
    index_count = 0

    with open(sample_file, 'r') as file:
        lines = file.readlines()

        for line in lines[1:]:
            row = line.strip().split("\t")

            sampleid = row[1]
            
            values = [sampleid]
            
            if len(row)> 2:
                for i in range(2,len(row)):
                    samplegroup = row[i]
                    
                    values.append(samplegroup)
            
            if sampleid in subsample_list:

                sample_dict[sampleid] = values

                sample_index_dict[sampleid] = index_count
                index_count += 1
            
    return sample_dict, sample_index_dict


def add_sample_to_haploblock_dict(sample_name, sample_haploblock_dict):
    '''
    A sample's haploblock info is loaded into the haploblock_dictionary. 
    '''

    haplofile = sample_name.join(args.haploblock_format.split("sample_ID"))

    haploblock_count = 0

    with open(haplofile, 'r') as file:
        
        lines = file.readlines()

        for line in lines:

            row = line.strip().split("\t")

            chrom = row[0]
            block_start = int(row[3])
            block_end = int(row[4])

            if args.chromosomes:

                if chrom not in args.chromosomes.split(","):

                    continue

            if args.only_autosomes and chrom in ["chrX", "chrY", "chrM"]:

                continue

            if chrom not in sample_haploblock_dict:

                sample_haploblock_dict[chrom] = {}

            if sample_name not in sample_haploblock_dict[chrom]:

                sample_haploblock_dict[chrom][sample_name] = [] 

            sample_haploblock_dict[chrom][sample_name].append((f'{sample_name}_haploblock_{haploblock_count}', chrom, block_start, block_end))

            haploblock_count += 1
    
    return sample_haploblock_dict


def create_sample_haploblock_dict(): # Calls add_sample_to_haploblock_dict
    '''
    All sample haploblock information is added to the haploblock dictionary.
    '''

    sample_haploblock_dict = {}

    for sample_name in sample_meta_dict:

        sample_haploblock_dict = add_sample_to_haploblock_dict(sample_name, sample_haploblock_dict)

    return sample_haploblock_dict


def find_start_pos(chrom, sample_haploblock_idx):
    '''
    Finds the next starting point for a common haploblock, by finding the left-most sample haploblock.
    '''
        
    first = True
    start_pos = None
    end_pos = None
    min_sample_idx = None
    end_found = True

    for sample_name in sample_haploblock_dict[chrom]:

        # If the sample name is not in the sample_haplotype_idx, then we should look at the first 
        if sample_name not in sample_haploblock_idx:

            sample_haploblock_idx[sample_name] = 0

        # FIXME: not necessarily sorted... but in my case it is.
        haplo_idx = sample_haploblock_idx[sample_name]

        # This if-statement makes sure that there are still more haploblocks from this sample. If we are at the end of the haploblock list, then we known that all haploblocks are to the right of the current start_pos
        if haplo_idx < len(sample_haploblock_dict[chrom][sample_name]):
            
            # Here we define that samples next haploblock to check if it is more to the left than the current start pos

            sample_haploblock = sample_haploblock_dict[chrom][sample_name][haplo_idx]
            
            # If it is the first then we just define this haploblock as the left most start position
            if first:

                first = False
                start_pos = sample_haploblock[2]
                end_pos = sample_haploblock[3]
                min_sample_idx = sample_index_dict[sample_name]

            # Here we check if the haploblock we are looking at now is more to the left than the current start pos. If it is, then it is defined as the now start pos to compare other haploclocks to
            elif start_pos > sample_haploblock[2] or (start_pos == sample_haploblock[2] and end_pos < sample_haploblock[3]):

                start_pos = sample_haploblock[2]
                end_pos = sample_haploblock[3]
                min_sample_idx = sample_index_dict[sample_name]
            
            end_found = False

    return sample_haploblock_idx, start_pos, end_pos, min_sample_idx, end_found


def find_borders(chrom, end_pos, min_sample_idx, sample_haploblock_idx):
    '''
    Here we find the end position of the common haploblock by finding the right-most sample haploblock, which is still connected to the left-most sample haploblock. 
    '''

    chrom_haploblocks = [None] * len(sample_index_dict)

    end_found = False

    while not end_found:

        for sample_name in sample_haploblock_dict[chrom]:

            sample_idx = sample_index_dict[sample_name]

            if min_sample_idx == sample_idx:

                haplo_idx = sample_haploblock_idx[sample_name]

                chrom_haploblocks[sample_idx] = sample_haploblock_dict[chrom][sample_name][haplo_idx]
            
            else:
                    
                sample_end_found = False

                # Finds the last haploblock that is still within the leftmost haploblock
                while not sample_end_found:

                    haplo_idx = sample_haploblock_idx[sample_name]

                    # Makes sure we are not at the end of the haploblock list
                    if haplo_idx < len(sample_haploblock_dict[chrom][sample_name]):
                        
                        sample_haploblock = sample_haploblock_dict[chrom][sample_name][haplo_idx]
                        
                    
                    else:
                        
                        break

                    # If the start pos of this haplotype is greater than the end pos of the leftmost haplotype, then they don't overlap.
                    if end_pos > sample_haploblock[2]: 

                        chrom_haploblocks[sample_idx] = sample_haploblock

                        sample_haploblock_idx[sample_name] += 1
                    
                    # This means the next haploblock is completely to the right of the current leftmost haploblock, wherefore the whileloop should break
                    else:

                        sample_end_found = True
                
                if chrom_haploblocks[sample_idx] != None:

                    sample_haploblock_idx[sample_name] -= 1

        new_min_idx = min_sample_idx
        
        for sample_idx in range(len(chrom_haploblocks)):

            if chrom_haploblocks[sample_idx] != None:

                if end_pos < chrom_haploblocks[sample_idx][3]:

                    end_pos = chrom_haploblocks[sample_idx][3]
                    new_min_idx = sample_idx
        
        if new_min_idx == min_sample_idx:

            end_found = True
        
        else:

            min_sample_idx = new_min_idx
            
            chrom_haploblocks = [None] * len(sample_index_dict)

    return sample_haploblock_idx, end_pos, chrom_haploblocks


def write_haploblock_bed(haploblock_dict):
    '''
    Writes a bed file with all common haploblocks. This can be loaded again if the process stopped prematurely, thereby skipping the step of creating the haploblock dict.
    '''

    infile = f'{args.outdir}/HMM_haploblocks.bed'

    with open(infile, "w") as f:

        for chrom in haploblock_dict:
            
            for haploblock in haploblock_dict[chrom]:

                f.write("{}\n".format("\t".join(str(element) for element in haploblock)))
            
    return


def make_haploblocks(): # Calls find_start_pos, find_borders, and write_haploblock_bed
    '''
    Information from all sample haploblocks are aggregated to form the common haploblocks. A common haploblock is defined as any sequence of CpG positions with an overlap of at least one sample haploblock.
    '''

    haploblocks_dict = {}

    haploblock_count = 1

    for chrom in sample_haploblock_dict:

        sample_haploblock_idxs = {}
        
        sample_haploblock_idxs, start_pos, end_pos, min_sample_idx, end_found = find_start_pos(chrom, sample_haploblock_idxs)

        while not end_found:
            
            sample_haploblock_idxs, end_pos, chrom_haploblocks = find_borders(chrom, end_pos, min_sample_idx, sample_haploblock_idxs)

            if chrom not in haploblocks_dict:

                haploblocks_dict[chrom] = []

            haploblocks_dict[chrom].append((chrom, int(start_pos), int(end_pos), f'haploblock_{haploblock_count}'))
            haploblock_count += 1

            for sample_name in sample_haploblock_idxs:

                sample_idx = sample_index_dict[sample_name]
                
                if chrom_haploblocks[sample_idx] != None:

                    sample_haploblock_idxs[sample_name] += 1

            sample_haploblock_idxs, start_pos, end_pos, min_sample_idx, end_found = find_start_pos(chrom, sample_haploblock_idxs)

    write_haploblock_bed(haploblocks_dict)

    return haploblocks_dict


def site_in_haploblock(chrom, start_pos, sample_name):

    for haploblock in sample_haploblock_dict[chrom][sample_name]:

        if haploblock[2] <= start_pos and start_pos <= haploblock[3]:

            for haploblock in haploblock_dict[chrom]:

                if haploblock[1] <= start_pos and start_pos <= haploblock[2]:

                    return True, haploblock[3]
    
    return False, None


def create_subdata_dict(sample_name, haplotype):

    process_start_time = time.time()

    # print(f'Creating subdict for {sample_name} Hap{haplotype}')

    data_dict = {}

    infile = haplotype.join(sample_name.join(args.infile_format.split("sample_ID")).split("haplotype"))

    with open(infile, 'r') as file:

        lines = file.readlines()

        for line in lines:

            row = line.strip().split("\t")
            
            chrom = row[0]
            start_pos = int(row[1])
            total_count = int(row[4])
            mod_count = int(row[11])

            if args.chromosomes:

                if chrom not in args.chromosomes.split(","):

                    continue

            if args.only_autosomes and chrom in ["chrX", "chrY", "chrM"]:

                continue

            if sample_name in sample_haploblock_dict[chrom]:

                start_pos_in_haploblock, haploblock = site_in_haploblock(chrom, start_pos, sample_name)

                if start_pos_in_haploblock:

                    if chrom not in data_dict:

                        data_dict[chrom] = {}

                    if haploblock not in data_dict[chrom]:

                        data_dict[chrom][haploblock] = {}
                    
                    if start_pos not in data_dict[chrom][haploblock]:

                        data_dict[chrom][haploblock][start_pos] = {}

                    if sample_name not in data_dict[chrom][haploblock][start_pos]:

                        data_dict[chrom][haploblock][start_pos][sample_name] = {}

                    data_dict[chrom][haploblock][start_pos][sample_name][haplotype] = (int(total_count), int(mod_count))
    
    elapsed_time = time.time() - process_start_time

    # print(f'Finished creating subdict for {sample_name} Hap{haplotype} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')

    return data_dict


def create_data_dict():

    process_start_time = time.time()

    haplotypes = []
    samples = []

    for sample_name in sample_meta_dict:

        for haplotype in ["1", "2"]:

            haplotypes.append(haplotype)
            samples.append(sample_name)
    
    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            create_subdata_dict,
            samples, haplotypes
        )

    elapsed_time = time.time() - process_start_time

    print(f'Loading subdicts took {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')

    process_start_time = time.time()

    first = True

    for result in results:

        forloop_start_time = time.time()

        if first:

            data_dict = result

            first = False
        
        else:

            merge(data_dict,result)
        
        # print(f'Merging one dict took {time.strftime("%H:%M:%S", time.gmtime(time.time() - forloop_start_time))}')
        
    elapsed_time = time.time() - process_start_time

    print(f'Merging subdicts took {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')

    return data_dict


def write_site_lookup_dict_to_tsv():

    infile = f'{args.outdir}/Observable_data/CpG_site_lookup'

    with open(infile, "w") as f:

        for chrom in data_dict:

            for haploblock in data_dict[chrom]:

                for site in data_dict[chrom][haploblock]:
                    
                    f.write("{}\n".format("\t".join(str(element) for element in [chrom, haploblock, site])))
            
    return


def write_data_dict_subfile(chrom, haploblock):

    infile = f'{args.outdir}/Observable_data/{haploblock}'

    # print(infile)
    
    with open(infile, "w") as f:

        for site in data_dict[chrom][haploblock]:

            for sample in data_dict[chrom][haploblock][site]:

                for haplotype in data_dict[chrom][haploblock][site][sample]:

                    mod_data = data_dict[chrom][haploblock][site][sample][haplotype]

                    f.write("{}\n".format("\t".join(str(element) for element in [chrom, haploblock, site, sample, haplotype, mod_data[0], mod_data[1]])))


def write_data_dict_to_tsv():

    if not os.path.exists(f'{args.outdir}/Observable_data'):
        os.makedirs(f'{args.outdir}/Observable_data')

    for chrom in data_dict:

        for haploblock in data_dict[chrom]:

            infile = f'{args.outdir}/Observable_data/{haploblock}'

            # print(infile)

            with open(infile, "w") as f:

                for site in data_dict[chrom][haploblock]:

                    for sample in data_dict[chrom][haploblock][site]:

                        for haplotype in data_dict[chrom][haploblock][site][sample]:

                            mod_data = data_dict[chrom][haploblock][site][sample][haplotype]

                            f.write("{}\n".format("\t".join(str(element) for element in [chrom, haploblock, site, sample, haplotype, mod_data[0], mod_data[1]])))


def create_CpG_site_lookup():

    CpG_site_lookup_dict = {}
    sample_count_pre_filter = [0] * len(sample_meta_dict)
    sample_count_post_filter = [0] * len(sample_meta_dict)

    for chrom in data_dict:

        CpG_site_lookup_dict[chrom] = {}

        haploblocks_to_remove = []

        for haploblock in data_dict[chrom]:

            site_list = []
            sites_to_remove = []

            for site in data_dict[chrom][haploblock]:

                site_sample_count = len(data_dict[chrom][haploblock][site])

                sample_count_pre_filter[site_sample_count - 1] += 1
                samples_to_remove = []

                for sample_name in data_dict[chrom][haploblock][site]:
                    
                    # This section makes sure that the region is covered at both haplotypes.

                    if len(data_dict[chrom][haploblock][site][sample_name]) == 2:
                        
                        # This section makes sure that each haplotype has the minimum coverage

                        both_hap_has_minimum_coverage = True

                        for haplotype in data_dict[chrom][haploblock][site][sample_name]:
                            
                            if args.min_read_coverage > data_dict[chrom][haploblock][site][sample_name][haplotype][0]:
                                
                                both_hap_has_minimum_coverage = False
                        
                        if not both_hap_has_minimum_coverage:

                            samples_to_remove.append(sample_name)
                    
                    else:
                        
                        samples_to_remove.append(sample_name)

                # If a sample doesn't have coverage above the minimum in both haplotype, the CpG is excluded from the data.
                
                for sample_name in samples_to_remove:
                    
                    del data_dict[chrom][haploblock][site][sample_name]

                site_sample_count = len(data_dict[chrom][haploblock][site])
                
                if site_sample_count >= args.min_sample_coverage:

                    site_list.append(site)
                    sample_count_post_filter[site_sample_count - 1] += 1
                
                else:

                    sites_to_remove.append(site)
                    
            # Removes the sites that do not have the minimum sample coverage

            for site in sites_to_remove:
                
                del data_dict[chrom][haploblock][site]

            CpG_site_lookup_dict[chrom][haploblock] = sorted(site_list, key = lambda x_elem : int(x_elem))

            # print(haploblock, len(CpG_site_lookup_dict[chrom][haploblock]))

            if len(CpG_site_lookup_dict[chrom][haploblock]) < 2:

                haploblocks_to_remove.append(haploblock)
        
        for haploblock in haploblocks_to_remove:
                
            del data_dict[chrom][haploblock]
            del CpG_site_lookup_dict[chrom][haploblock]
    
    os.makedirs(f'{args.outdir}/Observable_data')

    write_site_lookup_dict_to_tsv()
    write_data_dict_to_tsv()

    return CpG_site_lookup_dict, sample_count_pre_filter, sample_count_post_filter


def load_haploblock_dict(infile):

    haploblocks_dict = {}

    with open(infile, 'r') as file:

        lines = file.readlines()

        for line in lines:

            row = line.strip().split("\t")
            
            chrom = row[0]
            start_pos = int(row[1])
            end_pos = int(row[2])
            haploblock = row[3]

            if chrom not in haploblocks_dict:

                haploblocks_dict[chrom] = []

            haploblocks_dict[chrom].append((chrom, start_pos, end_pos, haploblock))

    return haploblocks_dict


def load_CpG_site_lookup_table():

    infile = infile = f'{args.outdir}/Observable_data/CpG_site_lookup'

    site_lookup_dict = {}

    with open(infile, 'r') as file:

        lines = file.readlines()

        for line in lines:

            row = line.strip().split("\t")
            
            chrom = row[0]
            haploblock = row[1]
            start_pos = row[2]

            if chrom not in site_lookup_dict:

                site_lookup_dict[chrom] = {}

            if haploblock not in site_lookup_dict[chrom]:

                site_lookup_dict[chrom][haploblock] = []
            
            site_lookup_dict[chrom][haploblock].append(start_pos)

    for chrom in site_lookup_dict:
        
        for haploblock in site_lookup_dict[chrom]:
            
            site_lookup_dict[chrom][haploblock] = sorted(site_lookup_dict[chrom][haploblock], key = lambda x_elem : int(x_elem))

    return site_lookup_dict


def load_obs_subdict(infile):

    data_dict = {}

    with open(infile, 'r') as file:

        lines = file.readlines()

        for line in lines:

            row = line.strip().split("\t")

            chrom = row[0]
            haploblock = row[1]
            site = row[2]
            sample = row[3]
            haplotype = row[4]
            total_count = int(row[5])
            mod_count = int(row[6])

            if site not in data_dict:

                data_dict[site] = {}

            if sample not in data_dict[site]:

                data_dict[site][sample] = {}

            data_dict[site][sample][haplotype] = (total_count, mod_count)
    
    site_list = CpG_site_lookup_dict[chrom][haploblock]

    obs_matrix = np.zeros((2, 2, len(sample_meta_dict), len(site_list)))

    for sample_name in sample_meta_dict:

        sample_idx = sample_index_dict[sample_name]

        for i in range(len(site_list)):

            site = site_list[i]

            if sample_name in data_dict[site]:

                site_dict = data_dict[site][sample_name]

                obs_matrix[0, 0, sample_idx, i] = site_dict["1"][0]
                obs_matrix[1, 0, sample_idx, i] = site_dict["2"][0]
                obs_matrix[0, 1, sample_idx, i] = site_dict["1"][1]
                obs_matrix[1, 1, sample_idx, i] = site_dict["2"][1]

    return obs_matrix, chrom, haploblock


def load_obs_dict():

    # FIXME: baser på CpG_lookup_dict i stedet
    infiles = [f'{args.outdir}/Observable_data/{haploblock}' for haploblock in os.listdir(f'{args.outdir}/Observable_data')]

    infiles.remove(f"{args.outdir}/Observable_data/CpG_site_lookup")
    
    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            load_obs_subdict,
            infiles
        )

    obs_dict = {}

    for result in results:

        obs_submatrix, chrom, haploblock = result

        if chrom not in obs_dict:

            obs_dict[chrom] = {}
        
        obs_dict[chrom][haploblock] = obs_submatrix

    return obs_dict


def add_sample_to_obs_dict(chrom, haploblock, sample_name, site_list, obs_dict):

    sample_idx = sample_index_dict[sample_name]

    for i in range(0,len(site_list)):

        site = site_list[i]

        if sample_name in data_dict[chrom][haploblock][site]:

            site_dict = data_dict[chrom][haploblock][site][sample_name]

            obs_dict[chrom][haploblock][0, 0, sample_idx, i] = site_dict["1"][0]
            obs_dict[chrom][haploblock][1, 0, sample_idx, i] = site_dict["2"][0]
            obs_dict[chrom][haploblock][0, 1, sample_idx, i] = site_dict["1"][1]
            obs_dict[chrom][haploblock][1, 1, sample_idx, i] = site_dict["2"][1]

    return obs_dict


# FIXME: Is not parallelized
def create_obs_dict():

    obs_dict = {}

    for chrom in CpG_site_lookup_dict:

        obs_dict[chrom] = {}

        for haploblock in CpG_site_lookup_dict[chrom]:

            site_list = CpG_site_lookup_dict[chrom][haploblock]

            obs_dict[chrom][haploblock] = np.zeros((2, 2, len(sample_meta_dict), len(site_list)))

            for sample_name in sample_meta_dict:

                obs_dict = add_sample_to_obs_dict(chrom, haploblock, sample_name, site_list, obs_dict)

            del data_dict[chrom][haploblock]
    
    return obs_dict


def load_parameters(parameter_file):

    with open(parameter_file) as f:

        last_line = f.readlines()[-1].strip("\n").split("\t")
        
        transition_probs = []

        for i in range(1, 1 + n_states):
            transition_probs.append(np.fromstring(last_line[i], dtype = float, sep=','))
        
        start_iteration = int(last_line[0])

        transition_probs = np.array(transition_probs)
        emission_probs = np.array(np.fromstring(last_line[1+n_states], dtype = float, sep=','))
        start_probs = np.array(np.fromstring(last_line[2+n_states], dtype = float, sep=','))
        end_probs = np.array(np.fromstring(last_line[3+n_states], dtype = float, sep=','))

    return start_iteration, transition_probs, emission_probs, start_probs, end_probs


def logsumexp(x):
    c = x.max()
    return c + np.log(np.sum(np.exp(x-c)))


@functools.lru_cache(maxsize=1000000)
def cached_log_pmf(s,n,p):
    return binom.logpmf(s,n,p)


def cached_pmf(s,n,p):
    return np.exp(cached_log_pmf(s,n,p))


def calculate_emission_prob(Hap1, Hap2, prob_less, prob_more):

    emission_prob = 0

    for sample_name in sample_index_dict:

        sample_index = sample_index_dict[sample_name]
        
        emission_prob += cached_log_pmf(Hap1[1, sample_index], Hap1[0, sample_index], emission_probs[prob_less]) + cached_log_pmf(Hap2[1, sample_index], Hap2[0, sample_index], emission_probs[prob_more])
        
    return emission_prob


def calculate_emission_prob_S4(Hap1, Hap2, prob_less, prob_more):

    emission_prob = 0

    for sample_name in sample_index_dict:

        sample_index = sample_index_dict[sample_name]

        emission_prob += np.log(cached_pmf(Hap1[1, sample_index], Hap1[0, sample_index], emission_probs[prob_less]) * cached_pmf(Hap2[1, sample_index], Hap2[0, sample_index], emission_probs[prob_more]) * 0.5 + cached_pmf(Hap1[1, sample_index], Hap1[0, sample_index], emission_probs[prob_more]) * cached_pmf(Hap2[1, sample_index], Hap2[0, sample_index], emission_probs[prob_less]) * 0.5)

    return emission_prob


def forward_probs(obs_list, n_obs):
    
    # Makes the initial k * N table 
    alpha_table = np.zeros((n_states, n_obs))
    
    # Makes a vector that is the length of the alpha table plus 1. The last G is exclusively for calculating the total probability of the data given the parameter space.
    G = np.zeros((n_obs + 1))

    for t in range(n_obs):
        
        # Sets the observable state for the t'th position for Hap_more and Hap_less
        
        Hap1 = obs_list[0, :, :, t]
        Hap2 = obs_list[1, :, :, t]

        probs = np.zeros(n_states)
        
        for i in range(n_states):

            a = calculate_emission_prob(Hap1, Hap2, i, i)

            emission_l_prob = np.array([
                a,a,a,
                calculate_emission_prob_S4(Hap1, Hap2, i, i + 1)
                ])

            # Fills in the first column of the table
            if t == 0:

                # The probability of the first column is calculated by alpha_i(1) = pi_i * b_i(y_1).
                alpha_table[i, t] = np.log(start_probs[i]) + emission_l_prob[i]

                probs[i] = alpha_table[i, t]
            
            # Fills in all other columns of the table
            else:
            
                values = np.log(alpha_table[:, t - 1]) + np.log(transition_probs[:, i]) + emission_l_prob
                alpha_table[i, t] = logsumexp(values)
                
                probs[i] = alpha_table[i, t]
        
        # After a column of the alpha table is filled the values are scaled to equal 1 by dividing by the scaling factor G

        G[t] = logsumexp(probs)

        alpha_table[:, t] = np.exp(alpha_table[:, t] - G[t])
    

    # This fills in the last position of the G vector. Since the end probabilities are not included in the alpha table we times the last column of the alpha table with the end probabilites.

    probs = np.log(alpha_table[:, n_obs - 1]) + np.log(end_probs)
    
    G[n_obs] = logsumexp(probs)

    # This is the sum of probabilities of all paths and therefore tells how well the model fits overall.

    total_log_prob = np.sum(G)
    
    return alpha_table, total_log_prob


def backward_probs(obs_list, n_obs):
    
    # Makes the initial k * N table 
    beta_table = np.zeros((n_states, n_obs))
    
    # Makes a vector that is the length of the beta table plus 1. The last G is exclusively for calculating the total probability of the data given the parameter space.
    G = np.zeros((n_obs+1))
    
    for t in range(-1,-(n_obs+2), -1):

        # Fills the last column with only end probabilities
        if t == -1:
            
            beta_table[:, t] = end_probs

            # Since the sum of the first column is 1 and log(1) is 0, it does not contribute to the total forward probability
            G[t] = 0
        
        else:
                
            # Sets the observable state for the i'th position for Hap_less and Hap_more
            Hap1 = obs_list[0, :, :, t + 1]
            Hap2 = obs_list[1, :, :, t + 1]

            # Calculated the emission probability array used in in all calculations
            emission_l_prob = np.array([
                calculate_emission_prob(Hap1, Hap2, 0, 0),
                calculate_emission_prob(Hap1, Hap2, 1, 1),
                calculate_emission_prob(Hap1, Hap2, 2, 2),
                calculate_emission_prob_S4(Hap1, Hap2, 3, 4),
                ])
            
            # Calculated the prob that is used to fill in the last position of the G vector.
            if t == -(n_obs+1):

                probs = np.log(beta_table[:, 0]) + np.log(start_probs) + emission_l_prob

                G[t] = logsumexp(probs)
            
            # Both calculates the columns probabilities but also scales them
            else:

                for i in range(n_states):

                    values = np.log(beta_table[:, t + 1]) + np.log(transition_probs[i, :]) + emission_l_prob

                    beta_table[i, t] = logsumexp(values)
                
                probs = beta_table[:, t]
            
                G[t] = logsumexp(probs)

                beta_table[:, t] = np.exp(beta_table[:, t] - G[t])
        
    # This is the sum of probabilities of all paths and therefore tells how well the model fits overall.
    total_log_prob = np.sum(G)
    
    return beta_table, total_log_prob


# FIXME: Snak med Simon om vectorisering
def si_probs(obs_list, n_obs, alpha_table, beta_table):
    
    # We make a three-dimensional table, where there is a k times k matrix for each of the N observations.
    
    si_table = np.zeros((n_obs-1, n_states, n_states))

    for t in range(n_obs-1):
            
        # Sets the observable state for the i'th position for Hap1 and Hap2

        Hap1 = obs_list[0, :, :, t + 1]
        Hap2 = obs_list[1, :, :, t + 1]
        
        # We add op all probabilities within a k times k matrix and normalize the matrix with that sum
        
        probs = []
        
        # Fills up the k times k matrix for a particular t
        
        for i in range(n_states):
            for j in range(n_states):
                
                # alpha_i(t) * a_ij(t) * b_j(O(t+1)) * beta_i(t+1)
                
                if j <= 2:

                    emission_prob = calculate_emission_prob(Hap1, Hap2, j, j)

                else:

                    emission_prob = calculate_emission_prob_S4(Hap1, Hap2, j, j + 1)

                si_table[t, i, j] = np.log(alpha_table[i,t]) + np.log(transition_probs[i,j]) + emission_prob + np.log(beta_table[j,t+1])
                    
                probs.append(si_table[t, i, j])
        
        # Normalizes the matrix with the sum.

        s = logsumexp(np.array(probs))
        
        for i in range(n_states):
            for j in range(n_states):
                
                si_table[t, i, j] = np.exp(si_table[t, i, j] - s)
                
    return si_table


def gamma_probs(n_obs, alpha_table, beta_table):
    
    gamma_table = np.zeros((n_states, n_obs))

    for t in range(n_obs):

        # gamma_i(t) = alpha_i(t) * beta_i(t)
        gamma_table[:,t] = alpha_table[:,t] * beta_table[:,t]
        
        total_prob = np.sum(gamma_table[:,t])
        
        gamma_table[:,t] = np.divide(gamma_table[:,t], total_prob)
    
    return gamma_table


def run_alpha_table_thread(chrom, haploblock):  

    process_start_time = time.time()

    obs_list = obs_dict[chrom][haploblock]
    n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    # print(f"Started running alpha table for {haploblock}")

    alpha_table, total_forward_prob = forward_probs(obs_list, n_obs)

    process_end_time = time.time()

    # elapsed_time = process_end_time - process_start_time

    # print(f'Finished running alpha table for {haploblock} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')
    # print(
    #     f"THREAD: {threading.get_ident()}",
    #     f"PROCESS: {os.getpid()}",
    #     f"CORE_ID: {psutil.Process().cpu_num()}"
    # )

    thread_outfile_line = "{}\n".format("\t".join(str(element) for element in [psutil.Process().cpu_num(), threading.get_ident(), os.getpid(), process_start_time, process_end_time, haploblock, "alpha_table"]))

    return alpha_table, total_forward_prob, thread_outfile_line


def run_beta_table_thread(chrom, haploblock):

    process_start_time = time.time()
    
    obs_list = obs_dict[chrom][haploblock]
    n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    # print(f"Started running beta table for {haploblock}")

    beta_table, total_backward_prob = backward_probs(obs_list, n_obs)

    process_end_time = time.time()

    elapsed_time = process_end_time - process_start_time

    # print(f'Finished running beta table for {haploblock} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')
    # print(
    #     f"THREAD: {threading.get_ident()}",
    #     f"PROCESS: {os.getpid()}",
    #     f"CORE_ID: {psutil.Process().cpu_num()}"
    # )

    thread_outfile_line = "{}\n".format("\t".join(str(element) for element in [psutil.Process().cpu_num(), threading.get_ident(), os.getpid(), process_start_time, process_end_time, haploblock, "beta_table"]))

    return beta_table, total_backward_prob, thread_outfile_line


def run_si_table_thread(chrom, haploblock, alpha, beta):

    process_start_time = time.time()

    obs_list = obs_dict[chrom][haploblock]
    n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    # print(f"Started running si table for {haploblock}")

    si_table = si_probs(obs_list, n_obs, alpha, beta)

    process_end_time = time.time()

    # elapsed_time = process_end_time - process_start_time

    # print(f'Finished running si table for {haploblock} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')
    # print(
    #     f"THREAD: {threading.get_ident()}",
    #     f"PROCESS: {os.getpid()}",
    #     f"CORE_ID: {psutil.Process().cpu_num()}"
    # )

    thread_outfile_line = "{}\n".format("\t".join(str(element) for element in [psutil.Process().cpu_num(), threading.get_ident(), os.getpid(), process_start_time, process_end_time, haploblock, "si_table"]))

    return si_table, thread_outfile_line


def run_gamma_table_thread(chrom, haploblock, alpha, beta):

    process_start_time = time.time()

    n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    # print(f"Started running gamma table for {haploblock}")

    gamma_table = gamma_probs(n_obs, alpha, beta)

    process_end_time = time.time()

    # elapsed_time = process_end_time - process_start_time

    # print(f'Finished running gamma table for {haploblock} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')
    # print(
    #     f"THREAD: {threading.get_ident()}",
    #     f"PROCESS: {os.getpid()}",
    #     f"CORE_ID: {psutil.Process().cpu_num()}"
    # )

    thread_outfile_line = "{}\n".format("\t".join(str(element) for element in [psutil.Process().cpu_num(), threading.get_ident(), os.getpid(), process_start_time, process_end_time, haploblock, "gamma_table"]))

    return gamma_table, thread_outfile_line


def fill_the_HMM_tables_for_haploblock(args):

    with process.ProcessPoolExecutor(max_workers=2) as multiprocessing_executor:

        process_start_time = time.time()

        haploblock_dict = {}

        thread_outfile_line_list = []

        chrom, haploblock = args

        if first:

            alpha = multiprocessing_executor.submit(run_alpha_table_thread, chrom, haploblock)
            beta = multiprocessing_executor.submit(run_beta_table_thread, chrom, haploblock)

            haploblock_dict["alpha"], haploblock_dict["forward_prob"], thread_outfile_line = alpha.result()
            thread_outfile_line_list.append(thread_outfile_line)

            haploblock_dict["beta"], haploblock_dict["backward_prob"], thread_outfile_line = beta.result()
            thread_outfile_line_list.append(thread_outfile_line)
        
        else:

            haploblock_dict["alpha"] = table_dict[chrom][haploblock]["alpha"]
            haploblock_dict["forward_prob"] = table_dict[chrom][haploblock]["forward_prob"]
            haploblock_dict["beta"] = table_dict[chrom][haploblock]["beta"]
            haploblock_dict["backward_prob"] = table_dict[chrom][haploblock]["backward_prob"]

        si = multiprocessing_executor.submit(run_si_table_thread, chrom, haploblock, haploblock_dict["alpha"], haploblock_dict["beta"])
        gamma = multiprocessing_executor.submit(run_gamma_table_thread, chrom, haploblock, haploblock_dict["alpha"], haploblock_dict["beta"])

        haploblock_dict["si"], thread_outfile_line = si.result()
        thread_outfile_line_list.append(thread_outfile_line)
        haploblock_dict["gamma"], thread_outfile_line = gamma.result()
        thread_outfile_line_list.append(thread_outfile_line)

        process_end_time = time.time()

        elapsed_time = process_end_time - process_start_time

        # print(f'\nFinished all tables for {haploblock} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}\n')
        print(f'Finished all tables for {haploblock} in {time.strftime("%H:%M:%S", time.gmtime(elapsed_time))}')

        thread_outfile_line_list.append("{}\n".format("\t".join(str(element) for element in [psutil.Process().cpu_num(), threading.get_ident(), os.getpid(), process_start_time, process_end_time, haploblock, "filling tables"])))

    return args, haploblock_dict, thread_outfile_line_list


def fill_the_HMM_tables(haploblocks):

    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            fill_the_HMM_tables_for_haploblock,
            haploblocks
        )
        
        total_forward_prob = 0
        total_backward_prob = 0

        with open(thread_time_outfile, "a") as f:

            for result in results:

                args, haploblock_dict, line_list = result

                chrom, haploblock = args

                table_dict[chrom][haploblock] = haploblock_dict

                total_forward_prob += haploblock_dict["forward_prob"]
                total_backward_prob += haploblock_dict["backward_prob"]

                for line in line_list:
            
                    f.write(line)
    
    return total_forward_prob, total_backward_prob


def fill_alpha_table_for_haploblock(args):

    thread_outfile_line_list = []

    with process.ProcessPoolExecutor(2) as multiprocessing_executor:

        chrom, haploblock = args

        alpha = multiprocessing_executor.submit(run_alpha_table_thread, chrom, haploblock)
        beta = multiprocessing_executor.submit(run_beta_table_thread, chrom, haploblock)

        alpha_table, total_forward_prob, thread_outfile_line = alpha.result()
        thread_outfile_line_list.append(thread_outfile_line)

        beta_table, total_backward_prob, thread_outfile_line = beta.result()
        thread_outfile_line_list.append(thread_outfile_line)
    
    return args, alpha_table, total_forward_prob, beta_table, total_backward_prob, thread_outfile_line_list


def fill_alpha_tables(haploblocks):

    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            fill_alpha_table_for_haploblock,
            haploblocks
        )

        total_forward_prob = 0

        with open(thread_time_outfile, "a") as f:

            for result in results:

                args, alpha_table, forward_prob, beta_table, backward_prob, line_list = result

                chrom, haploblock = args

                table_dict[chrom][haploblock]["alpha"] = alpha_table
                table_dict[chrom][haploblock]["forward_prob"] = forward_prob
                table_dict[chrom][haploblock]["beta"] = beta_table
                table_dict[chrom][haploblock]["backward_prob"] = backward_prob

                total_forward_prob += forward_prob

                for line in line_list:

                    f.write(line)

            #process_end_time = time.time()

            #f.write("{}\n".format("\t".join(str(element) for element in [psutil.Process().cpu_num(), threading.get_ident(), os.getpid(), process_start_time, process_end_time, haploblock, "calculating new forward probabilities"])))
            
    return total_forward_prob


def write_gamma_to_outfile(gamma_outfile, region_outfile): 
    
    with open(gamma_outfile, "w") as gf, open(region_outfile, "w") as rf:

        count = 0
        
        for chrom in table_dict:

            count += 1
            region_name = f'{chrom}_region_{count}'

            for haploblock in table_dict[chrom]:

                gamma_table = table_dict[chrom][haploblock]["gamma"]
                
                max_state = 0
                CpG_count = 0
                
                for t in range(len(CpG_site_lookup_dict[chrom][haploblock])):
                    
                    start_pos = CpG_site_lookup_dict[chrom][haploblock][t]
                        
                    line_list = [chrom, start_pos, int(start_pos) + 1]

                    # state_list = []
                    
                    # for i in range(n_states):
                        
                    #     state_list.append(gamma_table[i,t])

                    state_list = gamma_table[:,t].tolist()
                    max_index = state_list.index(max(state_list)) + 1
                    
                    if max_index != max_state:

                        if t != 0:

                            region_end_pos = int(CpG_site_lookup_dict[chrom][haploblock][t-1]) + 1

                            rf.write("{}\n".format("\t".join(str(element) for element in [chrom, region_start_pos, region_end_pos, max_state, region_name, CpG_count])))

                        count += 1
                        region_name = f'{chrom}_region_{count}'
                        
                        CpG_count = 1
                        max_state = max_index
                        region_start_pos = CpG_site_lookup_dict[chrom][haploblock][t]

                    else:
                        
                        CpG_count += 1
                    
                    line_list += state_list + [region_name, max_index]

                    gf.write("{}\n".format("\t".join(str(element) for element in line_list)))

                region_end_pos = int(CpG_site_lookup_dict[chrom][haploblock][t]) + 1

                rf.write("{}\n".format("\t".join(str(element) for element in [chrom, region_start_pos, region_end_pos, max_state, region_name, CpG_count])))

    return count


def write_parameter_values_to_outfile(outfile, total_prob, count, iteration):
    
    line_list = [iteration]
    
    with open(outfile, "a") as f:
        
        for i in range(len(transition_probs)):
            
            string = np.array2string(transition_probs[:,i], separator=',', formatter={'all': lambda x: str(x)}).replace("\n ", "")
            
            line_list.append(string[1:-1])

        string = np.array2string(emission_probs, separator=',', formatter={'all': lambda x: str(x)}).replace("\n ", "")
        line_list.append(string[1:-1])
        
        string = np.array2string(start_probs, separator=',', formatter={'all': lambda x: str(x)}).replace("\n ", "")
        line_list.append(string[1:-1])
        
        string = np.array2string(end_probs, separator=',', formatter={'all': lambda x: str(x)}).replace("\n ", "")
        line_list.append(string[1:-1])
        
        line_list.append(total_prob)
        line_list.append(count)

        f.write("{}\n".format("\t".join(str(element) for element in line_list)))

    return 


def update_start_probs():

    first = True

    for chrom in table_dict:

        for haploblock in table_dict[chrom]:

            if first:

                probs = np.array([table_dict[chrom][haploblock]["gamma"][i, 0] for i in range(n_states)])

                count = 1

                first = False

            else:

                np.add(probs, np.array([table_dict[chrom][haploblock]["gamma"][i, 0] for i in range(n_states)]))

                count += 1

    return np.divide(probs, np.sum(probs))


def update_end_probs():

    first = True

    for chrom in table_dict:

        for haploblock in table_dict[chrom]:

            if first:

                probs = np.array([table_dict[chrom][haploblock]["gamma"][i, -1] for i in range(n_states)])
                
                count = 1

                first = False

            else:

                np.add(probs, np.array([table_dict[chrom][haploblock]["gamma"][i, -1] for i in range(n_states)]))

                count += 1

    return np.divide(probs, np.sum(probs))


def update_transition_probs():
    
    # a_table = np.zeros((n_states, n_states))
    
    # for i in range(n_states):
        
    #     for j in range(n_states): 
            
    #         numerator = 0
    #         denominator = 0

    #         for chrom in table_dict:

    #             for haploblock in table_dict[chrom]:

    #                 n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    #                 si_table = table_dict[chrom][haploblock]["si"]
    #                 gamma_table = table_dict[chrom][haploblock]["gamma"]
                    
    #                 for t in range(n_obs-1):
                        
    #                     numerator += si_table[t,i,j]
    #                     denominator += gamma_table[i,t]
            
    #         a_table[i,j] = numerator / denominator
    
    # return a_table

    a_table = np.zeros((n_states, n_states))

    i_list = []
    j_list = []     
    
    for i in range(n_states):
        for j in range(n_states): 

            i_list.append(i)
            j_list.append(j)

    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            update_transition_probs_matrix_position,
            i_list, j_list
        )

        for result in results:

            i, j, value = result

            a_table[i,j] = value
    
    return a_table


def update_transition_probs_matrix_position(i, j):

    chrom_list = []
    haploblock_list = []  
    i_list = []
    j_list = []    
    
    for chrom in table_dict:
        
        for haploblock in table_dict[chrom]:

            chrom_list.append(chrom)
            haploblock_list.append(haploblock)
            i_list.append(i)
            j_list.append(j)

    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            update_transition_probs_subsection,
            i_list, j_list, chrom_list, haploblock_list
        )

        total_numerators = 0
        total_denominators = 0

        for result in results:

            numerators, denominators = result

            total_numerators += numerators
            total_denominators += denominators

    total_value = np.divide(total_numerators, total_denominators)
    
    return i, j, total_value


def update_transition_probs_subsection(i, j, chrom, haploblock):

    si_table = table_dict[chrom][haploblock]["si"]
    gamma_table = table_dict[chrom][haploblock]["gamma"]
    
    numerator = np.sum(si_table[:,i,j])
    denominator = np.sum(gamma_table[i,:])

    return numerator, denominator


def update_emission_probs(): # FIXME
    
    # b_table = np.zeros(n_states + 1)

    # numerators = np.zeros(n_states + 1)
    # denominators = np.zeros(n_states + 1)

    # for chrom in table_dict:

    #     for haploblock in table_dict[chrom]:
            
    #         obs_list = obs_dict[chrom][haploblock]
    #         n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    #         gamma_table = table_dict[chrom][haploblock]["gamma"]

    #         for sample_name in sample_meta_dict:

    #             sample_idx = sample_index_dict[sample_name]

    #             for t in range(n_obs):
    #                 for i in range(n_states):
    #                     for haplotype in range(2):
                        
    #                         if i <= 2:
                                
    #                             numerators[i] += obs_list[haplotype, 1, sample_idx, t] * gamma_table[i,t] 
    #                             denominators[i] += obs_list[haplotype, 0, sample_idx, t] * gamma_table[i,t]
                            
    #                         # FIXME
    #                         # else:
                                
    #                         #     numerators[3] += Hap_less_s_list[t] * gamma_table[i,t]
    #                         #     numerators[4] += Hap_more_s_list[t] * gamma_table[i,t]
                            
    #                         #     denominators[3] += Hap_less_n_list[t] * gamma_table[i,t]
    #                         #     denominators[4] += Hap_more_n_list[t] * gamma_table[i,t]

    # numerators[3] = numerators[0]
    # numerators[4] = numerators[2]

    # denominators[3] = denominators[0]
    # denominators[4] = denominators[2]

    # for k in range(n_states + 1):
        
    #     b_table[k] = numerators[k] / denominators[k]

    # return b_table    

    chrom_list = []
    haploblock_list = []  
    
    for chrom in table_dict:
        
        for haploblock in table_dict[chrom]:

            chrom_list.append(chrom)
            haploblock_list.append(haploblock)

    b_table = np.zeros(n_states + 1)
    
    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        results = multiprocessing_executor.map(
            update_emission_probs_subsection,
            chrom_list, haploblock_list
        )

        total_numerators = np.zeros(n_states + 1)
        total_denominators = np.zeros(n_states + 1)

        for result in results:

            numerators, denominators = result

            total_numerators = np.add(total_numerators, numerators)
            total_denominators = np.add(total_denominators, denominators)
    
    b_table = np.divide(total_numerators, total_denominators)

    b_table[3] = b_table[0]
    b_table[4] = b_table[2]

    return b_table


def update_emission_probs_subsection(chrom, haploblock):

    numerators = np.zeros(n_states + 1)
    denominators = np.zeros(n_states + 1)
            
    obs_list = obs_dict[chrom][haploblock]
    n_obs = len(CpG_site_lookup_dict[chrom][haploblock])

    gamma_table = table_dict[chrom][haploblock]["gamma"]

    for sample_name in sample_meta_dict:

        sample_idx = sample_index_dict[sample_name]

        for t in range(n_obs):
            for haplotype in range(2):
                for i in range(n_states):
                
                    if i <= 2:
                        
                        numerators[i] += obs_list[haplotype, 1, sample_idx, t] * gamma_table[i,t] 
                        denominators[i] += obs_list[haplotype, 0, sample_idx, t] * gamma_table[i,t]
                    
                    # FIXME
                    # else:
                        
                    #     numerators[3] += Hap_less_s_list[t] * gamma_table[i,t]
                    #     numerators[4] += Hap_more_s_list[t] * gamma_table[i,t]
                    
                    #     denominators[3] += Hap_less_n_list[t] * gamma_table[i,t]
                    #     denominators[4] += Hap_more_n_list[t] * gamma_table[i,t]

    return numerators, denominators


def update_parameter_values():

    with process.ProcessPoolExecutor(max_workers=thread_count) as multiprocessing_executor:

        start_probs = multiprocessing_executor.submit(update_start_probs)
        end_probs = multiprocessing_executor.submit(update_end_probs)
        a_table = multiprocessing_executor.submit(update_transition_probs)
        b_table = multiprocessing_executor.submit(update_emission_probs)
        
    return start_probs.result(), end_probs.result(), a_table.result(), b_table.result()



####################################################################################################
# SET CONSTANTS                                                                                    #
####################################################################################################

states = ["methylated", "partially methylated", "unmethylated", "differentially methylated"]

n_states = len(states)

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

first = True


### Loading existing parameter values

if os.path.isfile(f'{args.outdir}/{args.parameters}'):

    start_iteration, transition_probs, emission_probs, start_probs, end_probs = load_parameters(f'{args.outdir}/{args.parameters}')

    first_run = False

else:

    start_iteration = 0
    first_run = True

    ### Generating initial probabilities ###

    # Transition probabilities, also called the a table
    transition_probs = np.array([[0.7, 0.1, 0.1, 0.1],
                                 [0.1, 0.7, 0.1, 0.1],
                                 [0.1, 0.1, 0.7, 0.1],
                                 [0.1, 0.1, 0.1, 0.7]])

    # Emission binomial distribution probabilities, also called the b table
    emission_probs = np.array([0.95, 0.5, 0.05, 0.95, 0.05])

    # End state probabilities
    end_probs = np.array([0.25, 0.25, 0.25, 0.25])

    # Beginning state probabilities
    start_probs = np.array([0.25, 0.25, 0.25, 0.25])


####################################################################################################
# THE CODE                                                                                         #
####################################################################################################

################################################
# LOADING OR CREATING THE OBSERVABLE SEQUENCES #
################################################

### Loading observation data ###

# print("\n\nLoading sample data", file=sys.stderr)

# process_start_time = time.time()

if args.subsample == "all":

    sample_meta_dict, sample_index_dict = create_sample_dict(args.sample_meta)

else:

    subsample_list = args.subsample.split(",")

    sample_meta_dict, sample_index_dict = create_sample_dict_with_subsample(args.sample_meta, subsample_list)

# elapsed_time = time.time() - process_start_time

# print("Sample data was loaded in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)

print(f"\nThis dataset consists of the {len(sample_meta_dict)} samples: {list(sample_meta_dict.keys())}", file=sys.stderr)

# If the algorithm stopped prematurely, the following will make sure the training is continued where it left off.

if os.path.exists(f'{args.outdir}/Observable_data'):

    ### Loading CpG site lookup table ###

    print("\nLoading CpG site lookup table", file=sys.stderr)

    process_start_time = time.time()
    CpG_site_lookup_dict = load_CpG_site_lookup_table()
    elapsed_time = time.time() - process_start_time

    print("CpG site lookup table was loaded in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)


    ### Loading observation dictionary ###

    print("\nLoading observation dictionary", file=sys.stderr)

    process_start_time = time.time()
    obs_dict = load_obs_dict()
    elapsed_time = time.time() - process_start_time

    print("Observation dictionary was loaded in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)


else:
    
    ### Creating a haploblock dictionary ###

    print("\nCreating a haploblock dictionary", file=sys.stderr)

    process_start_time = time.time()
    sample_haploblock_dict = create_sample_haploblock_dict()
    haploblock_dict = make_haploblocks()
    elapsed_time = time.time() - process_start_time

    print("Haploblock dictionary was created in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)


    ### Add sample modification count data to observation dictionary

    print("\nCreating a data dictionary", file=sys.stderr)

    process_start_time = time.time()
    data_dict = create_data_dict()
    elapsed_time = time.time() - process_start_time

    print("Data dictionary was created in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)


    ### Creating a CpG site lookup table ###

    print("\nCreating a CpG site lookup table", file=sys.stderr)

    process_start_time = time.time()
    CpG_site_lookup_dict, sample_count_pre_filter, sample_count_post_filter = create_CpG_site_lookup()
    elapsed_time = time.time() - process_start_time

    print("CpG site lookup table was created in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)
    
    print(f"\nThe distribution of sample coverage with the format: {list(range(1, len(sample_meta_dict) + 1))}", file=sys.stderr)
    print(f"pre filtering: {sample_count_pre_filter}\n", file=sys.stderr)
    print(f"pre filtering: {[round(float(i)/sum(sample_count_pre_filter), 3) for i in sample_count_pre_filter]}", file=sys.stderr)
    print(f"post filtering: {sample_count_post_filter}\n", file=sys.stderr)
    print(f"post filtering: {[round(float(i)/sum(sample_count_post_filter), 3) for i in sample_count_post_filter]}", file=sys.stderr)


    ### Creating an observation dictionary ###

    print("\nCreating an observation dictionary", file=sys.stderr)

    process_start_time = time.time()
    obs_dict = create_obs_dict()
    elapsed_time = time.time() - process_start_time

    print("Observation dictionary was created in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)


############################
# RUNNING HMM SEGMENTATION #
############################

thread_time_outfile = f'{args.outdir}/thread_usage_info.tsv'

if os.path.isfile(thread_time_outfile):

    with open(thread_time_outfile, "w") as f:
        
        header_list = ["core", "thread", "process", "start", "end", "haploblock", "action"]
        
        f.write("{}\n".format("\t".join(str(element) for element in header_list)))

### Makes a list of all the different inputs for threading/multiprocessing ###

table_dict = {}
args_list = []

haploblock_count = 0

for chrom in obs_dict:

    table_dict[chrom] = {}

    for haploblock in obs_dict[chrom]:

        table_dict[chrom][haploblock] = {}

        args_list.append((chrom, haploblock))

        haploblock_count += 1

print(f'\nThere are {haploblock_count} haploblock(s)', file=sys.stderr)


### Performing iterations until convergence ###

for iteration in range(start_iteration, args.max_iterations): 

    print('\n\nIteration No: ', iteration + 1, file=sys.stderr)
    print('\nStart probabilities: \n', start_probs.round(decimals=4), file=sys.stderr)
    print('\nEnd probabilities: \n', end_probs.round(decimals=4), file=sys.stderr)
    print('\nTransition probabilities:\n ', np.matrix(transition_probs.round(decimals=4)), file=sys.stderr)
    print('\nEmission probabilities: \n', emission_probs.round(decimals=4), "\n", file=sys.stderr)


    ### Filling alpha, beta, si and gamma tables ###

    print("\nFilling HMM tables", file=sys.stderr)

    process_start_time = time.time()
    total_forward_prob, total_backward_prob = fill_the_HMM_tables(args_list)
    elapsed_time = time.time() - process_start_time

    print("Filled all tables in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), "\n", file=sys.stderr)

    print(f'The total forward probability is {total_forward_prob}', file=sys.stderr)
    print(f'The total backward probability is {total_backward_prob}', file=sys.stderr)


    ### Write gamma values and regions into file ###
    
    process_start_time = time.time()

    gamma_outfile = f'{args.outdir}/gamma_table_iteration_{iteration}.bed'
    region_outfile = f'{args.outdir}/HMM_regions_iteration_{iteration}.bed'

    print("\nWriting the gamma probabilities to", gamma_outfile, file=sys.stderr)
    
    if not os.path.isfile(gamma_outfile):

        with open(gamma_outfile, "w") as f:
            
            header_list = ["iteration", "chrom", "S1_probs", "S2_probs", "S3_probs", "S4_probs", "region_name", "most_likely_state"]
            
            f.write("{}\n".format("\t".join(str(element) for element in header_list)))

    count = write_gamma_to_outfile(gamma_outfile, region_outfile)
    elapsed_time = time.time() - process_start_time
    
    print("Writing gamma probabilities to file in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)
    
    print(f"\nGrouped into {count} regions", file=sys.stderr)
    
    
    ### Writing parameter values to outfile ###

    outfile = f'{args.outdir}/{args.parameters}'
    
    if first_run:
        
        with open(outfile, "w") as f:
            
            header_list = ["iteration"]
            
            for i in range(len(transition_probs)):

                header_list.append(f'transition_probs_{i+1}')

            header_list = header_list + ["emission_probs","start_probs", "end_probs", "total_forward_prob", "region_count"]

            f.write("{}\n".format("\t".join(str(element) for element in header_list)))
        
        write_parameter_values_to_outfile(outfile, total_forward_prob, count, iteration)

        first_run = False
    
    elif iteration != start_iteration:

        write_parameter_values_to_outfile(outfile, total_forward_prob, count, iteration)
    
    
    ### Update parameters ###
    
    print("\nUpdating parameter values", file=sys.stderr) 
    
    process_start_time = time.time()
    start_probs, end_probs, a_table, b_table = update_parameter_values()
    elapsed_time = time.time() - process_start_time
    
    print("Parameters were updated in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)
    
    print('\nNew start probs:\n', start_probs.round(decimals=4), file=sys.stderr)
    print('\nNew end probs:\n', end_probs.round(decimals=4), file=sys.stderr)
    print('\nMatrix a:\n', np.matrix(a_table.round(decimals=4)), file=sys.stderr)
    print('\nMatrix b:\n', b_table.round(decimals=4), file=sys.stderr)
    
    transition_probs = a_table
    emission_probs = b_table

    print("\nProcessing the forward probabilities", file=sys.stderr) 
    
    process_start_time = time.time()
    new_total_forward_prob = fill_alpha_tables(args_list)
    elapsed_time = time.time() - process_start_time
    
    print("\nForward Probs were processed in:", time.strftime("%H:%M:%S", time.gmtime(elapsed_time)), file=sys.stderr)

    print('\nOld forward probability:', total_forward_prob, file=sys.stderr)    
    print('New forward probability:', new_total_forward_prob, file=sys.stderr)
    
    # FIXME: Snak med Simon om hvordan man gør det i log space.
    diff =  total_forward_prob - new_total_forward_prob
    
    print('\nDifference in forward probability:', diff, file=sys.stderr)

    total_forward_prob = new_total_forward_prob

    if first:

        first = False
    
    if (abs(diff) < args.stop_conditions and iteration > args.min_iterations):

        break