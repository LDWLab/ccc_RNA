#!/usr/bin/env python3
'''Search for twincons changes between two groups of sequences.'''
#Import required modules
import pandas as pd
import csv, os, sys, argparse
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import determine_universality
sys.path.append('/c/Users/Ishihito/Dropbox (GaTech)/Programs/Score/bin')
from PhyMeas import slice_by_name

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment', help='Sequence alignment with all species', type=str)
    parser.add_argument('twc_data', help='CSV with twc from all species', type=str)
    basepair_arg_group = parser.add_mutually_exclusive_group()
    basepair_arg_group.add_argument('-bprv','--base_pairs_rv', help='FR3D csv data in RiboVision format with base pairing information (default use ./data/contacts/PyFu_LSU_BasePairs.csv)', type=str)
    basepair_arg_group.add_argument('-bpfr3d','--base_pairs_fr3d', help='FR3D format csv with base pairing information', type=str)
    parser.add_argument('-chains','--fr3d_chains', help='Chains to use for parsing fr3d format, renumbers sequence by the order of the chains. Required with -bpfr3d', nargs='+', type=str)
    parser.add_argument('-anno','--annotation_sequence', help='Name of annotation sequence (default: LSUa_186497_PYRFU/1-3048).', default='LSUa_186497_PYRFU/1-3048', type=str)
    parser.add_argument('-tlc','--twc_low_cutoff', help='Cutoff for filtering nucleotides/bps (default: -1).', default=-1, type=float)
    parser.add_argument('-o','--output_name', help='Path for output files (default ./output/).', default='./output/', type=str)
    parser.add_argument('-traln','--generate_tr_alignment', help='Generate truncated alignment from low scoring positions.', action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def load_cons_data(csv_location):
    '''Loads twincons data into provided dict
    '''
    cons_data_dict = dict()
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                cons_data_dict[row[0]] = row[1]
            line_count+=1
    return cons_data_dict

def gap_to_nogap_construct(alignIO_out_annotated):
    '''Construct gap to nogap translation dictionary for given sequence id
    '''
    nogap_ix=0
    nogap_to_gap_ix={}
    for ix,position in enumerate(alignIO_out_annotated):
        if position == '-':
            continue
        nogap_ix+=1
        nogap_to_gap_ix[ix+1]=nogap_ix
    return(nogap_to_gap_ix)

def get_annotation_sequence_index(seq_name,alignment):
    '''Get the sequence id of provided sequence name and alignment.
    '''
    for index, sequence in enumerate(alignment):
        if seq_name in sequence.id:
            annot_seq_ix = index
            break
    return annot_seq_ix

def parse_base_pairs(csv_location):
    '''Create a list of base-pairs in the form of nucl,nucl tuples from a csv.
    '''
    base_pairs = list()
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count < 2:
                line_count+=1
                continue
            if row[3] != 'cWW':
                continue
            base_pairs.append((int(row[1]),int(row[2])))
    return base_pairs

def generate_length_adjuster(csv_file, chainList):
    output_dict = dict()
    max_chain_length = 0
    for chain in chainList:
        output_dict[chain] = max_chain_length
        temp_resi_store = list()
        csv_file.seek(0)
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count < 1:
                line_count+=1
                continue
            if (row[0].split("|")[2] != chain) or (row[1].split("|")[2] != chain):
                continue
            if row[0].split("|")[2] == chain:
                temp_resi_store.append(int(row[0].split("|")[4]))
            if row[1].split("|")[2] == chain:
                temp_resi_store.append(int(row[1].split("|")[4]))
        if len(temp_resi_store) > 0:
            max_chain_length += max(temp_resi_store)
    return output_dict

def parse_fr3d_base_pairs(csv_location, chainList):
    base_pairs = list()
    position_adjuster_by_length = dict()
    with open(csv_location, 'r+') as csv_file:
        position_adjuster_by_length = generate_length_adjuster(csv_file, chainList)
        csv_file.seek(0)
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count < 1:
                line_count+=1
                continue
            if row[2] != 'cWW':
                continue
            if (row[0].split("|")[2] not in chainList) or (row[1].split("|")[2] not in chainList):
                continue
            bp1 = int(row[0].split("|")[4])+position_adjuster_by_length[row[0].split("|")[2]]
            bp2 = int(row[1].split("|")[4])+position_adjuster_by_length[row[1].split("|")[2]]
            base_pairs.append((int(bp1),int(bp2)))
    return base_pairs


def construct_nucl_to_twc(cons_data_dict, gap_to_nogap):
    '''Returns dict with key resi number and value the twc conservation.
    '''
    nucl_to_twc = dict()
    for aln_ix in sorted(cons_data_dict.keys()):
        if gap_to_nogap.get(int(aln_ix)) is None:
            continue
        nucl_to_twc[gap_to_nogap[int(aln_ix)]] = float(cons_data_dict[aln_ix])
    return nucl_to_twc

def construct_bp_to_twc(nucl_to_twc, base_pairs):
    '''Given nucleotide to twc mapping and list of base-pairs, produces 
    base-pair to twc mapping. Both keys and value are two member tupple.
    '''
    bp_to_twc = dict()
    for base_pair in base_pairs:
        if (base_pair[0] not in nucl_to_twc.keys()) or (base_pair[1] not in nucl_to_twc.keys()):
            continue
        bp_to_twc[base_pair] = (nucl_to_twc[base_pair[0]],nucl_to_twc[base_pair[1]])
    return bp_to_twc

def filter_bps_on_twc(bp_to_twc, low_thr, high_thr):
    '''Returns 3 lists of filtered bps based on twc.
    First list is bps with both nucleotides bellow low_thr;
    Second list is bps with both nucleotides between thresholds;
    Third list is bps with both nucleotides above high_thr
    '''
    lowtwc_bps = []
    randomtwc_bps = []
    hightwc_bps = []
    for base_pair in bp_to_twc.keys():
        if (bp_to_twc[base_pair][0] <= low_thr) and (bp_to_twc[base_pair][1] <= low_thr):
            lowtwc_bps.append(base_pair)
            continue
        if (bp_to_twc[base_pair][0] >= high_thr) and (bp_to_twc[base_pair][1] >= high_thr):
            hightwc_bps.append(base_pair)
            continue
        if (low_thr < bp_to_twc[base_pair][0] < high_thr) and (low_thr < bp_to_twc[base_pair][1] < high_thr):
            randomtwc_bps.append(base_pair)
            continue
    return lowtwc_bps,randomtwc_bps, hightwc_bps

def filter_nucl_on_twc(nucl_to_twc, low_thr, high_thr):
    '''Returns 3 lists of filtered nucleotides based on twc.
    First list is nucleotides bellow low_thr;
    Second list is nucleotides between thresholds;
    Third list is nucleotides above high_thr
    '''
    lowtwc_nucl = []
    randomtwc_nucl = []
    hightwc_nucl = []
    for nucl in nucl_to_twc.keys():
        if nucl_to_twc[nucl] <= low_thr:
            lowtwc_nucl.append(nucl)
            continue
        if nucl_to_twc[nucl] >= high_thr:
            hightwc_nucl.append(nucl)
            continue
        if low_thr < nucl_to_twc[nucl] < high_thr:
            randomtwc_nucl.append(nucl)
            continue
    return lowtwc_nucl,randomtwc_nucl,hightwc_nucl

def truncate_aln(alignment_obj, index_positions, *args, **kwargs):
    aln_anchor_map = kwargs.get('aln_anchor_map', None)
    truncated_aln = alignment_obj[:,1:1]
    if aln_anchor_map is not None:
        for index in index_positions:
            truncated_aln+=alignment_obj[:,aln_anchor_map[index]-1:aln_anchor_map[index]]
    else:
        for index in index_positions:
            truncated_aln+=alignment_obj[:,index-1:index]
    return truncated_aln

def freq_iterator(column):
    '''Calculates frequency of each AA in the column.'''
    col_aalist=[]
    for aa in ['G', 'C', 'A', 'U']:
        #print(aa, column.count(aa)/len(column.replace("-", "")))
        col_aalist.append(column.count(aa)/len(column.replace("-", "")))
    #print()
    return col_aalist

def calculate_nucl_frequencies_from_bps(sliced_alns, basepairs, alngroups):
    i = 0
    bp_to_freq = dict()
    for bp in basepairs:
        bp_to_freq[bp] = dict()
        bp_to_freq[bp][alngroups[0]] = ((freq_iterator(sliced_alns[alngroups[0]][:,i]),freq_iterator(sliced_alns[alngroups[0]][:,i+1])))
        bp_to_freq[bp][alngroups[1]] = ((freq_iterator(sliced_alns[alngroups[1]][:,i]),freq_iterator(sliced_alns[alngroups[1]][:,i+1])))
        i+=2
    return bp_to_freq

def calculate_nucl_frequencies_from_nucl(sliced_alns, nucleotides, alngroups):
    i = 0
    nuc_to_freq = dict()
    for nucl in nucleotides:
        print(i, nucl)
        nuc_to_freq[nucl] = dict()
        nuc_to_freq[nucl][alngroups[0]] = freq_iterator(sliced_alns[alngroups[0]][:,i])
        nuc_to_freq[nucl][alngroups[1]] = freq_iterator(sliced_alns[alngroups[1]][:,i])
        i+=1
    return nuc_to_freq


def build_df_from_bps(low_bps, bp_to_freq, alngroups):
    gc_list=list()
    au_list=list()
    group_list=list()
    bps_list=list()
    slope_list=list()
    slope_dict=dict()

    for bp in low_bps:
        gc_sum1 = bp_to_freq[bp][alngroups[0]][0][0]+bp_to_freq[bp][alngroups[0]][0][1]+bp_to_freq[bp][alngroups[0]][1][0]+bp_to_freq[bp][alngroups[0]][1][1]
        au_sum1 = bp_to_freq[bp][alngroups[0]][0][2]+bp_to_freq[bp][alngroups[0]][0][3]+bp_to_freq[bp][alngroups[0]][1][2]+bp_to_freq[bp][alngroups[0]][1][3]
        gc_list.append(((gc_sum1*-1)/2)*100)
        au_list.append((au_sum1/2)*100)
        group_list.append(alngroups[0])

        gc_sum2 = bp_to_freq[bp][alngroups[1]][0][0]+bp_to_freq[bp][alngroups[1]][0][1]+bp_to_freq[bp][alngroups[1]][1][0]+bp_to_freq[bp][alngroups[1]][1][1]
        au_sum2 = bp_to_freq[bp][alngroups[1]][0][2]+bp_to_freq[bp][alngroups[1]][0][3]+bp_to_freq[bp][alngroups[1]][1][2]+bp_to_freq[bp][alngroups[1]][1][3]
        gc_list.append((gc_sum2/2)*100)
        au_list.append((au_sum2/2)*100)
        group_list.append(alngroups[1])
        
        bps_list.append(bp)
        bps_list.append(bp)

        try:
            slope_calc = (au_sum1/2-au_sum2/2)/(((gc_sum1*-1)/2)-gc_sum2/2)
        except ZeroDivisionError:
            slope_calc = 0
        slope_dict[bp] = slope_calc
        slope_list.append(slope_calc)
        slope_list.append(slope_calc)

    df = pd.DataFrame(
        {
        'GC': gc_list,
        'AU': au_list,
        'Group': group_list,
        'Slope': slope_list,
        'Base pairs or nucleotides': bps_list})
    return df, slope_dict

def build_df_from_nucl(low_nucs, nucl_to_freq, alngroups):
    gc_list=list()
    au_list=list()
    group_list=list()
    nucs_list=list()
    slope_list=list()
    slope_dict=dict()

    for nucl in low_nucs:
        gc_sum1 = nucl_to_freq[nucl][alngroups[0]][0]+nucl_to_freq[nucl][alngroups[0]][1]
        au_sum1 = nucl_to_freq[nucl][alngroups[0]][2]+nucl_to_freq[nucl][alngroups[0]][3]
        gc_list.append(((gc_sum1*-1))*100)
        au_list.append((au_sum1)*100)
        group_list.append(alngroups[0])

        gc_sum2 = nucl_to_freq[nucl][alngroups[1]][0]+nucl_to_freq[nucl][alngroups[1]][1]
        au_sum2 = nucl_to_freq[nucl][alngroups[1]][2]+nucl_to_freq[nucl][alngroups[1]][3]
        gc_list.append((gc_sum2)*100)
        au_list.append((au_sum2)*100)
        group_list.append(alngroups[1])
        
        nucs_list.append(nucl)
        nucs_list.append(nucl)

        try:
            slope_calc = (au_sum1-au_sum2)/(((gc_sum1*-1))-gc_sum2)
        except ZeroDivisionError:
            slope_calc = 0
        slope_dict[nucl] = slope_calc
        slope_list.append(slope_calc)
        slope_list.append(slope_calc)

    df = pd.DataFrame(
        {
        'GC': gc_list,
        'AU': au_list,
        'Group': group_list,
        'Slope': slope_list,
        'Base pairs or nucleotides': nucs_list})
    return df, slope_dict



def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    alignment = AlignIO.read(open(comm_args.alignment), "fasta")
    gap_to_nogap = gap_to_nogap_construct(alignment[get_annotation_sequence_index(comm_args.annotation_sequence,alignment)])
    nogap_to_gap = {v: k for k, v in gap_to_nogap.items()}

    alnpos_to_twc = load_cons_data(comm_args.twc_data)
    nucl_to_twc = construct_nucl_to_twc(alnpos_to_twc, gap_to_nogap)

    if comm_args.base_pairs_fr3d or comm_args.base_pairs_rv:
        if comm_args.base_pairs_fr3d:
            if not comm_args.fr3d_chains:
                raise argparse.ArgumentTypeError('When using fr3d formated csv, you must define at least one chain to parse!')
            base_pairs = parse_fr3d_base_pairs(comm_args.base_pairs_fr3d, comm_args.fr3d_chains)
        else:
            base_pairs = parse_base_pairs(comm_args.base_pairs_rv)
        
        bp_to_twc = construct_bp_to_twc(nucl_to_twc, base_pairs)
        low_bps_or_nucs, rand_bps_or_nucs, high_bps_or_nucs = filter_bps_on_twc(bp_to_twc, comm_args.twc_low_cutoff, 5)
        index_unique_bps = list(set(list(sum(low_bps_or_nucs, ()))))
        truncated_aln = truncate_aln(alignment, index_unique_bps, aln_anchor_map=nogap_to_gap)
        if comm_args.generate_tr_alignment:
            with open(comm_args.output_name+'.csv','w') as result_file:
                wr = csv.writer(result_file)
                wr.writerows(map(lambda x: [x], index_unique_bps))
            AlignIO.write(truncated_aln, comm_args.output_name+".fa", "fasta")
        sliced_alns = slice_by_name(truncated_aln)
        alngroups = list(sorted(sliced_alns.keys()))
        bp_to_freq = calculate_nucl_frequencies_from_bps(sliced_alns, low_bps_or_nucs, alngroups)
        df, slope_dict = build_df_from_bps(low_bps_or_nucs, bp_to_freq, alngroups)
    else:
        low_bps_or_nucs, rand_bps_or_nucs, high_bps_or_nucs = filter_nucl_on_twc(nucl_to_twc, comm_args.twc_low_cutoff, 5)
        #manual_nucs = [31, 55, 114, 183, 190, 451, 524, 533, 689, 792, 843, 854, 1039, 1089, 1100, 1176, 1204, 1211, 1231, 1238, 1257, 1647,1655, 1783, 1806, 1842, 1864,1888, 1934, 1939, 1956, 1963, 2128, 2133, 2153, 2156, 2179, 2461, 2519, 2540, 2727]
        #print(sorted(set(manual_nucs)-set(low_bps_or_nucs)))
        truncated_aln = truncate_aln(alignment,low_bps_or_nucs, aln_anchor_map=nogap_to_gap)
        if comm_args.generate_tr_alignment:
            with open(comm_args.output_name+'.csv','w') as result_file:
                wr = csv.writer(result_file)
                wr.writerows(map(lambda x: [x], low_bps_or_nucs))
            AlignIO.write(truncated_aln, comm_args.output_name+".fa", "fasta")
        sliced_alns = slice_by_name(truncated_aln)
        alngroups = list(sorted(sliced_alns.keys()))
        nucl_to_freq = calculate_nucl_frequencies_from_nucl(sliced_alns, low_bps_or_nucs, alngroups)
        df, slope_dict = build_df_from_nucl(low_bps_or_nucs, nucl_to_freq, alngroups)

    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(df)

    #Plotting

    sns_plot = sns.FacetGrid(data=df.head(20), hue='Base pairs or nucleotides', height=10, palette="Set2")
    sns_plot.map(plt.plot, "GC", "AU")

    sns.scatterplot(x='GC', y='AU', data=df.head(20), hue='Group', palette="Set1")

            #[X1, X2], [Y1, Y2]
    plt.plot([0, 0], [0, 100], linewidth=0.75, color='black')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(0, 100)
    plt.xlim(-100, 100)
    plt.grid()
    sns_plot.savefig(comm_args.output_name+"_".join(alngroups)+"_slopes.png")

    plt.clf()
    x = np.arange(len(slope_dict))
    y = slope_dict.values()
    #plt.grid()
    plt.bar(x, y)
    plt.xticks(x + 0.5, slope_dict.keys(), rotation='vertical')
    sns_plot.savefig(comm_args.output_name+"_".join(alngroups)+"_divergence.png")

    df.to_csv(comm_args.output_name+"_".join(alngroups)+".csv",index=False)



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))