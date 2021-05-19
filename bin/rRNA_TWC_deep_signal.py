#!/usr/bin/env python3
'''Search for twincons changes between two groups of sequences.
Find whether sequences from a third group are closer to either one of the other two groups.
'''
#Import required modules
import csv, os, sys, argparse, importlib.util, Bio.Align, re, textwrap
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from itertools import product
import pandas as pd
import matplotlib.pyplot as plt


from twincons.TwinCons import slice_by_name
from twincons.TwinCons import main as twc_main
import twincons as TwinCons

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('alignment', help='Sequence alignment with all species', type=str)
    parser.add_argument('output_dir', help='Output directory for csvs', type=str)
    parser.add_argument('test_group', help='Test every sequence from this group.', type=str)
    parser.add_argument('-g1','--group_one', nargs='+', help='First list of named groups for twc calculation', type=str)
    parser.add_argument('-g2','--group_two', nargs='+', help='Second list of named groups for twc calculation', type=str)
    
    parser.add_argument('-tfg','--twc_filter_groups', nargs=2, help=textwrap.dedent('    Two groups for twc calculation. Each group can be made up from multiple groups.\n\
    Enter compound groups between quotes separated by | like so:\n\
    -tfg "LSUa|LSUASG" LSUe\n\
    The first group is compound from LSUa and LSUASG and the second group is LSUe'), type=str)
    
    parser.add_argument('-anno','--annotation_sequence', help='Name of annotation sequence (default: LSUa_186497_PYRFU/1-3048).', default='LSUa_186497_PYRFU/1-3048', type=str)
    parser.add_argument('-tlc','--twc_low_cutoff', help='Cutoff for filtering diverging nucleotides (default: -2).', default=-2, type=float)
    parser.add_argument('-thc','--twc_high_cutoff', help='Cutoff for similar nucleotides within diverging ones (default: 1).', default=1, type=float)
    parser.add_argument('-o','--output_name', help='Path for output files (default ./output/).', default='./output/', type=str)
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

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
    annot_seq_ix = ''
    for index, sequence in enumerate(alignment):
        if seq_name in sequence.id:
            annot_seq_ix = index
            break
    if annot_seq_ix == '':
        raise ValueError("Provided annotation sequence name "+str(seq_name)+" is not present in the provided alignment!")
    return annot_seq_ix

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
        if nucl_to_twc[nucl][0] <= low_thr:
            lowtwc_nucl.append(nucl)
            continue
        if nucl_to_twc[nucl][0] >= high_thr:
            hightwc_nucl.append(nucl)
            continue
        if low_thr < nucl_to_twc[nucl][0] < high_thr:
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

def temp_seqrecord(oldseqrecord, newname):
    temp = SeqRecord(oldseqrecord.seq, id=oldseqrecord.id, name=oldseqrecord.name, description=oldseqrecord.description)
    temp.id = newname
    return temp

def construct_aln_for_twc(filter_groups, sliced_alns):
    aln_for_calc = Bio.Align.MultipleSeqAlignment([])
    for fgroup in filter_groups:
        if re.search('\\|', fgroup):
            i = 0
            for f in fgroup.split('|'):
                if i == 0:
                    for seq in sliced_alns[f]:
                        aln_for_calc.append(seq)
                else:
                    for seq in sliced_alns[f]:
                        aln_for_calc.append(temp_seqrecord(seq, str(seq.id).replace(f, fgroup.split('|')[0])))
                i+=1
        else:
            for seq in sliced_alns[fgroup]:
                aln_for_calc.append(seq)
    return aln_for_calc

def construct_aln_for_twc2(filter_group, sliced_alns, output_aln):
    i = 0
    for fgroup in filter_group:
        if len(filter_group) > 1:
            if i == 0:
                for seq in sliced_alns[fgroup]:
                    output_aln.append(seq)
            else:
                for seq in sliced_alns[fgroup]:
                    output_aln.append(temp_seqrecord(seq, str(seq.id).replace(fgroup, filter_group[0])))
        else:
            for seq in sliced_alns[fgroup]:
                output_aln.append(seq)
        i+=1
    return output_aln


def construct_aln_for_twc_with_one_seq_vs_one_group(single_seq_and_single_group, sliced_alns, trunc_aln):
    aln_for_calc = Bio.Align.MultipleSeqAlignment([])

    i = 0
    for seq in trunc_aln:
        if seq.id == single_seq_and_single_group[0]:
            if i > 0:
                raise ValueError ("Alignment has two sequence entries with identical id:\n"+seq.id)
            aln_for_calc.append(seq)
            i+=1

    for seq in sliced_alns[single_seq_and_single_group[1]]:
        aln_for_calc.append(seq)
    return aln_for_calc

def barplot(df):
    f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
    
    df.pivot("Test sequence", "Group", "Number positions above TWC 1").plot(kind='bar',ax=ax)
    df.pivot("Test sequence", "Group", "Number positions above TWC 1").plot(kind='bar',ax=ax2)
    plt.tight_layout()
    plt.rcParams["figure.figsize"] = (20,60)
    plt.rcParams["figure.autolayout"] = True
    plt.subplots_adjust(hspace=0.08)

    ax2.get_legend().remove()
    ax.xaxis.set_visible(False)
    # zoom-in / limit the view to different portions of the data
    ax.set_ylim(150, 200)  # outliers only
    ax2.set_ylim(0, 40)  # most of the data
    ax2.yaxis.grid()
    
    # hide the spines between ax and ax2
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.tick_params(labelbottom='off')  # don't put tick labels at the bottom
    ax2.xaxis.tick_bottom()

    d = .005  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    return plt

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    alignment = AlignIO.read(open(comm_args.alignment), "fasta")
    gap_to_nogap = gap_to_nogap_construct(alignment[get_annotation_sequence_index(comm_args.annotation_sequence,alignment)])
    nogap_to_gap = {v: k for k, v in gap_to_nogap.items()}
    sliced_alns_for_calc = slice_by_name(alignment)
    
    aln_for_calc = Bio.Align.MultipleSeqAlignment([])
    aln_for_calc = construct_aln_for_twc2(comm_args.group_one, sliced_alns_for_calc, aln_for_calc)
    aln_for_calc = construct_aln_for_twc2(comm_args.group_two, sliced_alns_for_calc, aln_for_calc)
    
    #aln_for_calc = construct_aln_for_twc(comm_args.twc_filter_groups, sliced_alns_for_calc)
    
    list_for_phymeas = ['-as',aln_for_calc.format("fasta"), '-r', '-nc', '-mx', 'blastn']
    alnindex_score,sliced_alns,number_of_aligned_positions, gp_mapping = twc_main(list_for_phymeas)
    
    low_alnpos, rand_alnpos, hig_alnpos = filter_nucl_on_twc(alnindex_score, 
                        comm_args.twc_low_cutoff, comm_args.twc_high_cutoff)

    truncated_aln = truncate_aln(alignment, low_alnpos)
    
    trunc_aln_ix_to_annoseq_ix = dict()
    i = 0
    for alnpos in low_alnpos:
        i+=1
        trunc_aln_ix_to_annoseq_ix[i] = gap_to_nogap[alnpos]
        #print(alnpos, gap_to_nogap[alnpos], '\t', alnindex_score[alnpos])
    #print(truncated_aln.format("fasta"))
    #print(trunc_aln_ix_to_annoseq_ix)

    output_aln_name = '+'.join(comm_args.group_one)+'_'+'+'.join(comm_args.group_two)
    AlignIO.write(truncated_aln, comm_args.output_dir+output_aln_name+".fa", "fasta")
    ########################################################################
    sliced_trunc_aln = slice_by_name(truncated_aln)
    #Do for each sequence and aggregate data.

    testgroup_sequence_ids = list()
    for aln in truncated_aln:
        if re.search(comm_args.test_group,aln.id):
            testgroup_sequence_ids.append(aln.id)
    
    calcgroup_ids = list()
    for group in sliced_alns_for_calc.keys():
        if group != comm_args.test_group:
            calcgroup_ids.append(group)

    high_list_for_plot=list()
    i = 0
    with open(comm_args.output_dir+comm_args.test_group+".csv", mode='w') as output_csv:
        csv_writer = csv.writer(output_csv, delimiter=',')
        csv_writer.writerow(["Test sequence",
                            "Group",
                            "Number positions above "+str(comm_args.twc_high_cutoff)+" TWC", 
                            "Number positions bellow "+str(comm_args.twc_low_cutoff)+" TWC", 
                            "Number positions between",
                            "Identities of high positions by annotation seq "+comm_args.annotation_sequence])
        for group_comb in list(product(testgroup_sequence_ids, calcgroup_ids)):
            trunc_aln_for_calc = construct_aln_for_twc_with_one_seq_vs_one_group(group_comb, sliced_trunc_aln, truncated_aln)
            list_for_phymeas = ['-as',trunc_aln_for_calc.format("fasta"), '-r', '-nc', '-mx', 'blastn', '-gt', '0.9']
            alnindex_score, sliced_alns, number_of_aligned_positions, gp_mapping = twc_main(list_for_phymeas)
            low_trunc_alnpos, rand_trunc_alnpos, hig_trunc_alnpos = filter_nucl_on_twc(alnindex_score, 
                                                comm_args.twc_low_cutoff, comm_args.twc_high_cutoff)
            highpositions_by_annotation = list()
            for highposid in hig_trunc_alnpos:
                highpositions_by_annotation.append(str(trunc_aln_ix_to_annoseq_ix[highposid]))
            csv_writer.writerow([group_comb[0],
                                group_comb[1],
                                len(hig_trunc_alnpos),
                                len(low_trunc_alnpos),
                                len(rand_trunc_alnpos),
                                ";".join(highpositions_by_annotation)])
            seq_name = group_comb[0].split('/')[0]
            seq_name = seq_name.replace(comm_args.test_group+"_",'')
            if re.search(':', seq_name):
                seq_name = seq_name.split("-")[0]
            if re.search('Frey', seq_name):
                seq_name = "Freyarchaeota"
            if re.search('WNEK', seq_name):
                seq_name = "F3H4-B5_LOKI"
            if re.search('KP308697', seq_name):
                seq_name = "Uncharacterized"
            high_list_for_plot.append([seq_name,group_comb[1],len(hig_trunc_alnpos)])

    # for f in high_list_for_plot:
    #     print(f)

    df = pd.DataFrame(high_list_for_plot,columns=['Test sequence','Group','Number positions above TWC 1'])
    plt = barplot(df)
    plt.savefig(comm_args.output_dir+comm_args.test_group+".png", dpi=300)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))