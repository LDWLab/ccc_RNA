#!/usr/bin/env python3

#%%
import pandas as pd
import csv
from Bio import AlignIO

#%%
def load_cons_data(csv_location,cons_data_dict):
    '''Loads twincons data into provided dict
    '''
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count != 0:
                cons_data_dict[row[0]] = row[1]
            line_count+=1
    return cons_data_dict

#%%
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

#%%
def get_annotation_sequence_index(seq_name,alignment):
    '''Get the sequence id of provided sequence name and alignment.
    '''
    for index, sequence in enumerate(alignment):
        if seq_name in sequence.id:
            annot_seq_ix = index
            break
    return annot_seq_ix


#%%
alignment = AlignIO.read(open('data/cons/LSU_AB.fas'), "fasta")


gap_to_nogap_ecoli = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUb_511145_ECOLI/1-2904',alignment)])
gap_to_nogap_pyrfu = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUa_186497_PYRFU/1-3048',alignment)])
#print(gap_to_nogap_pyrfu)
cons_data_dict = dict()
cons_data_dict = load_cons_data('data/cons/twc_AB_LSU.csv', cons_data_dict)

for aln_ix in cons_data_dict.keys():
    try:
        print('Alignment position:', aln_ix)
        print ('ECOLI:', gap_to_nogap_ecoli[int(aln_ix)], cons_data_dict[aln_ix])
        print ('PYRFU:', gap_to_nogap_pyrfu[int(aln_ix)], cons_data_dict[aln_ix])
    except:
        pass

