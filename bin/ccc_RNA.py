#!/usr/bin/env python3

#%%
import pandas as pd
import csv
from Bio import AlignIO

import determine_universality

#%%
def load_chain_data(csv_location):
    '''Load chain -> polymer name dictionaries from Nick.
    '''
    out_dict = dict()
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            out_dict[row[1]] = row[0]
    return out_dict

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
def load_contacts(contacts_location, rrna_chain, dist_thr):
    contacts_output = dict()
    with open(contacts_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        first_line = True
        for row in csv_reader:
            if first_line:
                first_line = False
                continue
            if float(row[3]) > dist_thr:
                continue
            if (str(row[0].split(':')[0]) == rrna_chain) and (len(row[2].split(':')[3]) == 3) and (len(row[0].split(':')[3]) == 1):
                if row[0] not in contacts_output.keys():
                    contacts_output[row[0].split(":")[1]] = []
                contacts_output[row[0].split(":")[1]].append((str(row[2].split(":")[0]),int(row[2].split(":")[1])))
                #contacts_list.append([row[0],row[2],row[3],row[6],row[7]])
            if (str(row[2].split(':')[0]) == rrna_chain) and (len(row[0].split(':')[3]) == 3) and (len(row[2].split(':')[3]) == 1):
                if row[2] not in contacts_output.keys():
                    contacts_output[row[2].split(":")[1]] = []
                contacts_output[row[2].split(":")[1]].append((str(row[0].split(":")[0]),int(row[0].split(":")[1])))
                #contacts_list.append([row[2],row[0],row[3],row[6],row[7]])
    return contacts_output

#%%
alignment = AlignIO.read(open('data/cons/LSU_AB.fas'), "fasta")
thet8_chaindict = load_chain_data('data/pdb/4qcn_chain_dict.csv')
pyrfu_chaindict = load_chain_data('data/pdb/3j2l-3j21_chain_dict.csv')

gap_to_nogap_thet8 = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUb_262724_THET2/1-2912',alignment)])
gap_to_nogap_pyrfu = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUa_186497_PYRFU/1-3048',alignment)])
#print(gap_to_nogap_pyrfu)
cons_data_dict = dict()
cons_data_dict = load_cons_data('data/cons/twc_AB_LSU.csv', cons_data_dict)

thet8_contacts = load_contacts('data/contacts/4qcn_edges.txt', 'A', 3)
pyrfu_contacts = load_contacts('data/contacts/3j2l-3j21_merge_edges.txt', '1', 3)

for aln_ix in cons_data_dict.keys():
    try:
        print('Alignment position:', aln_ix)
        tt_resi_id = thet8_contacts[str(gap_to_nogap_thet8[int(aln_ix)])][0][1] #Limits to only one interactor with nucleotide# Fix here
        pf_resi_id = pyrfu_contacts[str(gap_to_nogap_pyrfu[int(aln_ix)])][0][1]
        tt_chain = thet8_contacts[str(gap_to_nogap_thet8[int(aln_ix)])][0][0]
        pf_chain = pyrfu_contacts[str(gap_to_nogap_pyrfu[int(aln_ix)])][0][0]
        print ('THET8:', determine_universality.main('THET8',tt_chain,tt_resi_id), cons_data_dict[aln_ix])
        print ('PYRFU:', determine_universality.main('PYRFU',pf_chain,pf_resi_id), cons_data_dict[aln_ix])
    except:
        pass
