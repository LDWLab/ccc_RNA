#!/usr/bin/env python3

#%%
#Import required modules

import pandas as pd
import csv
import os
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

import determine_universality

#Hacky way of fixing jupyter path weirdness:
def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')

if not is_interactive():
    os.chdir('../')
#%%
#Initialize functions for loading and parsing data.

def load_chain_data(csv_location):
    '''Load chain -> polymer name dictionaries from Nick.
    '''
    out_dict = dict()
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            out_dict[row[1]] = row[0]
    return out_dict

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

def base_pairs(csv_location):
    '''Create a list of base-pairs in the form of nucl,nucl tuples from a csv.
    '''
    base_pairs=[]
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count < 2:
                line_count+=1
                continue
            base_pairs.append((int(row[1]),int(row[2])))
    return base_pairs

def load_contacts(contacts_location, rrna_chain, dist_thr):
    '''Given a RING2.0 contacts csv, a chain of interest and threshold distance,
    returns dictionary of interactions with key residue number and value list of contacts.
    '''
    contacts_output = dict()
    with open(contacts_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        first_line = True
        for row in csv_reader:
            if first_line:
                first_line = False
                continue
            if float(row[3]) > float(dist_thr):
                continue
            if (str(row[0].split(':')[0]) == rrna_chain) and (len(row[2].split(':')[3]) == 3) and (len(row[0].split(':')[3]) == 1):
                if int(row[0].split(":")[1]) not in contacts_output.keys():
                    contacts_output[int(row[0].split(":")[1])] = []
                contacts_output[int(row[0].split(":")[1])].append((str(row[2].split(":")[0]),int(row[2].split(":")[1]),float(row[3])))
                #contacts_list.append([row[0],row[2],row[3],row[6],row[7]])
            if (str(row[2].split(':')[0]) == rrna_chain) and (len(row[0].split(':')[3]) == 3) and (len(row[2].split(':')[3]) == 1):
                if int(row[2].split(":")[1]) not in contacts_output.keys():
                    contacts_output[int(row[2].split(":")[1])] = []
                contacts_output[int(row[2].split(":")[1])].append((str(row[0].split(":")[0]),int(row[0].split(":")[1]),float(row[3])))
                #contacts_list.append([row[2],row[0],row[3],row[6],row[7]])
    return contacts_output

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

def construct_bp_to_contact(base_pairs, prot_contacts):
    '''Given a list of base-pairs and a interaction dictionary, produces
    base-pair to interactions mapping. Keys are tuples, values are lists.
    Looses the information which nucleotide of the bp is doing which interaction.
    '''
    bp_to_contact = dict()
    for base_pair in base_pairs:
        if (base_pair[0] not in prot_contacts.keys()) and (base_pair[1] not in prot_contacts.keys()):
            continue
        if base_pair not in bp_to_contact.keys():
            bp_to_contact[base_pair] = []
        if (base_pair[0] in prot_contacts.keys()):
            for contact in prot_contacts[base_pair[0]]:
                bp_to_contact[base_pair].append(contact)
        if (base_pair[1] in prot_contacts.keys()):
            for contact in prot_contacts[base_pair[1]]:
                bp_to_contact[base_pair].append(contact)
    return bp_to_contact

def filter_on_twc(base_pairs, bp_to_twc, low_thr, high_thr):
    '''Returns 3 lists of filtered bps based on twc.
    First list is bps with both nucleotides bellow low_thr;
    Second list is bps with both nucleotides between thresholds;
    Third list is bps with both nucleotides above high_thr
    '''
    lowtwc_bps = []
    randomtwc_bps = []
    hightwc_bps = []
    for base_pair in base_pairs:
        if base_pair not in bp_to_twc.keys():
            continue
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

def filter_on_specificity(base_pairs, bp_to_contact, species, chaindict):
    '''Returns a list of filtered bps that have at least 1 domain specific contact.
    Looses information about number and type of contacts with bps.
    '''
    bps_with_spec_contact = []
    for base_pair in base_pairs:
        if base_pair not in bp_to_contact.keys():
            continue
        for contact in bp_to_contact[base_pair]:
            if determine_universality.main(species, chaindict[str(contact[0])], int(contact[1])) == 'uni':
                continue
            if determine_universality.main(species, chaindict[str(contact[0])], int(contact[1])) == 'spec':
                bps_with_spec_contact.append(base_pair)
                break
    return bps_with_spec_contact

def unique_tup_list(tup_list):
    list_of_sorted_value_tups = []
    for tup in tup_list:
        list_of_sorted_value_tups.append((sorted(tup)[0],sorted(tup)[1]))
    return set(list_of_sorted_value_tups)


#%%
#Load and parse data

thet8_chaindict = load_chain_data('data/pdb/4qcn_chain_dict.csv')
pyrfu_chaindict = load_chain_data('data/pdb/3j2l-3j21_chain_dict.csv')

alignment = AlignIO.read(open('data/cons/LSU_AB.fas'), "fasta")
gap_to_nogap_thet8 = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUb_262724_THET2/1-2912',alignment)])
gap_to_nogap_pyrfu = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUa_186497_PYRFU/1-3048',alignment)])

alnpos_to_twc = load_cons_data('data/cons/twc_AB_LSU.csv')

tt_base_pairs = base_pairs('data/contacts/TT_LSU_3D_BasePairs.csv') #Completely wrong...
pf_base_pairs = base_pairs('data/contacts/PyFu_LSU_BasePairs.csv')  #Manually fixed by adding 1 to all in the datafile

thet8_contacts = load_contacts('data/contacts/4qcn_edges.txt', 'A', 4)
pyrfu_contacts = load_contacts('data/contacts/3j2l-3j21_merge_edges.txt', '1', 4)

#Carefully check the following functions
ttnucl_to_twc = construct_nucl_to_twc(alnpos_to_twc, gap_to_nogap_thet8)
pfnucl_to_twc = construct_nucl_to_twc(alnpos_to_twc, gap_to_nogap_pyrfu)

ttbp_to_twc = construct_bp_to_twc(ttnucl_to_twc, tt_base_pairs)
pfbp_to_twc = construct_bp_to_twc(pfnucl_to_twc, pf_base_pairs)

ttbp_to_contacts = construct_bp_to_contact(tt_base_pairs, thet8_contacts)
pfbp_to_contacts = construct_bp_to_contact(pf_base_pairs, pyrfu_contacts)

ttspec_bps = filter_on_specificity(tt_base_pairs, ttbp_to_contacts, 'THET8', thet8_chaindict)
pfspec_bps = filter_on_specificity(pf_base_pairs, pfbp_to_contacts, 'PYRFU', pyrfu_chaindict)
#%%

ttlow_bps, ttrand_bps, tthigh_bps = filter_on_twc(tt_base_pairs, ttbp_to_twc, -1.5, 5)
pflow_bps, pfrand_bps, pfhigh_bps = filter_on_twc(pf_base_pairs, pfbp_to_twc, -1.5, 5)

pflowspec_bps, pfrandspec_bps, pfhighspec_bps = filter_on_twc(pfspec_bps, pfbp_to_twc, -1.5, 5)
#%%
print("P(A):", len(pfspec_bps)/len(pf_base_pairs))
print("P(B):", len(pflow_bps)/len(pf_base_pairs))
print("P(B|A)", len(pflowspec_bps)/len(pfspec_bps))

#%%
print(len(pfspec_bps), len(pflow_bps), len(pflowspec_bps), len(pfrandspec_bps), len(pfhighspec_bps))
print(len(pfspec_bps), len(pflow_bps), len(pflow_bps), len(pfrand_bps), len(pfhigh_bps))
print(len(pflowspec_bps)/len(pflow_bps), len(pfrandspec_bps)/len(pfrand_bps), len(pfhighspec_bps)/len(pfhigh_bps))

#%%Uniques

print(len(unique_tup_list(pfspec_bps)), len(unique_tup_list(pflow_bps)), \
    len(unique_tup_list(pflowspec_bps)), len(unique_tup_list(pfrandspec_bps)), len(unique_tup_list(pfhighspec_bps)))

print(len(unique_tup_list(pflowspec_bps))/len(unique_tup_list(pflow_bps)), \
    len(unique_tup_list(pfrandspec_bps))/len(unique_tup_list(pfrand_bps)), \
    len(unique_tup_list(pfhighspec_bps))/len(unique_tup_list(pfhigh_bps)))
#%%
for bps in unique_tup_list(pflow_bps):
    print(bps, pfbp_to_twc[bps], bool(bps in pflowspec_bps))

#%%

pflowspec_nucl = []
pflownonspec_nucl = []
for nucl,twc in pfnucl_to_twc.items():
    if twc > -1.5:
        continue
    if (any(nucl in i for i in pflow_bps)):
        continue
    if nucl not in pyrfu_contacts.keys():
        pflownonspec_nucl.append(nucl)
        continue
    if (any('spec' in determine_universality.main('PYRFU', pyrfu_chaindict[str(i[0])], int(i[1])) for i in pyrfu_contacts[nucl])):
        pflowspec_nucl.append(nucl)
        continue
    pflownonspec_nucl.append(nucl)

print(pflowspec_nucl)
print(pflownonspec_nucl)


# %%

lowperc = []
randperc = []
highperc = []
distance_list = []
lspec_num = []
rspec_num = []
hspec_num = []
for i in range(1,7):
    distance_list.append(i)
    pyrfu_contacts = load_contacts('data/contacts/3j2l-3j21_merge_edges.txt', '1', i)
    pfnucl_to_twc = construct_nucl_to_twc(alnpos_to_twc, gap_to_nogap_pyrfu)
    pfbp_to_twc = construct_bp_to_twc(pfnucl_to_twc, pf_base_pairs)
    pfbp_to_contacts = construct_bp_to_contact(pf_base_pairs, pyrfu_contacts)
    pfspec_bps = filter_on_specificity(pf_base_pairs, pfbp_to_contacts, 'PYRFU', pyrfu_chaindict)
    pflow_bps, pfrand_bps, pfhigh_bps = filter_on_twc(pf_base_pairs, pfbp_to_twc, -1.5, 5)
    pflowspec_bps, pfrandspec_bps, pfhighspec_bps = filter_on_twc(pfspec_bps, pfbp_to_twc, -1.5, 5)
    lowperc.append(len(pflowspec_bps)/len(pflow_bps))
    randperc.append(len(pfrandspec_bps)/len(pfrand_bps))
    highperc.append(len(pfhighspec_bps)/len(pfhigh_bps))
    lspec_num.append(len(pflowspec_bps))
    rspec_num.append(len(pfrandspec_bps))
    hspec_num.append(len(pfhighspec_bps))

#%%
spec_num = hspec_num+lspec_num+rspec_num

df = pd.DataFrame(
    {'Distance': distance_list,
    'Low': lowperc,
    'Rand': randperc,
    'High': highperc})
df
df = df.melt('Distance', var_name='TwinCons',  value_name='Percentage specifics')
df['Absolute numbers'] = spec_num
df
#%%
%matplotlib inline
sns.scatterplot(x='Distance', y = 'Percentage specifics', data=df, hue='TwinCons',size='Absolute numbers')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)