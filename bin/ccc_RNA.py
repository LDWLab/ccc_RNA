#!/usr/bin/env python3

#%%
import pandas as pd
import csv
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

import bin.determine_universality

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
def base_pairs(csv_location):
    base_pairs=[]
    with open(csv_location) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count > 1:
                base_pairs.append((int(row[1])+1,int(row[2])+1))
            line_count+=1
    return base_pairs

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

#%%
#Initiate lists for dataframe and fill them with data
nucl_pos = []
universality = []
nucl_twc = []
distance_contact = []
nucl_to_lowtwc = dict()

for aln_ix in sorted(cons_data_dict.keys()):
    if gap_to_nogap_thet8.get(int(aln_ix)) is None:
        continue
    if thet8_contacts.get(gap_to_nogap_thet8[int(aln_ix)]) is None:
        continue
    if float(cons_data_dict[aln_ix]) < -0.5:
        nucl_to_lowtwc[gap_to_nogap_thet8[int(aln_ix)]] = float(cons_data_dict[aln_ix])
    for contact in thet8_contacts[gap_to_nogap_thet8[int(aln_ix)]]:
        tt_polymer = thet8_chaindict[str(contact[0])]
        nucl_ix = gap_to_nogap_thet8[int(aln_ix)]
        contact_universality = bin.determine_universality.main('THET8', tt_polymer, int(contact[1]))
        #print(nucl_ix, contact_universality, cons_data_dict[aln_ix], contact[2])
        nucl_pos.append(nucl_ix)
        universality.append(contact_universality)
        nucl_twc.append(float(cons_data_dict[aln_ix]))
        distance_contact.append(contact[2])


#%%
#Initiate lists for dataframe and fill them with data for pyrfu
nucl_pos = []
universality = []
nucl_twc = []
distance_contact = []
nucl_to_lowtwc = dict()
pyrfu_chaindict['6'] = 'eS24'
pyrfu_chaindict['4'] = 'aL08'
for aln_ix in sorted(cons_data_dict.keys()):
    if gap_to_nogap_pyrfu.get(int(aln_ix)) is None:
        continue
    if pyrfu_contacts.get(gap_to_nogap_pyrfu[int(aln_ix)]) is None:
        continue
    if float(cons_data_dict[aln_ix]) < -1:
        nucl_to_lowtwc[gap_to_nogap_pyrfu[int(aln_ix)]] = float(cons_data_dict[aln_ix])
    for contact in pyrfu_contacts[gap_to_nogap_pyrfu[int(aln_ix)]]:
        pf_polymer = pyrfu_chaindict[str(contact[0])]
        nucl_ix = gap_to_nogap_pyrfu[int(aln_ix)]
        contact_universality = bin.determine_universality.main('PYRFU', pf_polymer, int(contact[1]))
        #print(nucl_ix, contact_universality, cons_data_dict[aln_ix], contact[2])
        nucl_pos.append(nucl_ix)
        universality.append(contact_universality)
        nucl_twc.append(float(cons_data_dict[aln_ix]))
        distance_contact.append(contact[2])


#%%
#Build dataframe
tt_df = pd.DataFrame({
    'Nucleotide index': nucl_pos,
    'Universality of contact': universality,
    'TwinCons conservation': nucl_twc,
    'Contact distance': distance_contact
})

#%%
#Print dataframe
#print(tt_df)
#tt_df[tt_df['Universality of contact'].str.contains('uni', na=False)]

tt_base_pairs = base_pairs('data/contacts/TT_LSU_3D_BasePairs.csv')
pf_base_pairs = base_pairs('data/contacts/PyFu_LSU_BasePairs.csv')
print(pf_base_pairs)

#%%TT

# for nucl_ix, nucl_twc in nucl_to_lowtwc.items():
#     for contact in thet8_contacts[nucl_ix]:
#         tt_polymer = thet8_chaindict[str(contact[0])]
#         contact_universality = bin.determine_universality.main('THET8', tt_polymer, int(contact[1]))
#         print(nucl_ix, contact_universality, nucl_twc)

for base_pair in tt_base_pairs:
    if (base_pair[0] not in nucl_to_lowtwc.keys()):
        continue
    #if (base_pair[1] not in nucl_to_lowtwc.keys()):
    #    continue
    print(base_pair)
    for contact in thet8_contacts[base_pair[1]]:
        tt_polymer = thet8_chaindict[str(contact[0])]
        contact_universality = bin.determine_universality.main('THET8', tt_polymer, int(contact[1]))
        print (base_pair[1], contact_universality, nucl_to_lowtwc[base_pair[1]])
    for contact in thet8_contacts[base_pair[0]]:
        tt_polymer = thet8_chaindict[str(contact[0])]
        contact_universality = bin.determine_universality.main('THET8', tt_polymer, int(contact[1]))
        print (base_pair[0], contact_universality, nucl_to_lowtwc[base_pair[0]])

#%%PF
for base_pair in pf_base_pairs:
    if (base_pair[0] not in nucl_to_lowtwc.keys()):
        continue
    if (base_pair[1] not in nucl_to_lowtwc.keys()):
        continue
    print(base_pair)
    for contact in pyrfu_contacts[base_pair[1]]:
        tt_polymer = pyrfu_chaindict[str(contact[0])]
        contact_universality = bin.determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
        print (base_pair[1], tt_polymer,contact_universality, nucl_to_lowtwc[base_pair[1]], contact[2])
    for contact in pyrfu_contacts[base_pair[0]]:
        tt_polymer = pyrfu_chaindict[str(contact[0])]
        contact_universality = bin.determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
        print (base_pair[0], tt_polymer,contact_universality, nucl_to_lowtwc[base_pair[0]], contact[2])


#%%
#Plot things
plt.figure(figsize=(16,12))
g = sns.distplot(tt_df[tt_df['Universality of contact'].str.contains('uni', na=False)]['TwinCons conservation'], bins=20, )

#%%
unidata = list(tt_df[tt_df['Universality of contact'].str.contains('uni', na=False)]['TwinCons conservation'])
specdata = list(tt_df[tt_df['Universality of contact'].str.contains('spec', na=False)]['TwinCons conservation'])

fig, ax = plt.subplots()
for a in [unidata, specdata]:
    sns.distplot(a, ax=ax, kde=False, bins=10)

ax.legend(['Universal '+str(len(unidata)),'Specific '+str(len(specdata))])