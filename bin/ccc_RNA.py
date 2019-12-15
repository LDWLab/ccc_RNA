#!/usr/bin/env python3

#%%
#Import required modules

import pandas as pd
import csv
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')

import bin.determine_universality as determine_universality

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
            if line_count > 1:
                base_pairs.append((int(row[1]),int(row[2])))
            line_count+=1
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

#%%
#Load and parse data

thet8_chaindict = load_chain_data('data/pdb/4qcn_chain_dict.csv')
pyrfu_chaindict = load_chain_data('data/pdb/3j2l-3j21_chain_dict.csv')

alignment = AlignIO.read(open('data/cons/LSU_AB.fas'), "fasta")
gap_to_nogap_thet8 = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUb_262724_THET2/1-2912',alignment)])
gap_to_nogap_pyrfu = gap_to_nogap_construct(alignment[get_annotation_sequence_index('LSUa_186497_PYRFU/1-3048',alignment)])

alnpos_to_twc = load_cons_data('data/cons/twc_AB_LSU.csv')

tt_base_pairs = base_pairs('data/contacts/TT_LSU_3D_BasePairs.csv') #Completely wrong...
pf_base_pairs = base_pairs('data/contacts/PyFu_LSU_BasePairs.csv')  #Manually fixed by adding 1 to all

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

ttlow_bps, ttrand_bps, tthigh_bps = filter_on_twc(tt_base_pairs, ttbp_to_twc, -1, 5)
pflow_bps, pfrand_bps, pfhigh_bps = filter_on_twc(pf_base_pairs, pfbp_to_twc, -1, 5)

pflowspec_bps, pfrandspec_bps, pfhighspec_bps = filter_on_twc(pfspec_bps, pfbp_to_twc, -1, 5)
#%%
print("P(A):", len(pfspec_bps)/len(pf_base_pairs))
print("P(B):", len(pflow_bps)/len(pf_base_pairs))
print("P(B|A)", len(pflowspec_bps)/len(pfspec_bps))

#%%
print(len(pfspec_bps), len(pflow_bps), len(pflowspec_bps), len(pfrandspec_bps), len(pfhighspec_bps))
print(len(pfspec_bps), len(pflow_bps), len(pflow_bps), len(pfrand_bps), len(pfhigh_bps))
print(len(pflowspec_bps)/len(pflow_bps), len(pfrandspec_bps)/len(pfrand_bps), len(pfhighspec_bps)/len(pfhigh_bps))


#%%
for bps in pflow_bps:
    print(bps, pfbp_to_twc[bps], bool(bps in pflowspec_bps))


# #%%
# #Initiate lists for dataframe and fill them with data
# nucl_pos = []
# universality = []
# nucl_twc = []
# distance_contact = []
# ttnucl_to_lowtwc = dict()
# ttnucl_to_hightwc = dict()

# for aln_ix in sorted(alnpos_to_twc.keys()):
#     if gap_to_nogap_thet8.get(int(aln_ix)) is None:
#         continue
#     if thet8_contacts.get(gap_to_nogap_thet8[int(aln_ix)]) is None:
#         continue
#     if float(alnpos_to_twc[aln_ix]) < -0.5:
#         ttnucl_to_lowtwc[gap_to_nogap_thet8[int(aln_ix)]] = float(alnpos_to_twc[aln_ix])
#     if float(alnpos_to_twc[aln_ix]) > 5:
#         ttnucl_to_hightwc[gap_to_nogap_thet8[int(aln_ix)]] = float(alnpos_to_twc[aln_ix])
#     for contact in thet8_contacts[gap_to_nogap_thet8[int(aln_ix)]]:
#         tt_polymer = thet8_chaindict[str(contact[0])]
#         nucl_ix = gap_to_nogap_thet8[int(aln_ix)]
#         contact_universality = determine_universality.main('THET8', tt_polymer, int(contact[1]))
#         #print(nucl_ix, contact_universality, alnpos_to_twc[aln_ix], contact[2])
#         nucl_pos.append(nucl_ix)
#         universality.append(contact_universality)
#         nucl_twc.append(float(alnpos_to_twc[aln_ix]))
#         distance_contact.append(contact[2])


# #%%
# #Initiate lists for dataframe and fill them with data for pyrfu
# nucl_pos = []
# universality = []
# nucl_twc = []
# distance_contact = []
# pfnucl_to_lowtwc = dict()
# pfnucl_to_hightwc = dict()
# pfnucl_to_random = dict()
# for aln_ix in sorted(alnpos_to_twc.keys()):
#     if gap_to_nogap_pyrfu.get(int(aln_ix)) is None:
#         continue
#     if pyrfu_contacts.get(gap_to_nogap_pyrfu[int(aln_ix)]) is None:
#         continue
#     if float(alnpos_to_twc[aln_ix]) > 5:
#         pfnucl_to_hightwc[gap_to_nogap_pyrfu[int(aln_ix)]] = float(alnpos_to_twc[aln_ix])
#     if float(alnpos_to_twc[aln_ix]) < -0.5:
#         pfnucl_to_lowtwc[gap_to_nogap_pyrfu[int(aln_ix)]] = float(alnpos_to_twc[aln_ix])
#     if (float(alnpos_to_twc[aln_ix]) > -0.5) and (float(alnpos_to_twc[aln_ix]) < 0.5):
#         pfnucl_to_random[gap_to_nogap_pyrfu[int(aln_ix)]] = float(alnpos_to_twc[aln_ix])
#     for contact in pyrfu_contacts[gap_to_nogap_pyrfu[int(aln_ix)]]:
#         pf_polymer = pyrfu_chaindict[str(contact[0])]
#         nucl_ix = gap_to_nogap_pyrfu[int(aln_ix)]
#         contact_universality = determine_universality.main('PYRFU', pf_polymer, int(contact[1]))
#         #print(nucl_ix, contact_universality, alnpos_to_twc[aln_ix], contact[2])
#         nucl_pos.append(nucl_ix)
#         universality.append(contact_universality)
#         nucl_twc.append(float(alnpos_to_twc[aln_ix]))
#         distance_contact.append(contact[2])


# #%%
# #Build dataframe
# tt_df = pd.DataFrame({
#     'Nucleotide index': nucl_pos,
#     'Universality of contact': universality,
#     'TwinCons conservation': nucl_twc,
#     'Contact distance': distance_contact
# })

# #%%
# #Print dataframe
# #print(tt_df)
# #tt_df[tt_df['Universality of contact'].str.contains('uni', na=False)]



# #%%PF

# pfcontact_dict_neg = dict()
# for base_pair in pf_base_pairs:
#     if (base_pair[0] not in pfnucl_to_lowtwc.keys()):
#         continue
#     if (base_pair[1] not in pfnucl_to_lowtwc.keys()):
#         continue
#     if base_pair not in pfcontact_dict_neg.keys():
#         pfcontact_dict_neg[base_pair] = []
#     for contact in pyrfu_contacts[base_pair[1]]:
#         tt_polymer = pyrfu_chaindict[str(contact[0])]
#         contact_universality = determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
#         pfcontact_dict_neg[base_pair].append((contact_universality,pfnucl_to_lowtwc[base_pair[1]], contact[2]))
#         #print (base_pair[1], tt_polymer,contact_universality, pfnucl_to_lowtwc[base_pair[1]], contact[2])
#     for contact in pyrfu_contacts[base_pair[0]]:
#         tt_polymer = pyrfu_chaindict[str(contact[0])]
#         contact_universality = determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
#         pfcontact_dict_neg[base_pair].append((contact_universality,pfnucl_to_lowtwc[base_pair[0]], contact[2]))
#         #print (base_pair[0], tt_polymer,contact_universality, pfnucl_to_lowtwc[base_pair[0]], contact[2])

# #%%
# pfcontact_dict_rand = dict()
# for base_pair in pf_base_pairs:
#     if (base_pair[0] not in pfnucl_to_random.keys()):
#         continue
#     if (base_pair[1] not in pfnucl_to_random.keys()):
#         continue
#     if base_pair not in pfcontact_dict_rand.keys():
#         pfcontact_dict_rand[base_pair] = []
#     for contact in pyrfu_contacts[base_pair[1]]:
#         tt_polymer = pyrfu_chaindict[str(contact[0])]
#         contact_universality = determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
#         pfcontact_dict_rand[base_pair].append((contact_universality,pfnucl_to_random[base_pair[1]], contact[2]))
#     for contact in pyrfu_contacts[base_pair[0]]:
#         tt_polymer = pyrfu_chaindict[str(contact[0])]
#         contact_universality = determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
#         pfcontact_dict_rand[base_pair].append((contact_universality,pfnucl_to_random[base_pair[0]], contact[2]))

# #%%
# pfcontact_dict_pos = dict()
# for base_pair in pf_base_pairs:
#     if (base_pair[0] not in pfnucl_to_hightwc.keys()):
#         continue
#     if (base_pair[1] not in pfnucl_to_hightwc.keys()):
#         continue
#     if base_pair not in pfcontact_dict_pos.keys():
#         pfcontact_dict_pos[base_pair] = []
#     for contact in pyrfu_contacts[base_pair[1]]:
#         tt_polymer = pyrfu_chaindict[str(contact[0])]
#         contact_universality = determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
#         pfcontact_dict_pos[base_pair].append((contact_universality,pfnucl_to_hightwc[base_pair[1]], contact[2]))
#     for contact in pyrfu_contacts[base_pair[0]]:
#         tt_polymer = pyrfu_chaindict[str(contact[0])]
#         contact_universality = determine_universality.main('PYRFU', tt_polymer, int(contact[1]))
#         pfcontact_dict_pos[base_pair].append((contact_universality,pfnucl_to_hightwc[base_pair[0]], contact[2]))

# #%%
# simplified_contact_type_neg = dict()
# simplified_contact_type_pos = dict()
# neg_list=[]
# pos_list=[]
# pos_list=[]
# rand_list=[]
# for base_pair, contacts_list in pfcontact_dict_neg.items():
#     uni = 0
#     spec = 0
#     for cont in contacts_list:
#         if cont[0] == 'uni':
#             uni+=1
#         elif cont[0] == 'spec':
#             spec+=1
#         else:
#             pass
#     simplified_contact_type_neg[base_pair] = (uni,spec)
#     neg_list.append((uni,spec))

# for base_pair, contacts_list in pfcontact_dict_pos.items():
#     uni = 0
#     spec = 0
#     for cont in contacts_list:
#         if cont[0] == 'uni':
#             uni+=1
#         elif cont[0] == 'spec':
#             spec+=1
#         else:
#             pass
#     simplified_contact_type_pos[base_pair] = (uni,spec)
#     pos_list.append((uni,spec))

# for base_pair, contacts_list in pfcontact_dict_rand.items():
#     uni = 0
#     spec = 0
#     for cont in contacts_list:
#         if cont[0] == 'uni':
#             uni+=1
#         elif cont[0] == 'spec':
#             spec+=1
#         else:
#             pass
#     simplified_contact_type_pos[base_pair] = (uni,spec)
#     rand_list.append((uni,spec))

# #%%
# fig = plt.figure(figsize=(10,10))
# ax1 = fig.add_subplot(111)
# #ax2 = fig.add_subplot(132)
# #ax3 = fig.add_subplot(133)
# ax1.scatter(*zip(*pos_list), c='r', label = 'Base-pairs above twc 5')
# #ax1.scatter(*zip(*neg_list), c='b', label = 'Base-pairs bellow twc -0.5')
# #ax1.scatter(*zip(*rand_list), c='g', label = 'Base-pairs random twc')
# ax1.set_xlabel('Universal contacts')
# ax1.set_ylabel('Dom specific contacts')
# ax1.set(xlim=(-0.1,10), ylim=(-0.1,10))
# plt.legend()


# #%%
# x=zip(*pos_list)
# sns.kdeplot(x[0],x[1], cmap="Reds", shade=True)
