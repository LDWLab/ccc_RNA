#!/usr/bin/env python3

import os, json, glob

#In the future make into class to not load jsons every time it's executed.
#class JSONSpecies:
'''Class for loading JSON annotation of a species.'''
#Some data:
folder_dict = {
    'THET8':'Bacteria',
    'PYRFU':'Archaea',
    'HUMAN':'Eukarya'}


def load_jsons_from_folder(folder_path,spec_name):
    output_dict=dict()
    for json_file in glob.glob(folder_path+'/'+spec_name+'_*'):
        json_data=dict()
        with open(json_file) as F:
            json_data = json.load(F)
        output_dict[json_file.replace('.','_').split("_")[1]] = json_data
    return output_dict

def main(spec_name, polymer_name, resi_id):
    if polymer_name[0] != 'u':
        return 'spec'
    json_dicts = load_jsons_from_folder('data/json/'+folder_dict[spec_name],spec_name)
    for definition in json_dicts.keys():
        if (polymer_name not in json_dicts[definition].keys()):
            continue
        if not any(resi_id in range(r[0], r[1]) for r in json_dicts[definition][polymer_name]):
            continue
        if definition[0] in {'u','G'}:
            return 'uni'
        else:
            return 'spec'
        return None


#determine_universality('PYRFU','uL16',150)