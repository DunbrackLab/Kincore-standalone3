import sys
import tqdm
import json

"""
Created on July 17 2024 by JG
Last edited on Fri Feb 13, 2026 by JG
@author: joangizzio
"""

# Description: takes a text file containing the output of kinase_state.py and creates a json object

# Usage:
#   python make_json.py input.txt output.json

#The json file can be loaded and used as follows:
#   import json
#   with open('output.json', 'rb') as jsonfile:
#       kincore = json.load(jsonfile)
#   print(kincore['CMGC_CDK2_HUMAN_1QMZA']['Labels']['State']) #Print out the activity state of CMGC_CDK2_HUMAN_1QMZA
#> 'Active'

db = {}
with open(sys.argv[1], 'r') as f:
    for line in tqdm.tqdm(f.readlines()):
        if line == '\n':
            if prev_set != 'Autoinhibit': #Check if the Autoinhibit entry was present in the previous PDB
                db[l[1]]['Autoinhibit'] = {}
                db[l[1]]['Autoinhibit']['HRD3'] = 'None'
                db[l[1]]['Autoinhibit']['Pseudosub'] = 'None'
                db[l[1]]['Autoinhibit']['Distance'] = 'None'
        else:
            l = line.strip().split()
            if l[1] not in db.keys():
                db[l[1]] = {}
            if l[5] == 'Labels':
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                db[l[1]][l[5]]['Family'] = l[7]
                db[l[1]][l[5]]['Hmm'] = l[9]
                db[l[1]][l[5]]['Score'] = l[11]
                db[l[1]][l[5]]['State'] = l[12]
                db[l[1]][l[5]]['DFG'] = l[13]
                db[l[1]][l[5]]['XDF'] = l[14]
                db[l[1]][l[5]]['Chelix'] = l[15]
                db[l[1]][l[5]]['Saltbridge'] = l[16]
                db[l[1]][l[5]]['NT'] = l[17]
                db[l[1]][l[5]]['CT'] = l[18]
                db[l[1]][l[5]]['Spine'] = l[19]
                db[l[1]][l[5]]['APE10-dihe'] = l[20]
                db[l[1]][l[5]]['APE9-dihe'] = l[21]
                db[l[1]][l[5]]['APE8-dihe'] = l[22]
                db[l[1]][l[5]]['APE8-rot'] = l[23]
                db[l[1]][l[5]]['APE67-dihe'] = l[24]
                db[l[1]][l[5]]['APE12-dist'] = l[25]
                db[l[1]][l[5]]['APE11-dist'] = l[26]
                db[l[1]][l[5]]['APE10-dist'] = l[27]
                db[l[1]][l[5]]['APE9-dist'] = l[28]
                db[l[1]][l[5]]['HRD'] = l[29]

            if l[5] == 'Residues':
                prev_set = l[5]
                key = ''
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                for e,entry in enumerate(('.'.join(l[6:])).split('.')):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[1]][l[5]][key] = entry

            if l[5] == 'Distances':
                prev_set = l[5]
                key = ''
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                for e,entry in enumerate(l[6:]):
                    if key == 'Spine':
                        db[l[1]][l[5]][key] = l[(6+e):]
                        key = ''
                        break
                    elif e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[1]][l[5]][key] = entry
            if l[5] == 'Dihedrals':
                prev_set = l[5]
                key = ''
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                for e,entry in enumerate(l[6:]):
                    try:
                        float(entry)
                        values.append(float(entry))
                    except ValueError:
                        if e > 0:
                            db[l[1]][l[5]][key] = values
                        key = entry
                        values = []
                db[l[1]][l[5]][key] = values

            if l[5] == 'Ligands':
                prev_set = l[5]
                key = ''
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                db[l[1]][l[5]]['ID'] = {entry.split(':')[0]: entry.split(':')[1] for entry in l[6].split(',') if ':' in entry}
                if not db[l[1]][l[5]]['ID']:
                    db[l[1]][l[5]]['ID']['No_ligand'] = 'None'

                db[l[1]][l[5]]['Type'] = {l[6].split(',')[e].split(':')[0]:entry for e,entry in enumerate(l[7].split(','))}
            if l[5] == 'Actloop':
                prev_set = l[5]
                key = ''
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                db[l[1]][l[5]]['Length'] = l[6]
                db[l[1]][l[5]]['Bfac'] = {}
                db[l[1]][l[5]]['Bfac']['avg'] = l[7]
                db[l[1]][l[5]]['Bfac']['min'] = l[8]
                db[l[1]][l[5]]['Bfac']['max'] = l[9]
            if l[5] == 'Autoinhibit':
                prev_set = l[5]
                key = ''
                if l[5] not in db[l[1]].keys():
                    db[l[1]][l[5]] = {}
                db[l[1]][l[5]]['HRD3'] = l[7]
                db[l[1]][l[5]]['Pseudosub'] = l[6]
                db[l[1]][l[5]]['Distance'] = l[9]


with open(str(sys.argv[2]), 'w') as jsonfile:
    json.dump(db, jsonfile, indent=4)
