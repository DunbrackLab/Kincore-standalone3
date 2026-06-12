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
#   print(kincore['1QMZA']['Labels']['State']) #Print out the activity state of CMGC_CDK2_HUMAN_1QMZA
#> 'Active'

db = {}
with open(sys.argv[1], 'r') as f:
    for line in tqdm.tqdm(f.readlines()):
        if line == '\n':
            if prev_set != 'Autoinhibit': #Check if the Autoinhibit entry was present in the previous PDB
                db[l[0]][l[1]][l[2]]['Autoinhibit'] = {}
                db[l[0]][l[1]][l[2]]['Autoinhibit']['HRD3'] = 'None'
                db[l[0]][l[1]][l[2]]['Autoinhibit']['Pseudosub'] = 'None'
                db[l[0]][l[1]][l[2]]['Autoinhibit']['Distance'] = 'None'
        else:
            l = line.strip().split()
            if l[0] not in db.keys():
                db[l[0]] = {}
            if l[1] not in db[l[0]].keys():
                db[l[0]][l[1]] = {}
            if l[2] not in db[l[0]][l[1]].keys():
                db[l[0]][l[1]][l[2]] = {}
            if l[4] == 'Labels':
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                db[l[0]][l[1]][l[2]][l[4]]['Family'] = l[6]
                db[l[0]][l[1]][l[2]][l[4]]['Hmm'] = l[8]
                db[l[0]][l[1]][l[2]][l[4]]['Score'] = l[10]
                db[l[0]][l[1]][l[2]][l[4]]['State'] = l[11]
                db[l[0]][l[1]][l[2]][l[4]]['DFG'] = l[12]
                db[l[0]][l[1]][l[2]][l[4]]['DFGAsp-rot'] = l[13]
                db[l[0]][l[1]][l[2]][l[4]]['XDF'] = l[14]
                db[l[0]][l[1]][l[2]][l[4]]['Chelix'] = l[15]
                db[l[0]][l[1]][l[2]][l[4]]['Saltbridge'] = l[16]
                db[l[0]][l[1]][l[2]][l[4]]['NT'] = l[17]
                db[l[0]][l[1]][l[2]][l[4]]['CT'] = l[18]
                db[l[0]][l[1]][l[2]][l[4]]['Spine'] = l[19]
                db[l[0]][l[1]][l[2]][l[4]]['APE10-dihe'] = l[20]
                db[l[0]][l[1]][l[2]][l[4]]['APE9-dihe'] = l[21]
                db[l[0]][l[1]][l[2]][l[4]]['APE8-dihe'] = l[22]
                db[l[0]][l[1]][l[2]][l[4]]['APE8-rot'] = l[23]
                db[l[0]][l[1]][l[2]][l[4]]['APE67-dihe'] = l[24]
                db[l[0]][l[1]][l[2]][l[4]]['APE12-dist'] = l[25]
                db[l[0]][l[1]][l[2]][l[4]]['APE11-dist'] = l[26]
                db[l[0]][l[1]][l[2]][l[4]]['APE10-dist'] = l[27]
                db[l[0]][l[1]][l[2]][l[4]]['APE9-dist'] = l[28]
                db[l[0]][l[1]][l[2]][l[4]]['HRD'] = l[29]

            if l[4] == 'Residues':
                prev_set = l[4]
                key = ''
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                for e,entry in enumerate(('.'.join(l[5:])).split('.')):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]][l[4]][key] = entry

            if l[4] == 'Distances':
                prev_set = l[4]
                key = ''
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                for e,entry in enumerate(l[5:]):
                    if key == 'Spine':
                        db[l[0]][l[1]][l[2]][l[4]][key] = l[(5+e):]
                        key = ''
                        break
                    elif e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]][l[4]][key] = entry
            if l[4] == 'Dihedrals':
                prev_set = l[4]
                key = ''
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                for e,entry in enumerate(l[5:]):
                    try:
                        float(entry)
                        values.append(float(entry))
                    except ValueError:
                        if e > 0:
                            db[l[0]][l[1]][l[2]][l[4]][key] = values
                        key = entry
                        values = []
                db[l[0]][l[1]][l[2]][l[4]][key] = values

            if l[4] == 'Ligands':
                prev_set = l[4]
                key = ''
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                db[l[0]][l[1]][l[2]][l[4]]['ID'] = {entry.split(':')[0]: entry.split(':')[1] for entry in l[5].split(',') if ':' in entry}
                if not db[l[0]][l[1]][l[2]][l[4]]['ID']:
                    db[l[0]][l[1]][l[2]][l[4]]['ID']['No_ligand'] = 'None'

                db[l[0]][l[1]][l[2]][l[4]]['Type'] = {l[5].split(',')[e].split(':')[0]:entry for e,entry in enumerate(l[6].split(','))}
            if l[4] == 'Actloop':
                prev_set = l[4]
                key = ''
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                db[l[0]][l[1]][l[2]][l[4]]['Length'] = l[5]
                db[l[0]][l[1]][l[2]][l[4]]['Bfac'] = {}
                db[l[0]][l[1]][l[2]][l[4]]['Bfac']['avg'] = l[6]
                db[l[0]][l[1]][l[2]][l[4]]['Bfac']['min'] = l[7]
                db[l[0]][l[1]][l[2]][l[4]]['Bfac']['max'] = l[8]
            if l[4] == 'Autoinhibit':
                prev_set = l[4]
                key = ''
                if l[4] not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]][l[4]] = {}
                db[l[0]][l[1]][l[2]][l[4]]['HRD3'] = l[6]
                db[l[0]][l[1]][l[2]][l[4]]['Pseudosub'] = l[5]
                db[l[0]][l[1]][l[2]][l[4]]['Distance'] = l[8]


with open(str(sys.argv[2]), 'w') as jsonfile:
    json.dump(db, jsonfile, indent=4)
