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
            if l[4] == 'Sequence':
                prev_set = 'Sequence'
                if 'Sequence' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Sequence'] = {}
                db[l[0]][l[1]][l[2]]['Sequence']['Family'] = l[6]
                db[l[0]][l[1]][l[2]]['Sequence']['Hmm'] = l[8]
                db[l[0]][l[1]][l[2]]['Sequence']['Score'] = l[10]
                db[l[0]][l[1]][l[2]]['Sequence']['Length'] = l[12]
            if l[4] == 'Labels':
                prev_set = 'Labels'
                if 'Labels' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Labels'] = {}
                db[l[0]][l[1]][l[2]]['Labels']['State'] = l[5]
                db[l[0]][l[1]][l[2]]['Labels']['DFG'] = l[6]
                db[l[0]][l[1]][l[2]]['Labels']['DFGAsp-rot'] = l[7]
                db[l[0]][l[1]][l[2]]['Labels']['XDF'] = l[8]
                db[l[0]][l[1]][l[2]]['Labels']['Chelix'] = l[9]
                db[l[0]][l[1]][l[2]]['Labels']['Saltbridge'] = l[10]
                db[l[0]][l[1]][l[2]]['Labels']['NT'] = l[11]
                db[l[0]][l[1]][l[2]]['Labels']['CT'] = l[12]
                db[l[0]][l[1]][l[2]]['Labels']['HRD'] = l[13]
                db[l[0]][l[1]][l[2]]['Labels']['Spine'] = l[14]
            if l[4] == 'APEdihe_Labels':
                db[l[0]][l[1]][l[2]]['Labels']['APE10-dihe'] = l[5]
                db[l[0]][l[1]][l[2]]['Labels']['APE9-dihe'] = l[6]
                db[l[0]][l[1]][l[2]]['Labels']['APE8-dihe'] = l[7]
                db[l[0]][l[1]][l[2]]['Labels']['APE8-rot'] = l[8]
                db[l[0]][l[1]][l[2]]['Labels']['APE67-dihe'] = l[9]
            if l[4] == 'APEdist_Labels':
                db[l[0]][l[1]][l[2]]['Labels']['APE12-dist'] = l[5]
                db[l[0]][l[1]][l[2]]['Labels']['APE11-dist'] = l[6]
                db[l[0]][l[1]][l[2]]['Labels']['APE10-dist'] = l[7]
                db[l[0]][l[1]][l[2]]['Labels']['APE9-dist'] = l[8]
            if l[4] == 'NTdom_Residues':
                prev_set = 'Residues'
                key = ''
                if 'Residues' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Residues'] = {}
                for e,entry in enumerate(('.'.join(l[5:])).split('.')):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]]['Residues'][key] = entry
            if l[4] == 'CTdom_Residues':
                key = ''
                if 'Residues' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Residues'] = {}
                for e,entry in enumerate(('.'.join(l[5:])).split('.')):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]]['Residues'][key] = entry
            if l[4] == 'DFG_Residues':
                key = ''
                if 'Residues' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Residues'] = {}
                for e,entry in enumerate(('.'.join(l[5:])).split('.')):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]]['Residues'][key] = entry
            if l[4] == 'APE_Residues':
                key = ''
                if 'Residues' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Residues'] = {}
                for e,entry in enumerate(('.'.join(l[5:])).split('.')):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]]['Residues'][key] = entry
            if l[4] == 'HRD_bbDihedrals':
                prev_set = 'Dihedrals'
                key = ''
                if 'Dihedrals' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Dihedrals'] = {}
                for e,entry in enumerate(l[5:]):
                    try:
                        float(entry)
                        values.append(float(entry))
                    except ValueError:
                        if e > 0:
                            db[l[0]][l[1]][l[2]]['Dihedrals'][key] = values
                        key = entry
                        values = []
                db[l[0]][l[1]][l[2]]['Dihedrals'][key] = values
            if l[4] == 'XDFG_bbDihedrals':
                key = ''
                for e,entry in enumerate(l[5:]):
                    try:
                        float(entry)
                        values.append(float(entry))
                    except ValueError:
                        if e > 0:
                            db[l[0]][l[1]][l[2]]['Dihedrals'][key] = values
                        key = entry
                        values = []
                db[l[0]][l[1]][l[2]]['Dihedrals'][key] = values
            if l[4] == 'APE_bbDihedrals':
                key = ''
                for e,entry in enumerate(l[5:]):
                    try:
                        float(entry)
                        values.append(float(entry))
                    except ValueError:
                        if e > 0:
                            db[l[0]][l[1]][l[2]]['Dihedrals'][key] = values
                        key = entry
                        values = []
                db[l[0]][l[1]][l[2]]['Dihedrals'][key] = values
            if l[4] == 'scDihedrals':
                key = ''
                for e,entry in enumerate(l[5:]):
                    try:
                        float(entry)
                        values.append(float(entry))
                    except ValueError:
                        if e > 0:
                            db[l[0]][l[1]][l[2]]['Dihedrals'][key].extend(values)
                        key = entry
                        values = []
                db[l[0]][l[1]][l[2]]['Dihedrals'][key].extend(values)

            if l[4] == 'Distances':
                prev_set = 'Distances'
                key = ''
                if 'Distances' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Distances'] = {}
                for e,entry in enumerate(l[5:]):
                    if e % 2 == 0: #Even entries are keys
                        key = entry
                    else:
                        db[l[0]][l[1]][l[2]]['Distances'][key] = entry

            if l[4] == 'Ligands':
                prev_set = 'Ligands'
                key = ''
                if 'Ligands' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Ligands'] = {}
                db[l[0]][l[1]][l[2]]['Ligands']['ID'] = {entry.split(':')[0]: entry.split(':')[1] for entry in l[5].split(',') if ':' in entry}
                if not db[l[0]][l[1]][l[2]]['Ligands']['ID']:
                    db[l[0]][l[1]][l[2]]['Ligands']['ID']['No_ligand'] = 'None'

                db[l[0]][l[1]][l[2]]['Ligands']['Type'] = {l[5].split(',')[e].split(':')[0]:entry for e,entry in enumerate(l[6].split(','))}
            if l[4] == 'Actloop':
                prev_set = 'Actloop'
                key = ''
                if 'Actloop' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Actloop'] = {}
                db[l[0]][l[1]][l[2]]['Actloop']['Missing'] = l[6]
                db[l[0]][l[1]][l[2]]['Actloop']['Length'] = l[8]
                db[l[0]][l[1]][l[2]]['Actloop']['Bfac'] = {}
                db[l[0]][l[1]][l[2]]['Actloop']['Bfac']['avg'] = l[10]
                db[l[0]][l[1]][l[2]]['Actloop']['Bfac']['min'] = l[12]
                db[l[0]][l[1]][l[2]]['Actloop']['Bfac']['max'] = l[14]
            if l[4] == 'Autoinhibit':
                prev_set = 'Autoinhibit'
                key = ''
                if 'Autoinhibit' not in db[l[0]][l[1]][l[2]].keys():
                    db[l[0]][l[1]][l[2]]['Autoinhibit'] = {}
                db[l[0]][l[1]][l[2]]['Autoinhibit']['HRD3'] = l[6]
                db[l[0]][l[1]][l[2]]['Autoinhibit']['Pseudosub'] = l[5]
                db[l[0]][l[1]][l[2]]['Autoinhibit']['Distance'] = l[8]


with open(str(sys.argv[2]), 'w') as jsonfile:
    json.dump(db, jsonfile, indent=4)
