#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:13:13 2020

@author: vivekmodi
added 2025 crieria for active kinases. RLD 12/12/2025
"""

import math
import numpy as np

# angle within range function
def dihe_in_range(left, right, value):
    if value > 900: return None
    if value >= left and value <= right:
        return True
    else:
        return False
    
def make_apelabel(apelabel):    
    apedihe="APEdihe_"
    if "APE10-dihe-none" in apelabel: apedihe = apedihe + 'n'
    if "APE10-dihe-in" in apelabel: apedihe = apedihe + 'i'
    if "APE10-dihe-na" in apelabel: apedihe = apedihe + 'a'
    if "APE10-dihe-out" in apelabel: apedihe = apedihe + 'o'
    if "APE9-dihe-none" in apelabel: apedihe = apedihe + 'n'
    if "APE9-dihe-in" in apelabel: apedihe = apedihe + 'i'
    if "APE9-dihe-na" in apelabel: apedihe = apedihe + 'a'
    if "APE9-dihe-out" in apelabel: apedihe = apedihe + 'o'
    if "APE8-dihe-none" in apelabel: apedihe = apedihe + 'n'
    if "APE8-dihe-in" in apelabel: apedihe = apedihe + 'i'
    if "APE8-dihe-na" in apelabel: apedihe = apedihe + 'a'
    if "APE8-dihe-out" in apelabel: apedihe = apedihe + 'o'
    if "APE8-rot-none" in apelabel: apedihe = apedihe + 'n'
    if "APE8-rot-in" in apelabel: apedihe = apedihe + 'i'
    if "APE8-rot-na" in apelabel: apedihe = apedihe + 'a'
    if "APE8-rot-out" in apelabel: apedihe = apedihe + 'o'
    if "APE67-dihe-none" in apelabel: apedihe = apedihe + 'n'
    if "APE67-dihe-in" in apelabel: apedihe = apedihe + 'i'
    if "APE67-dihe-na" in apelabel: apedihe = apedihe + 'a'
    if "APE67-dihe-out" in apelabel: apedihe = apedihe + 'o'

    apedist="APEdist_"
    if "APE12-dist-none" in apelabel: apedist =apedist+'n'
    if "APE12-dist-in" in apelabel: apedist =apedist+'i'
    if "APE12-dist-na" in apelabel: apedist =apedist+'a'
    if "APE12-dist-out" in apelabel: apedist =apedist+'o'
    if "APE11-dist-none" in apelabel: apedist =apedist+'n'
    if "APE11-dist-in" in apelabel: apedist =apedist+'i'
    if "APE11-dist-na" in apelabel: apedist =apedist+'a'
    if "APE11-dist-out" in apelabel: apedist =apedist+'o'
    if "APE10-dist-none" in apelabel: apedist =apedist+'n'
    if "APE10-dist-in" in apelabel: apedist =apedist+'i'
    if "APE10-dist-na" in apelabel: apedist =apedist+'a'
    if "APE10-dist-out" in apelabel: apedist =apedist+'o'
    if "APE9-dist-none" in apelabel: apedist =apedist+'n'
    if "APE9-dist-in" in apelabel: apedist =apedist+'i'
    if "APE9-dist-na" in apelabel: apedist =apedist+'a'
    if "APE9-dist-out" in apelabel: apedist =apedist+'o'

    return(apedihe, apedist)

def active_labels(pwd, df):
    print('Assigning active labels...')
    hbondcutoff=3.6

    OUTPUT=open(f'{pwd}/Active_labels.txt','w')
    OUTPUT.write("tyr kinase pdb activity dfg dihe SB CH_SB HRD NT CT NT_CT SB_NT_CT SNC APEdihe APEdist dihe_10 dihe_9 dihe_8 rot_8 dihe_67 dist_12 dist_11  dist_10  dist_9 his hisphi hispsi arg argphi argpsi ape10 ape10phi ape10psi ape9 ape9phi ape9psi ape8 ape8phi ape8psi ape8chi1 ape7 ape7phi ape7psi ape6 ape6phi ape6psi sb sbdist x xhrd_dist apedist ape12dist ape11dist ape10dist ape9dist\n")

    for i in df.index:

        # APEtype: TYR, nonTYR, noAPE
        group=df.at[i, 'hmm']
        if group == "TYR": APEtype="TYR"
        else: APEtype="nonTYR"
        # BUB, HASP, PKDCC, TP53RK exception, no APE motif; check
        # PEAK3 has group=PEAK and score~170 and no APE motif (PEAK1 and PRAG1, >500)
        #            if (group=="BUB" or group=="HASP" or group=="PKDCC" or
        if (group=="HASP" or group=="PKDCC" or
            group=="TP53RK" or group=="PAN3" or group=="RNASEL" or
            (group=="PEAK" and df.at[i, 'Score']<300)): APEtype="noAPE"
        df.at[i, 'APEtype']=APEtype
        
        # Determine final group
        finalgroup=group
        if group != "None":
            if group=="CDC7": finalgroup="OTHER"
            if group=="PAN3": finalgroup="OTHER"
            if group=="TP53RK": finalgroup="OTHER"
            if group=="PIK3R4": finalgroup="OTHER"
            if group=="RNASEL": finalgroup="OTHER"
            if group=="PKDCC": finalgroup="OTHER"
            if group=="BUB": finalgroup="OTHER"
            if group=="ULK": finalgroup="OTHER"
            if group=="WNK": finalgroup="OTHER"
            if group=="PINK1": finalgroup="OTHER"
            if group=="EIF2AK41": finalgroup="OTHER"
            if group=="PEAK": finalgroup="OTHER"
            if group=="POMK": finalgroup="OTHER"
            if group=="HASP": finalgroup="OTHER"
            if group=="SCYL": finalgroup="OTHER"
            if group=="STK31": finalgroup="OTHER"
            if group=="PXK": finalgroup="OTHER"
            if group=="MOS": finalgroup="OTHER"
            if group=="TBCK": finalgroup="OTHER"
            if group=="RPS6KC1": finalgroup="OTHER"
            if group=="MAP3K123": finalgroup="TKL"
            if df.at[i, 'Score']<180.0:
                finalgroup="OTHER"  # some OTHER kinases have higher scores for CAMK, CMGC, etc than the OTHER HMM

        # structural parameters for active state determination
        dis_saltbridge=float(df.at[i,'LysNZ_GluOE_dis'])
        dis_xhrd=float(df.at[i,'DFG6_XHRD_dis'])
        dis_ape9=float(df.at[i,'APE9_Arg_dis'])
        dis_ape10=df.at[i, 'APE10_DFG4_dis']
        dis_ape11=df.at[i, 'APE11_DFG4_dis']
        dis_ape12=df.at[i, 'APE12_DFG4_dis']
        HRD_Phi=df.at[i, 'HRD_Phi'] ;     HRD_Psi=df.at[i, 'HRD_Psi']
        Arg_Phi=df.at[i, 'Arg_Phi'] ;     Arg_Psi=df.at[i, 'Arg_Psi']
        XDFG_Phi=df.at[i, 'XDFG_Phi'] ;       XDFG_Psi=df.at[i, 'XDFG_Psi']
        DFG_Phi=df.at[i, 'DFG_Phi'] ;     DFG_Psi=df.at[i, 'DFG_Psi']
        Phe_Phi=df.at[i, 'Phe_Phi'] ;     Phe_Psi=df.at[i, 'Phe_Psi'] ;   Phe_Chi1=df.at[i, 'Phe_Chi1']
        Gly_Phi=df.at[i, 'Gly_Phi'] ;     Gly_Psi=df.at[i, 'Gly_Psi']
        APE10_Phi=df.at[i, 'APE10_Phi'] ; APE10_Psi=df.at[i, 'APE10_Psi']
        APE9_Phi=df.at[i, 'APE9_Phi'] ;   APE9_Psi=df.at[i, 'APE9_Psi']
        APE8_Phi=df.at[i, 'APE8_Phi'] ;   APE8_Psi=df.at[i, 'APE8_Psi'] ; APE8_Chi1=df.at[i, 'APE8_Chi1']
        APE7_Phi=df.at[i, 'APE7_Phi'] ;   APE7_Psi=df.at[i, 'APE7_Psi']
        APE6_Phi=df.at[i, 'APE6_Phi'] ;   APE6_Psi=df.at[i, 'APE6_Psi']

        # SB residues have to be {RK} and {E}, else SB-na
        KEtype=df.at[i, 'LYSres'] + df.at[i, 'GLUres']
        if KEtype == 'KE' or KEtype == 'RE' :
            KEcutoff=hbondcutoff
            if dis_saltbridge<=KEcutoff:
                df.at[i, 'Saltbr_label']='Saltbr-in'
            elif dis_saltbridge>0 and dis_saltbridge<900:
                df.at[i, 'Saltbr_label']='Saltbr-out'
            else:
                df.at[i, 'Saltbr_label']="Saltbr-none"
        else:
            df.at[i, 'Saltbr_label']="Saltbr-na"

        # HRD requirement; skipping exception for AATK, LMTK2, LMTK3 for now
        # HRD: phi (-180,0); psi (-100,0)
        # Arg: phi (0,180); psi (-100,50)

        if df.at[i,'aFaspres'] != "D":
            HRD_label = "HRD-na" + df.at[i,'aFaspres']
        elif max( HRD_Phi, HRD_Psi, Arg_Phi, Arg_Psi)>900: # none if disordered
            HRD_label="HRD-none"
        else:
            if ( dihe_in_range(-180,  0, HRD_Phi) and  # HRD_Phi in (-180,0)
                 dihe_in_range(-100,  0, HRD_Psi) and  # HRD_Psi in (-100,0)
                 dihe_in_range(   0,180, Arg_Phi) and  # Arg_Phi in (0,180)
                 dihe_in_range(-100, 50, Arg_Psi) ):   # Arg_Psi in (-100,50)
                HRD_label="HRD-in"
            else:
                HRD_label="HRD-out"
        df.at[i, 'HRD_label']=HRD_label
        
        # ActLoopNT calculation
        if dis_xhrd <= hbondcutoff:
            df.at[i, 'ActLoopNT_label']='ActLoopNT-in'
        elif dis_xhrd>0 and dis_xhrd<900:
            df.at[i, 'ActLoopNT_label']='ActLoopNT-out'
        else:
            df.at[i, 'ActLoopNT_label']="ActLoopNT-none"

        # Dihedral angle criteria on APE6-APE10
        # APR67 dihedral label
        # original criteria
        #      APE6: ϕ ∈ (-100°, 0°), ψ ∈ (-80°, 25°) and APE7: ϕ ∈ (-180°, 0°), ψ ∈ (-100°, 0°)
        #      OR, in case of peptide flip –
        #      APE6: ϕ ∈ (60°, 150°), ψ ∈ (-80°, 25°) and APE7: ϕ ∈ (-180°, 0°), ψ ∈ (100°, 180°)
        if APEtype=="noAPE":
            APE67_dihe_label="APE67-dihe-na" # not applicable
        elif max( APE7_Phi, APE7_Psi, APE6_Phi, APE6_Psi)>900: # none if disordered
            APE67_dihe_label="APE67-dihe-none"
        else:
            if ( dihe_in_range(-180,  0, APE7_Phi) and  # APE7_Phi in (-180,0)
                 dihe_in_range(-100, 50, APE7_Psi) and  # APE7_Psi in (-100,0)
                 dihe_in_range(-180,  0, APE6_Phi) and  # APE6_Phi in (-180,0)
                 dihe_in_range(-100, 50, APE6_Psi) ):  # APE6_Psi in (-80,25)
                APE67_dihe_label="APE67-dihe-in"
            elif ( dihe_in_range(-180,  0, APE7_Phi) and  # APE7_Phi in (-180,0)
                 dihe_in_range( 100, 180, APE7_Psi) and  # APE7_Psi in (-100,0)
                 dihe_in_range(  0, 180,  APE6_Phi) and  # APE6_Phi in (-180,0)
                 dihe_in_range( -50, 100, APE6_Psi) ):  # APE6_Psi in (-80,25)
                APE67_dihe_label="APE67-dihe-inflip"
            else:
                APE67_dihe_label="APE67-dihe-out"
        df.at[i, 'APE67_dihe_label']=APE67_dihe_label
        
        # APR8 dihedral label
        # APE8: ϕ ∈ (-180°, 0°), ψ ∈ (100°, 180°), χ1[Thr/Ser] ∈ (-120°, 0°) in non-TYR kinases,
        # APE8: ϕ ∈ (-100°, -62°), ψ ∈ (60°, 115°), χ1[Pro] ∈ (0°, 180°) in TYR kinases
        if APEtype=="noAPE":
            APE8_dihe_label="APE8-dihe-na" # not applicable
        elif max( APE8_Phi, APE8_Psi )>900: # none if disordered
            APE8_dihe_label="APE8-dihe-none"
        else:
            APE8_dihe_label="APE8-dihe-out"
            if (dihe_in_range(-180, 0,APE8_Phi) and
                dihe_in_range(50, 180,APE8_Psi)):
                APE8_dihe_label="APE8-dihe-in"
        df.at[i, 'APE8_dihe_label']=APE8_dihe_label

        # APE8 rotamer
        if APEtype=="noAPE" or APEtype=="TYR":
            APE8_rot_label="APE8-rot-na" # not applicable
        elif APE8_Chi1 >900: # none if disordered
            APE8_rot_label="APE8-rot-none"
        else:
            APE8_rot_label="APE8-rot-out"
            if df.at[i,'APE8res'] in ("S","T"):
                if (dihe_in_range(240,360,APE8_Chi1)):  APE8_rot_label="APE8-rot-in"
            else:
                APE8_rot_label="APE8-rot-na"
        df.at[i, 'APE8_rot_label']=APE8_rot_label

        if APEtype=="TYR":
            if max(APE9_Phi, APE9_Psi)>900:
                APE9_dihe_label="APE9-dihe-none"
            elif (dihe_in_range(-180, 0,APE9_Phi) and dihe_in_range(50, 180,APE9_Psi)):
                APE9_dihe_label="APE9-dihe-in"
            else:
                APE9_dihe_label="APE9-dihe-out"
        else:
            APE9_dihe_label="APE9-dihe-na"
        df.at[i, 'APE9_dihe_label']=APE9_dihe_label
            
        # APE10
        # earlier criteria
        #      ψ ∈ (80°, 180°) in TYR kinases
        if APEtype=="TYR":
            if max(APE10_Phi,APE10_Psi) > 900:
                APE10_dihe_label="APE10-dihe-none"
            elif dihe_in_range(-180, 0, APE10_Phi) and dihe_in_range(50, 180, APE10_Psi):
                APE10_dihe_label="APE10-dihe-in"
            else:
                APE10_dihe_label="APE10-dihe-out"
        else:
            APE10_dihe_label="APE10-dihe-na" # not applicable
        df.at[i, 'APE10_dihe_label']=APE10_dihe_label


        # Distance criteria on APE9-APE12

        # APE9: APE9-Cα/hRd-O distance < 6 Å in non-TYR kinases
        # APE9-Cα/hRd-O distance < 8 Å
        if group == "TYR":
            ape9cutoff=8.0
        else:
            ape9cutoff=6.0

        if APEtype=="TYR" or APEtype=="nonTYR":
            if dis_ape9>900:
                APE9_dist_label="APE9-dist-none"
            elif dis_ape9 <= ape9cutoff:
                APE9_dist_label="APE9-dist-in"
            else:
                APE9_dist_label="APE9-dist-out"
        else:
            APE9_dist_label="APE9-dist-na"
        df.at[i, 'APE9_dist_label']=APE9_dist_label
                

        # APE10
        # APE10: APE10-Cβ/DFG4-Cα: distance ∈ (8 Å, 14 Å) in non-TYR kinases
        if APEtype=="nonTYR":
            if dis_ape10>900:
                APE10_dist_label="APE10-dist-none"
            elif dis_ape10 <= 8.0:
                APE10_dist_label="APE10-dist-in"
            else:
                APE10_dist_label="APE10-dist-out"
        else:
            APE10_dist_label="APE10-dist-na"
        df.at[i, 'APE10_dist_label']=APE10_dist_label
        print("APE10",df.at[i,'Kincore_name'], df.at[i, 'pdbchain'], APEtype,dis_ape10,APE10_dist_label)

        # APE11
        # APE11: APE11-Cβ/DFG4-Cα: distance ∈ (8 Å, 14 Å) in non-TYR kinases
        if APEtype=="nonTYR":
            if dis_ape11>900:
                APE11_dist_label="APE11-dist-none"
            elif (dis_ape11 >= 8.0 and dis_ape11 <= 14):
                APE11_dist_label="APE11-dist-in"
            else:
                APE11_dist_label="APE11-dist-out"
        else:
            APE11_dist_label="APE11-dist-na"
        df.at[i, 'APE11_dist_label']=APE11_dist_label

        
        # APE12
        # APE12: APE12-Cβ/DFG4-Cα: distance ∈ (7 Å, 14 Å) in non-TYR kinases
        if APEtype=="nonTYR":
            if dis_ape12>900:
                APE12_dist_label="APE12-dist-none"
            elif (dis_ape12 >= 7.0 and dis_ape12 <= 14):
                APE12_dist_label="APE12-dist-in"
            else:
                APE12_dist_label="APE12-dist-out"
        else:
            APE12_dist_label="APE12-dist-na"
        df.at[i, 'APE12_dist_label']=APE12_dist_label


        APE_label=f'{APE10_dihe_label}_{APE9_dihe_label}_{APE8_dihe_label}_{APE8_rot_label}_{APE67_dihe_label}'
        APE_label=f'{APE_label}_{APE12_dist_label}_{APE11_dist_label}_{APE10_dist_label}_{APE9_dist_label}'
        (apedihe, apedist)=make_apelabel(APE_label)
        df.at[i, 'APE_label']=APE_label
        df.at[i, 'APEdihe_label']=apedihe
        df.at[i, 'APEdist_label']=apedist

        ActLoopCT_label="ActLoopCT-in"
        if APEtype=="noAPE":
            ActLoopCT_label="ActLoopCT-na"
        else:
            if "none" in APE_label: ActLoopCT_label="ActLoopCT-none"
            if "out" in APE_label:  ActLoopCT_label="ActLoopCT-out"   # if any component is "out", then "out"

        df.at[i, 'ActLoopCT_label']=ActLoopCT_label
        print(df.at[i,'Saltbr_label'], df.at[i, 'ActLoopNT_label'], df.at[i, 'ActLoopCT_label'])

        # Determine activity state: active kinase must be DFGin-BLAminus-Actloop-in with KEsb<= 3.6 Angstroms
        label0=df.at[i,'Chelix_label']
        (tmp1,label1)=df.at[i,'Saltbr_label'].split("-")
        (tmp2,label2)=df.at[i,'ActLoopNT_label'].split("-")
        (tmp3,label3)=df.at[i,'ActLoopCT_label'].split("-")

        SNClabel=f'{label1}-{label2}-{label3}'
        if label1=="na": onelabel1="a"
        else: onelabel1=label1[0:1]
        if label2=="na": onelabel2="a"
        else: onelabel2=label2[0:1]
        if label3=="na": onelabel3="a"
        else: onelabel3=label3[0:1]
        
        SNCpymollabel=f'SNC{onelabel1}{onelabel2}{onelabel3}'
        df.at[i,'ChelixSaltbr_label']=f'{label0}-{label1}'
        df.at[i,'ActLoop_label']=f'{label2}-{label3}'
        df.at[i,'SNC_label']=SNClabel
        df.at[i,'SNCpymol_label']=SNCpymollabel

        ActivityState="Active"
        if df.at[i, 'Spatial_label'] != "DFGin":     ActivityState='Inactive' # must be DFGin not none
        if df.at[i, 'Dihedral_label'] != "BLAminus": ActivityState='Inactive' # must be BLAminus not none
        if group=="PKDCC" and df.at[i, 'Dihedral_label'] == "ABAminus": ActivityState="Active" # exception
        if df.at[i, 'HRD_label'] == "HRD-out" :      ActivityState='Inactive' # must be HRD-in or HRD-na
        if df.at[i, 'HRD_label'] == "HRD-none" :     ActivityState='Inactive' # must be HRD-in or HRD-na
        if "in" not in df.at[i,'ChelixSaltbr_label'] : ActivityState='Inactive' # must be at least Chelix-in
        
        # once inactive, then always inactive; but any "out" criterion makes kinase Inactive
        if ActivityState in ("Active","None") and "none" in df.at[i, 'Saltbr_label'] : ActivityState="None"
        if ActivityState in ("Active","None") and "out"  in df.at[i, 'Saltbr_label'] : ActivityState="Inactive"
        if ActivityState in ("Active","None") and "none" in df.at[i, 'ActLoopNT_label'] : ActivityState="None"
        if ActivityState in ("Active","None") and "out"  in df.at[i, 'ActLoopNT_label'] : ActivityState="Inactive"
        if ActivityState in ("Active","None") and "none" in df.at[i, 'ActLoopCT_label'] : ActivityState="None"
        if ActivityState in ("Active","None") and "out"  in df.at[i, 'ActLoopCT_label'] : ActivityState="Inactive"

        dihedral=df.at[i,'Dihedral_label']
        if dihedral is None:
            dihedral="None"
            df.at[i,'Dihedral_label']=dihedral

        if df.at[i,'Status']=='Pseudo': ActivityState="Pseudo"
        df.at[i,'Activity_label']=ActivityState
        apelabel=df.at[i,'APE_label']
        apelabel=apelabel.replace("_"," ")
        string=""
        string += f'{df.at[i,"APEtype"]} '
        string += f'{df.at[i,"Kincore_name"]} '
        string += f'{df.at[i,"pdbchain"]} '
        string += f'{df.at[i,"Activity_label"]} '
        string += f'{df.at[i,"Spatial_label"]} '
        string += f'{df.at[i,"Dihedral_label"]} '
        string += f'{df.at[i,"Chelix_label"]} '
        string += f'{df.at[i,"ChelixSaltbr_label"]} '
        string += f'{df.at[i,"HRD_label"]} '
        string += f'{df.at[i,"ActLoopNT_label"]} '
        string += f'{df.at[i,"ActLoopCT_label"]} '
        string += f'{df.at[i,"ActLoop_label"]} '
        string += f'{df.at[i,"SNC_label"]} '
        string += f'{df.at[i,"SNCpymol_label"]} '
        string += f'{df.at[i,"APEdihe_label"]} '
        string += f'{df.at[i,"APEdist_label"]} '
        string += f'{apelabel} '
        string += f'HRD {HRD_Phi} {HRD_Psi} '
        string += f'ARG {Arg_Phi} {Arg_Psi} '
        string += f'APE10 {APE10_Phi} {APE10_Psi} '
        string += f'APE9 {APE9_Phi} {APE9_Psi} '
        string += f'APE8 {APE8_Phi} {APE8_Psi} {APE8_Chi1} '
        string += f'APE7 {APE7_Phi} {APE7_Psi} '
        string += f'APE6 {APE6_Phi} {APE6_Psi} '
        string += f'SB {dis_saltbridge}  '
        string += f'XHRD {dis_xhrd} '
        string += f'APEDist {dis_ape12} '
        string += f'{dis_ape11} '
        string += f'{dis_ape10} '
        string += f'{dis_ape9} '
        
        OUTPUT.write(string + "\n")
        print(string)
    df=df.copy()
    return df



def spatial_labels(pwd,df):
    print('Assigning spatial labels...')
    OUTPUT=open(f'{pwd}/Spatial_labels.txt','w')
    for i in df.index:
      
        df.at[i,'Spatial_label']='None'
        dis_phe_glu=float(df.at[i,'Phe_Glu4_dis'])
        dis_phe_lys=float(df.at[i,'Phe_Lys_dis'])

        if dis_phe_glu<=11 and dis_phe_lys>=11 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial_label']='DFGin'

        elif dis_phe_glu>11 and dis_phe_lys<=14 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial_label']='DFGout'

        elif dis_phe_glu<=11 and dis_phe_lys<=11 and dis_phe_glu!=999 and dis_phe_lys!=999:
            df.at[i,'Spatial_label']='DFGinter'

            OUTPUT.write(f"{i} {df.at[i,'pdbchain']} {df.at[i,'Spatial_label']}\n")

    df=df.copy()
    return df

def cosine_dis_without_chi1(df,i,spatial,cutoff):
    xdfg_phi=float(df.at[i,'XDFG_Phi'])
    xdfg_psi=float(df.at[i,'XDFG_Psi'])
    dfg_phi=float(df.at[i,'DFG_Phi'])
    dfg_psi=float(df.at[i,'DFG_Psi'])
    phe_phi=float(df.at[i,'Phe_Phi'])
    phe_psi=float(df.at[i,'Phe_Psi'])
    min_spatial=999
    
    for clusters in spatial:    
        cosine_dis=(2/6)*((1-math.cos(math.radians(xdfg_phi-float(spatial[clusters][0]))))+(1-math.cos(math.radians(xdfg_psi-float(spatial[clusters][1]))))+\
                (1-math.cos(math.radians(dfg_phi-float(spatial[clusters][2]))))+(1-math.cos(math.radians(dfg_psi-float(spatial[clusters][3]))))+\
                (1-math.cos(math.radians(phe_phi-float(spatial[clusters][4]))))+(1-math.cos(math.radians(phe_psi-float(spatial[clusters][5])))))
        
        if cosine_dis<=min_spatial:
            df.at[i,'Dihedral_dis_NoChi1']=np.round(cosine_dis,2)
            min_spatial=cosine_dis
            
            if cosine_dis<=float(cutoff):
                return clusters      #Only Dihedral column name is used for final labeling without chi1
    
def dihedral_labels(pwd, df,cutoff):
    print('Assigning dihedral labels...')
    
    dfginter={'BABtrans':(-80.20,128.22,-117.47,23.76,-85.16,133.21,181.42)}
    dfgout={'BBAminus':(-138.56,-176.12,-144.35,103.66,-82.59,-9.03,290.59)}
    
    dfgin_minus={'BLAminus':(-128.64,178.67,61.15,81.21,-96.89,20.53,289.12),\
       'ABAminus':(-111.82,-7.64,-141.55,148.01,-127.79,23.32,296.17),\
       'BLBminus':(-134.79,175.48,60.44,65.35,-79.44,145.34,287.56)}
    
    dfgin_plus={'BLAplus':(-119.24,167.71,58.94,34.08,-89.42,-8.54,55.63),\
       'BLBplus':(-125.28,172.53,59.98,32.92,-85.51,145.28,49.01)}
    
    dfgin_trans={'BLBtrans':(-106.16,157.24,69.37,21.33,-61.73,134.56,215.23)}
    
    OUTPUT=open(f'{pwd}/Dihedral_labels.txt','w')
    for i in df.index:
        xdfg_phi=float(df.at[i,'XDFG_Phi']);xdfg_psi=float(df.at[i,'XDFG_Psi']);dfg_phi=float(df.at[i,'DFG_Phi']);
        dfg_psi=float(df.at[i,'DFG_Psi']);phe_phi=float(df.at[i,'Phe_Phi'])
        phe_psi=float(df.at[i,'Phe_Psi']);phe_chi1=float(df.at[i,'Phe_Chi1'])
        
        df.at[i,'Dihedral_label']='None'    #Default label is None
        df.at[i,'Dihedral_dis_NoChi1']=999
        
        if (xdfg_phi==999 or xdfg_psi==999 or dfg_phi==999 or dfg_psi==999 or phe_phi==999 or  
            phe_psi==999 or phe_chi1==999):
            df.at[i,'Dihedral_label']='None'
            df.at[i,'Dihedral_dis_NoChi1']=999
            OUTPUT.write(f"{i} {df.at[i,'pdbchain']} {df.at[i,'Dihedral_label']}\n")
            continue
        
        if df.at[i, 'Spatial_label']=='DFGin':     
            if phe_chi1 > 240 and phe_chi1 <= 360:
                df.at[i, 'Dihedral_label']=cosine_dis_without_chi1(df, i, dfgin_minus, cutoff)
            if phe_chi1>0 and phe_chi1<=120:
                df.at[i, 'Dihedral_label']=cosine_dis_without_chi1(df, i, dfgin_plus, cutoff)
            if phe_chi1>120 and phe_chi1<=240:
                df.at[i, 'Dihedral_label']=cosine_dis_without_chi1(df, i, dfgin_trans, cutoff)
            
        if df.at[i, 'Spatial_label']=='DFGinter': 
            if phe_chi1>120 and phe_chi1<=240:
                df.at[i, 'Dihedral_label']=cosine_dis_without_chi1(df, i, dfginter, cutoff)
        
        if df.at[i, 'Spatial_label']=='DFGout':        
            if phe_chi1>240 and phe_chi1<=360:
                df.at[i, 'Dihedral_label']=cosine_dis_without_chi1(df, i, dfgout, cutoff)

        OUTPUT.write(f"{i} {df.at[i, 'pdbchain']} {df.at[i,'Dihedral_label']}\n")
    df=df.copy()
    return df
