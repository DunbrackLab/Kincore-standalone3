#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 21:13:13 2020

@author: vivekmodi
added 2025 crieria for active kinases. RLD 12/12/2025
Last edited on Thu June 18, 2026 by JG
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
    
def make_apedihelabel(apelabel):
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
    return apedihe

def make_apedistlabel(apelabel):    
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
    return apedist


def active_labels(index, df):
    hbondcutoff=3.6
    # APEtype: TYR, nonTYR, noAPE
    group=df.at[index, 'Group']
    if group == "TYR": APEtype="TYR"
    else: APEtype="nonTYR"
    # BUB, HASP, PKDCC, TP53RK exception, no APE motif; check
    # PEAK3 has group=PEAK and score~170 and no APE motif (PEAK1 and PRAG1, >500)
    #            if (group=="BUB" or group=="HASP" or group=="PKDCC" or
    if (group=="HASP" or group=="PKDCC" or
        group=="TP53RK" or group=="PAN3" or group=="RNASEL" or
        (group=="PEAK" and df.at[index, 'Score']<300)): APEtype="noAPE"
    df.at[index, 'APEtype']=APEtype
    
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
        if df.at[index, 'Score']<180.0:
            finalgroup="OTHER"  # some OTHER kinases have higher scores for CAMK, CMGC, etc than the OTHER HMM

    # structural parameters for active state determination
    dis_saltbridge=float(df.at[index,'LysNZ-GluOE-dis'])
    dis_xhrd=float(df.at[index,'DFG6-XHRD-dis'])
    dis_ape9=float(df.at[index,'APE9-Arg-dis'])
    dis_ape10=df.at[index, 'APE10-DFG4-dis']
    dis_ape11=df.at[index, 'APE11-DFG4-dis']
    dis_ape12=df.at[index, 'APE12-DFG4-dis']
    HRD_Phi=df.at[index, 'HRD_Phi'] ;     HRD_Psi=df.at[index, 'HRD_Psi']
    Arg_Phi=df.at[index, 'Arg_Phi'] ;     Arg_Psi=df.at[index, 'Arg_Psi']
    Asp_Chi1=df.at[index, 'Asp_Chi1']
    Phe_Phi=df.at[index, 'Phe_Phi'] ;     Phe_Psi=df.at[index, 'Phe_Psi'] ;   Phe_Chi1=df.at[index, 'Phe_Chi1']
    Gly_Phi=df.at[index, 'Gly_Phi'] ;     Gly_Psi=df.at[index, 'Gly_Psi']
    APE10_Phi=df.at[index, 'APE10_Phi'] ; APE10_Psi=df.at[index, 'APE10_Psi']
    APE9_Phi=df.at[index, 'APE9_Phi'] ;   APE9_Psi=df.at[index, 'APE9_Psi']
    APE8_Phi=df.at[index, 'APE8_Phi'] ;   APE8_Psi=df.at[index, 'APE8_Psi'] ; APE8_Chi1=df.at[index, 'APE8_Chi1']
    APE7_Phi=df.at[index, 'APE7_Phi'] ;   APE7_Psi=df.at[index, 'APE7_Psi']
    APE6_Phi=df.at[index, 'APE6_Phi'] ;   APE6_Psi=df.at[index, 'APE6_Psi']

    # SB residues have to be {RK} and {E}, else SB-na
    KEtype=df.at[index, 'Lys_restype'] + df.at[index, 'Glu_restype']
    if KEtype == 'KE' or KEtype == 'RE' :
        KEcutoff=hbondcutoff
        if dis_saltbridge<=KEcutoff:
            df.at[index, 'Saltbr_label']='Saltbr-in'
        elif dis_saltbridge>0 and dis_saltbridge<900:
            df.at[index, 'Saltbr_label']='Saltbr-out'
        else:
            df.at[index, 'Saltbr_label']="Saltbr-none"
    else:
        df.at[index, 'Saltbr_label']="Saltbr-na"

    # HRD requirement; skipping exception for AATK, LMTK2, LMTK3 for now
    # HRD: phi (-180,0); psi (-100,0)
    # Arg: phi (0,180); psi (-100,50)
    if df.at[index,'aFasp_restype'] != "D":
        HRD_label = "HRD-na" + df.at[index,'aFasp_restype']
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
    df.at[index, 'HRD_label']=HRD_label
    
    # DFG Asp rotamer
    if (df.at[index,'Asp_restype'] != "D"):
        Asp_rot_label="DFGAsp-rot-na"
    elif (Asp_Chi1 > 900):
        Asp_rot_label="DFGAsp-rot-none"
    else:
        Asp_rot_label="DFGAsp-rot-out"
        if (dihe_in_range(120,240,Asp_Chi1)):  Asp_rot_label="DFGAsp-rot-in"
    df.at[index, 'Asp_rot_label']=Asp_rot_label

    # ActLoopNT calculation
    if dis_xhrd <= hbondcutoff:
        df.at[index, 'ActLoopNT_label']='ActLoopNT-in'
    elif dis_xhrd>0 and dis_xhrd<900:
        df.at[index, 'ActLoopNT_label']='ActLoopNT-out'
    else:
        df.at[index, 'ActLoopNT_label']="ActLoopNT-none"

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
    df.at[index, 'APE67_dihe_label']=APE67_dihe_label
    
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
    df.at[index, 'APE8_dihe_label']=APE8_dihe_label

    # APE8 rotamer
    if APEtype=="noAPE" or APEtype=="TYR":
        APE8_rot_label="APE8-rot-na" # not applicable
    elif APE8_Chi1 >900: # none if disordered
        APE8_rot_label="APE8-rot-none"
    else:
        APE8_rot_label="APE8-rot-out"
        if df.at[index,'APE8_restype'] in ("S","T"):
            if (dihe_in_range(240,360,APE8_Chi1)):  APE8_rot_label="APE8-rot-in"
        else:
            APE8_rot_label="APE8-rot-na"
    df.at[index, 'APE8_rot_label']=APE8_rot_label

    if APEtype=="TYR":
        if max(APE9_Phi, APE9_Psi)>900:
            APE9_dihe_label="APE9-dihe-none"
        elif (dihe_in_range(-180, 0,APE9_Phi) and dihe_in_range(50, 180,APE9_Psi)):
            APE9_dihe_label="APE9-dihe-in"
        else:
            APE9_dihe_label="APE9-dihe-out"
    else:
        APE9_dihe_label="APE9-dihe-na"
    df.at[index, 'APE9_dihe_label']=APE9_dihe_label
        
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
    df.at[index, 'APE10_dihe_label']=APE10_dihe_label


    # Distance criteria on APE9-APE12

    # APE9
    # APE9: APE9-Cα/hRd-O distance < 6 Å in non-TYR kinases
    if APEtype=="nonTYR":
        if dis_ape9>900:
            APE9_dist_label="APE9-dist-none"
        elif dis_ape9 <= 6.0:
            APE9_dist_label="APE9-dist-in"
        else:
            APE9_dist_label="APE9-dist-out"
    else:
        APE9_dist_label="APE9-dist-na"
    df.at[index, 'APE9_dist_label']=APE9_dist_label
            

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
    df.at[index, 'APE10_dist_label']=APE10_dist_label

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
    df.at[index, 'APE11_dist_label']=APE11_dist_label

    
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
    df.at[index, 'APE12_dist_label']=APE12_dist_label

    APEdihe_label=f'{APE10_dihe_label} {APE9_dihe_label} {APE8_dihe_label} {APE8_rot_label} {APE67_dihe_label}'
    APEdist_label=f'{APE12_dist_label} {APE11_dist_label} {APE10_dist_label} {APE9_dist_label}'
    apedihe=make_apedihelabel(APEdihe_label)
    apedist=make_apedistlabel(APEdist_label)
    
    APE_label = f'{APEdihe_label} {APEdist_label}'

    df.at[index, 'APE_label']=APE_label
    df.at[index, 'APEdihe_label']=APEdihe_label
    df.at[index, 'APEdist_label']=APEdist_label

    ActLoopCT_label="ActLoopCT-in"
    if APEtype=="noAPE":
        ActLoopCT_label="ActLoopCT-na"
    else:
        if "none" in APE_label: ActLoopCT_label="ActLoopCT-none"
        if "out" in APE_label:  ActLoopCT_label="ActLoopCT-out"   # if any component is "out", then "out"

    df.at[index, 'ActLoopCT_label']=ActLoopCT_label
# Determine activity state: active kinase must be DFGin-BLAminus-Actloop-in with KEsb<= 3.6 Angstroms
    label0=df.at[index,'Chelix_label']
    (tmp1,label1)=df.at[index,'Saltbr_label'].split("-")
    (tmp2,label2)=df.at[index,'ActLoopNT_label'].split("-")
    (tmp3,label3)=df.at[index,'ActLoopCT_label'].split("-")

    SNClabel=f'{label1}-{label2}-{label3}'
    if label1=="na": onelabel1="a"
    else: onelabel1=label1[0:1]
    if label2=="na": onelabel2="a"
    else: onelabel2=label2[0:1]
    if label3=="na": onelabel3="a"
    else: onelabel3=label3[0:1]

    df.at[index,'ChelixSaltbr_label']=f'{label0}-{label1}'

    ActivityState="Active"
    if df.at[index, 'Spatial_label'] != "DFGin":     ActivityState='Inactive' # must be DFGin not none
    if df.at[index, 'Dihedral_label'] != "BLAminus": ActivityState='Inactive' # must be BLAminus not none
    if df.at[index, 'Asp_rot_label'] != "DFGAsp-rot-in": ActivityState='Inactive'
    if group=="PKDCC" and df.at[index, 'Dihedral_label'] == "ABAminus": ActivityState="Active" # exception
    if df.at[index, 'HRD_label'] == "HRD-out" :      ActivityState='Inactive' # must be HRD-in or HRD-na
    if df.at[index, 'HRD_label'] == "HRD-none" :     ActivityState='Inactive' # must be HRD-in or HRD-na
    if ("in" not in df.at[index,'ChelixSaltbr_label']) and ("na" not in df.at[index,'ChelixSaltbr_label']) : ActivityState='Inactive' # must be at least Chelix-in
    
    # once inactive, then always inactive; but any "out" criterion makes kinase Inactive
    if ActivityState in ("Active","None") and "none" in df.at[index, 'Saltbr_label'] : ActivityState="None"
    if ActivityState in ("Active","None") and "out"  in df.at[index, 'Saltbr_label'] : ActivityState="Inactive"
    if ActivityState in ("Active","None") and "none" in df.at[index, 'ActLoopNT_label'] : ActivityState="None"
    if ActivityState in ("Active","None") and "out"  in df.at[index, 'ActLoopNT_label'] : ActivityState="Inactive"
    if ActivityState in ("Active","None") and "none" in df.at[index, 'ActLoopCT_label'] : ActivityState="None"
    if ActivityState in ("Active","None") and "out"  in df.at[index, 'ActLoopCT_label'] : ActivityState="Inactive"

    dihedral=df.at[index,'Dihedral_label']
    if dihedral is None:
        dihedral="None"
        df.at[index,'Dihedral_label']=dihedral

    # for MOS and PINK1 where Glu4=ALA
    if df.at[index,'Glu4_restype']=="A":
        df.at[index,'Spine-dis']=df.at[index,'Spine1-dis']
                
    if df.at[index,'Spine-dis']<0:
        df.at[index,'Spine_label']='Spine-none'
    elif df.at[index,'Spine-dis']<=5.0:
        df.at[index,'Spine_label']="Spine-in"
    elif df.at[index,'Spine-dis']<900:
        df.at[index,'Spine_label']="Spine-out"
    else:
        df.at[index,'Spine_label']="Spine-none"

#    if df.at[index,'Status']=='Pseudo': ActivityState="Pseudo" #Not yet implemented
    df.at[index,'Activity_label']=ActivityState        
    return df
