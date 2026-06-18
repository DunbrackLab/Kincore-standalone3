#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 12:27:13 2020
Edited on Wed Aug 30, 2023 by RLD
Last edited on Thu June 18, 2026 by JG
@author: vivekmodi
"""

def assign_default_values(index,conf_df):
    conf_df.at[index,"APEtype"] = 'None'
    conf_df.at[index,"ChelixSaltbr_label"] = 'None'
    conf_df.at[index,"HRD_label"] = 'None'
    conf_df.at[index,"APEdihe_label"] = 'None'
    conf_df.at[index,"APEdist_label"] = 'None'

    conf_df.at[index, 'Actloop_phos']='Unphos'
    conf_df.at[index,'Activity_label']='None'
    conf_df.at[index,'ActLoopNT_label']='None'
    conf_df.at[index,'ActLoopCT_label']='None'
    conf_df.at[index,'APE_label']='None' # Not used anywhere currently
    conf_df.at[index,'APEdihe_label']='None'
    conf_df.at[index,'APEdist_label']='None'

    conf_df.at[index, 'HRD_label'] = 'None'
    conf_df.at[index, 'vActLoopCT_label'] = 'XXXXXXXXXAPE' # "Verbose ActLoopCT label" not used anywhere currently
    conf_df.at[index,'Spine_label']='None'
    conf_df.at[index,'Saltbr_label']='None'
    conf_df.at[index,'Chelix_label']='None'
    conf_df.at[index,'Spatial_label']='None'
    conf_df.at[index,'Dihedral_label']='None'
    conf_df.at[index,'Ligand']='None'
    conf_df.at[index,'Ligand_label']='None'

    conf_df.at[index,'Sequence']='-'
    conf_df.at[index,'First_res']=0
    conf_df.at[index,'Group']='None'
    conf_df.at[index,'hmm']='None'

    conf_df.at[index,'APE3-aFasp3-dis']=999 #Not used anywhere currently 
    conf_df.at[index,'APE10-DFG4-dis']=999 #Added 11/20/2024 - JG
    conf_df.at[index,'APE11-DFG4-dis']=999 #Added 11/20/2024 - JG
    conf_df.at[index,'APE12-DFG4-dis']=999 #Added 11/20/2024 - JG
    conf_df.at[index,'DFG6-XHRD-dis']=999
    conf_df.at[index,'Spine1-dis']=999
    conf_df.at[index,'Spine2-dis']=999
    conf_df.at[index,'Spine3-dis']=999
    conf_df.at[index,'Spine-dis']=999
    conf_df.at[index,'Glu4-Phe-dis']=999
    conf_df.at[index,'APE9-Arg-dis']=999
    conf_df.at[index,'Lys-Phe-dis']=999
    conf_df.at[index,'Lys-Glu-dis']=999
    conf_df.at[index,'LysNZ-GluOE-dis']=999

    conf_df.at[index,'aFasp3_num']=0 #Not used anywhere currently 
    conf_df.at[index,'real_aFasp3_num']=0
    conf_df.at[index,'aFasp3_restype']='-'

    conf_df.at[index,'aFasp_num']=0 #Added 07/10/2024 - JG
    conf_df.at[index,'real_aFasp_num']=0
    conf_df.at[index,'aFasp_restype']='-'

    conf_df.at[index,'HRDAsp_num']=0
    conf_df.at[index,'real_HRDAsp_num']=0
    conf_df.at[index,'HRDAsp_restype']='-'

    conf_df.at[index,'Lys_num']=0
    conf_df.at[index,'real_Lys_num']=0
    conf_df.at[index,'Lys_restype']='-'

    conf_df.at[index,'Glu_num']=0
    conf_df.at[index,'real_Glu_num']=0
    conf_df.at[index,'Glu_restype']='-'

    conf_df.at[index,'Glu4_num']=0
    conf_df.at[index,'real_Glu4_num']=0
    conf_df.at[index,'Glu4_restype']='-'

    conf_df.at[index,'APE_num']=0
    conf_df.at[index,'real_APE_num']=0
    conf_df.at[index,'APE_restype']='-'

    conf_df.at[index,'APE3_num']=0
    conf_df.at[index,'real_APE3_num']=0
    conf_df.at[index,'APE3_restype']='-'

    conf_df.at[index,'APE6_num']=0
    conf_df.at[index,'real_APE6_num']=0
    conf_df.at[index,'APE6_restype']='-'

    conf_df.at[index,'APE7_num']=0
    conf_df.at[index,'real_APE7_num']=0
    conf_df.at[index,'APE7_restype']='-'

    conf_df.at[index,'APE8_num']=0
    conf_df.at[index,'real_APE8_num']=0
    conf_df.at[index,'APE8_restype']='-'

    conf_df.at[index,'APE9_num']=0
    conf_df.at[index,'real_APE9_num']=0
    conf_df.at[index,'APE9_restype']='-'

    conf_df.at[index,'APE10_num']=0
    conf_df.at[index,'real_APE10_num']=0
    conf_df.at[index,'APE10_restype']='-'

    conf_df.at[index,'APE11_num']=0
    conf_df.at[index,'real_APE11_num']=0
    conf_df.at[index,'APE11_restype']='-'

    conf_df.at[index,'APE12_num']=0
    conf_df.at[index,'real_APE12_num']=0
    conf_df.at[index,'APE12_restype']='-'

    conf_df.at[index,'APEp2_num']=0
    conf_df.at[index,'real_APEp2_num']=0
    conf_df.at[index,'APEp2_restype']='-'

    conf_df.at[index,'APEp3_num']=0
    conf_df.at[index,'real_APEp3_num']=0
    conf_df.at[index,'APEp3_restype']='-'

    conf_df.at[index,'DFG6_num']=0
    conf_df.at[index,'real_DFG6_num']=0
    conf_df.at[index,'DFG6_restype']='-'

    conf_df.at[index,'HPN7_num']=0   # for spine calculation
    conf_df.at[index,'real_HPN7_num']=0
    conf_df.at[index,'HPN7_restype']='-'

    conf_df.at[index,'HRD_num']=0
    conf_df.at[index,'real_HRD_num']=0
    conf_df.at[index,'HRD_restype']='-'

    conf_df.at[index,'Arg_num']=0
    conf_df.at[index,'real_Arg_num']=0
    conf_df.at[index,'Arg_restype']='-'

    conf_df.at[index,'XHRD_num']=0
    conf_df.at[index,'real_XHRD_num']=0
    conf_df.at[index,'XHRD_restype']='-'

    conf_df.at[index,'XDFG_num']=0
    conf_df.at[index,'real_XDFG_num']=0
    conf_df.at[index,'XDFG_restype']='-'

    conf_df.at[index,'Asp_num']=0
    conf_df.at[index,'real_Asp_num']=0
    conf_df.at[index,'Asp_restype']='-'

    conf_df.at[index,'Phe_num']=0
    conf_df.at[index,'real_Phe_num']=0
    conf_df.at[index,'Phe_restype']='-'

    conf_df.at[index,'Gly_num']=0
    conf_df.at[index,'real_Gly_num']=0
    conf_df.at[index,'Gly_restype']='-'

    conf_df.at[index,'DFG4_num']=0
    conf_df.at[index,'real_DFG4_num']=0
    conf_df.at[index,'DFG4_restype']='-'

    conf_df.at[index,'Hinge1_num']=0
    conf_df.at[index, 'Score']=0

    conf_df.at[index,'XHRD_Phi']=999
    conf_df.at[index,'XHRD_Psi']=999
    conf_df.at[index,'HRD_Phi']=999
    conf_df.at[index,'HRD_Psi']=999
    conf_df.at[index,'Arg_Phi']=999
    conf_df.at[index,'Arg_Psi']=999

    conf_df.at[index,'XDFG_Phi']=999
    conf_df.at[index,'XDFG_Psi']=999
    conf_df.at[index,'Asp_Phi']=999
    conf_df.at[index,'Asp_Psi']=999
    conf_df.at[index,'Asp_Chi1']=999
    conf_df.at[index,'Asp_Chi2']=999
    conf_df.at[index,'Phe_Phi']=999
    conf_df.at[index,'Phe_Psi']=999
    conf_df.at[index,'Phe_Chi1']=999
    conf_df.at[index,'Phe_Chi2']=999
    conf_df.at[index,'Gly_Phi']=999
    conf_df.at[index,'Gly_Psi']=999
    conf_df.at[index,'APE6_Phi']=999
    conf_df.at[index,'APE6_Psi']=999
    conf_df.at[index,'APE7_Phi']=999
    conf_df.at[index,'APE7_Psi']=999
    conf_df.at[index,'APE8_Phi']=999
    conf_df.at[index,'APE8_Psi']=999
    conf_df.at[index,'APE8_Chi1']=999
    conf_df.at[index,'APE9_Phi']=999
    conf_df.at[index,'APE9_Psi']=999
    conf_df.at[index,'APE10_Phi']=999
    conf_df.at[index,'APE10_Psi']=999

    return conf_df
