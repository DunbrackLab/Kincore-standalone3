#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Last edited on Wed Aug 30, 2023 by RLD                                                                                                                                                                                                                                    

#def chelix_conformation(index,conf_df):
#    dis_sb=conf_df.at[index,'Lys-Glu-dis']
#    
#    if dis_sb <= 10.0:
#        conf_df.at[index,'Chelix_label']='Chelix-in'
#    elif dis_sb<900:
#        conf_df.at[index,'Chelix_label']='Chelix-out'
#    else:
#        conf_df.at[index,'Chelix_label']='None'
#    return conf_df

def chelix_conformation(index,df):
    if float(df.at[index,'Lys-Glu-dis'])<=10.0:
        df.at[index,'Chelix_label']='Chelix-in'
    if float(df.at[index,'Lys-Glu-dis'])>10.0:
        df.at[index,'Chelix_label']='Chelix-out'
    if float(df.at[index,'Lys-Glu-dis'])==999:
        df.at[index,'Chelix_label']='Chelix-none'
    if df.at[index, 'hmm']=='WNK':
        df.at[index, 'Chelix_label'] = 'Chelix-na'
    return df
