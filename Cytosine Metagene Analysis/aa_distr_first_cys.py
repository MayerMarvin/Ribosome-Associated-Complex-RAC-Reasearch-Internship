# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:42:14 2022

@author: user

Metaanalysis Koji

Identify amino acid distribution around window of the first cytosin in yeast proteome
"""


"""
IMPORT AND LOAD PROTEOME
"""

import pickle
import pandas as pd
import plotly.express as px
import plotly.io as pio
pio.renderers.default='browser'

#what about cases where Cystein ist at front/ back border of sequence?
def load_window_into_aa_dict(window, window_aa):
    for pos, char in enumerate(window):
        if char =='*': # stop @ stop codon and dont display it
            break

        if pos not in window_aa.keys():
            window_aa[pos] = {}
            window_aa[pos][char] = 1
        elif pos in window_aa.keys():
            if char not in window_aa[pos].keys():
                window_aa[pos][char] = 1
            elif char in window_aa[pos].keys():
                window_aa[pos][char] += 1

    return window_aa

def import_data(path_to_file, window_size):
    with open(path_to_file, 'rb') as f:
        proteome = pickle.load(f)
        
    window_aa = {}

    for key, sequence in proteome.items():   
        aa=""
        for pos, aa in enumerate(sequence): 
            if aa == "C": 
                window = sequence[pos-window_size:pos+window_size+1]
                window_aa = load_window_into_aa_dict(window, window_aa)              
    
                break #only first cystein
    
    return proteome, window_aa



"""
NORMALIZE TO 1
"""

def normalize_dict(aa_dict):
    for pos in aa_dict.keys():
        total_count = 0
        
        for aa in aa_dict[pos].keys():
            total_count += aa_dict[pos][aa]

        for aa in aa_dict[pos].keys():
            aa_dict[pos][aa] = aa_dict[pos][aa] / total_count

    return aa_dict
    



"""
GET SIGNIFIKANT VALUES -> change to 0 if not
"""

def check_if_in_range(occurence, significant_percentage, normalized_value):
    lower = normalized_value - significant_percentage
    upper = normalized_value + significant_percentage
    #lower = normalized_value * (1 - significant_percentage)
    #upper = normalized_value * (1 + significant_percentage)
    if lower <= occurence <= upper:
        return True
    else:
        return False 
    
# 1.) get normal distribution of all aa in genome
def get_proteome_aa_distribution(sequence_data):
    normal_distribution = {}  
    total_aa_count = 0
    for key, sequence in sequence_data.items():   
        for aa in sequence: 
            if aa != "*":
                if aa not in normal_distribution:
                    normal_distribution[aa] = 1
                    total_aa_count += 1
                else:
                    normal_distribution[aa] += 1
                    total_aa_count += 1

    for aa in normal_distribution.keys():
        normal_distribution[aa] = normal_distribution[aa] / total_aa_count

    return normal_distribution



def check_for_significant_difference (normal_distribution, aa_dict, range_to_be_tested):
    for pos in aa_dict.keys():
        for aa in aa_dict[pos].keys():
            if check_if_in_range(aa_dict[pos][aa], range_to_be_tested, normal_distribution[aa]):
                aa_dict[pos][aa] = 0
    return aa_dict
 

def get_significant_aa(proteome, significant_percentage, normalized_aa_window):     
    proteome_aa_distribution = get_proteome_aa_distribution(proteome)    
    significant_aa_window = check_for_significant_difference(proteome_aa_distribution, 
                                                             normalized_aa_window, 
                                                             significant_percentage) 
    
    plotly_df = pd.DataFrame.from_dict(significant_aa_window, orient='index')
    plotly_df.index.name = 'position'
    return plotly_df


"""
COLOR ACCORDING TO AMINO ACID GROUP
"""

colors_aa_by_group= {
    # positively charged --> orange
    'R':'#FF8C00',
    'H':'#FFA500',
    'K':'#FF7F50',
    # negatively charged --> yellow
    'D':'#FFFF00',
    'E':'#FFD700',
    # polar, uncharged --> red
    'S':'#DC143C',
    'T':'#FF0000',
    'N':'#B22222',
    'Q':'#8B0000',
    # special cases --> green
    'C':'#32CD32',
    'U':'#00FF00',
    'G':'#9ACD32',
    'P':'#ADFF2F',
    # hydrophobic --> blue
    'A':'#00BFFF',
    'V':'#87CEEB',
    'I':'#1E90FF',
    'L':'#6495ED',
    'M':'#4682B4',
    'F':'#4169E1',
    'Y':'#0000FF',
    'W':'#0000CD'
    }


"""
PLOTTING
"""
def plotly_plotting(df, label=True):

    fig = px.bar(df, x=df.index, y=df.columns, 
                 title="AA distribution around first Cystein",
                 labels={'position':'Position', 'value':'Percentage'},
                 color_discrete_map=colors_aa_by_group
                 )
    if label:
        fig.add_annotation(text="positively charged --> orange",
                          xref="paper", yref="paper",
                          x=0.8, y=0.975, showarrow=False, xanchor="left", yanchor="top")
        fig.add_annotation(text="negatively charged --> yellow",
                          xref="paper", yref="paper",
                          x=0.8, y=0.95, showarrow=False, xanchor="left", yanchor="top")
        fig.add_annotation(text="polar, uncharged --> red",
                          xref="paper", yref="paper",
                          x=0.8, y=0.925, showarrow=False, xanchor="left", yanchor="top")
        fig.add_annotation(text="special cases --> green",
                          xref="paper", yref="paper",
                          x=0.8, y=0.9, showarrow=False, xanchor="left", yanchor="top")
        fig.add_annotation(text="hydrophobic --> blue",
                          xref="paper", yref="paper",
                          x=0.8, y=0.875, showarrow=False, xanchor="left", yanchor="top")
    fig.show()

    


window_size = 50
path_to_yeast_proteome = "C:/Users/user/Desktop/BPC_FP/programming/python/reference.pkl"
significant_percentage = 0.01

def main(window_size, path_to_yeast_proteome, significant_percentage):
    
    proteome, aa_window = import_data(path_to_yeast_proteome, window_size)
    normalized_aa_window = normalize_dict(aa_window)
    plotly_df = get_significant_aa(proteome, significant_percentage, normalized_aa_window)
    plotly_plotting(plotly_df)
    
main(window_size, path_to_yeast_proteome, significant_percentage)
