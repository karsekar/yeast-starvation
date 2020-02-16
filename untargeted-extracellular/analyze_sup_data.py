# -*- coding: utf-8 -*-
"""
@author: sekark
"""

import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt

IONS_FILE = 'data.xlsx'
ANNOTATION_FILE = 'annotation.xls'
LABELS_FILE = 'labels.xlsx' 

COMPOUND_NAMES = {'C00049': 'Aspartate', 'C00073': 'Methionine'}

COMPOUNDS_OF_INTEREST = list(COMPOUND_NAMES.keys())

DATA_DICT = {"SC": {}, "SCD": {}}

TP_REMOVE = {"SC": [], "SCD": [2,18]}


def load_data():
        ions_df = pd.read_excel(IONS_FILE)
        
        
        ann_df = pd.read_excel(ANNOTATION_FILE)
        
        
        labels_df = pd.read_excel(LABELS_FILE)
        
        compound_dict = {}
        
        compound_to_ion  = {}
        
        for compound in COMPOUNDS_OF_INTEREST:
            compound_dict[compound] = []
            for index,row in ann_df[ann_df.id == compound].iterrows():
                if '+' not in row['mod']:
                    compound_to_ion[compound] = np.array(row['ion'])
            
        
        
        for index,row in labels_df.iterrows():
            media = []
            if "SC" in row.Media and row.Media != 'SCD':
                media.append('SC')
            if "SCD" in row.Media:
                media.append('SCD')
            
            for medium in media:
                
                if row.Strain not in DATA_DICT[medium].keys():
                    DATA_DICT[medium][row.Strain] = {}
                if row.Timepoint not in DATA_DICT[medium][row.Strain].keys():
                    DATA_DICT[medium][row.Strain][row.Timepoint] = copy.deepcopy(compound_dict)
            
            
                for compound in compound_to_ion.keys():
                    temp_value = np.array(ions_df[int(compound_to_ion[compound])])[index]
            
                    if medium != 'na' and row.Timepoint != 'na' and row.Strain != 'na':
                        DATA_DICT[medium][row.Strain][row.Timepoint][compound].append(temp_value)
               
        

        

def analyze():
    
    medias =  ["SC", "SCD"]
    strains =  ["WT", "cbp"]
    
    for strain in strains:
         fig, axs = plt.subplots(1, len(COMPOUNDS_OF_INTEREST), figsize=(7, 3.5))
    
         for met in COMPOUNDS_OF_INTEREST:
             
             ax = axs[COMPOUNDS_OF_INTEREST.index(met)]
             for media in medias:
                 times = [int(key) for key in DATA_DICT[media][strain].keys()]
                 times.sort()
            
                 
                 current_values = [np.mean(DATA_DICT[media][strain][time][met]) for time in times]
                 current_stderr = [np.std(DATA_DICT[media][strain][time][met])/np.sqrt(len(DATA_DICT[media][strain][time][met])) for time in times]
        
                 ax.errorbar(times,current_values, yerr = current_stderr, label=media)
             ax.set_title(strain + " " + COMPOUND_NAMES[met])
             ax.set_ylabel("Measured Value")
             ax.legend()
         plt.savefig(strain + "_ec.svg")
        
if __name__ == "__main__":
    load_data()
    analyze()