# -*- coding: utf-8 -*-
"""
@author: sekark
"""

import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
import copy

EXCEL_FILE = "OverFourHoursLC-QQQ.xlsx"
DATA_DICT = {"SC": {}, "SCD": {}}

COMPOUNDS_OF_INTEREST = ['Aspartate', 'Glutamate', 'GlycerolP', 'R5P', 'S7P']
COMPOUNDS_OF_INTEREST+=["Acetyl-CoA","Adenine","ADP-pentose","cAMP","cGMP","Citrate_isocitrate"]
COMPOUNDS_OF_INTEREST+=["FAD","FBP"	,"G6P"	,"Glutamine","Glycerol-P","GMP","Guanine"]
COMPOUNDS_OF_INTEREST+=["Ms1P","Panthothenate","PEP","Phenylalanine","Tryptophane","Tyrosine","UDP-hexose"]

COMPOUNDS_OF_INTEREST+= ['ATP', 'AMP', 'ADP']

TIMES_OF_INTEREST = {"SCD": [0,4], "SC": [0, 1, 2, 3, 4]}

def load_data():
    ec_df = pd.read_excel(EXCEL_FILE)
    
    compound_dict={}
    
    for compound in COMPOUNDS_OF_INTEREST:
        compound_dict[compound] = []
    
    compound_dict['OD'] = []
    
    for index,row in ec_df.iterrows():
        media = []
        if "SC" in row.Media and row.Media != 'SCD':
            media.append('SC')
        if "SCD" in row.Media:
            media.append('SCD')
        for medium in media:
            if row.Strain not in DATA_DICT[medium].keys():
                DATA_DICT[medium][row.Strain] = {}
            if row.Time not in DATA_DICT[medium][row.Strain].keys():
                DATA_DICT[medium][row.Strain][row.Time] = copy.deepcopy(compound_dict)
            
            for compound in COMPOUNDS_OF_INTEREST:
                DATA_DICT[medium][row.Strain][row.Time][compound].append(row[compound])
                
            DATA_DICT[medium][row.Strain][row.Time]['OD'].append(row.OD)


    
def analyze(mets = COMPOUNDS_OF_INTEREST):
    
    medias =  ["SC", "SCD"]
    strains =  ["WT", "dcbp2"]
    #mets = ['AMP','ATP', 'ADP']
    fig, axs = plt.subplots(len(mets),2, figsize=(6,3*len(mets)))
    
    for met in mets:
        ymins=[]
        ymaxs=[]
        for strain in strains:
    
            ax = axs[(mets.index(met)),(strains.index(strain))]
            
            for media in medias:
            
              
                times=TIMES_OF_INTEREST[media]
                
                OD_values = [np.mean(DATA_DICT[media][strain][time]['OD']) for time in times]
                current_values = [np.mean(DATA_DICT[media][strain][time][met]) for time in times]
                current_stderr = [np.std(DATA_DICT[media][strain][time][met])/np.sqrt(len(DATA_DICT[media][strain][time][met])) for time in times]
                            
                ax.errorbar(times,np.divide(current_values, OD_values), yerr = np.divide(current_stderr, OD_values), label=media)
                (ymin,ymax)=ax.get_ylim()
                ymins.append(ymin)
                ymaxs.append(ymax)
                
            ax.set_title(strain + " " + met)
            ax.set_ylabel("Measured Value")
        ax1=axs[(mets.index(met)),0]
        ax2=axs[(mets.index(met)),1]
        ax1.set_ylim(bottom=min(ymins),top=max(ymaxs))
        ax2.set_ylim(bottom=min(ymins),top=max(ymaxs))
    fig.savefig("intracellularplotsVert.svg")

    
if __name__ == "__main__":
    load_data()
    analyze()
    
