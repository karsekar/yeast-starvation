# -*- coding: utf-8 -*-
"""
@author: sekark
"""

import pandas as pd
import numpy as np
from IPython.display import display
import matplotlib.pyplot as plt
import json
import copy
from scipy import stats

#define the file inputs
            
DATASET_EM = "inputs/matches_em.csv"
DATASET_PEAKLIST = "inputs/Yeast_starvation_norm.csv"

#these are jsons that were generated earlier for easier traversal
DATASET_DICTFILE = "inputs/dataset_dict.json"
DATASET_DICTFILE_REV = "inputs/dataset_dict_rev.json"

#define all the lipid classes per LIPIDMAPS
LIPID_CLASS = "inputs/lipidclass.txt"

FULL_LIPID_NAMES = ['Fatty Acyls [FA]', 'Glycerolipids [GL]', 
                    'Glycerophospholipids [GP]', 'Sphingolipids [SP]',
                    'Sterol Lipids [ST]', 'Prenol Lipids [PR]',
                    'Saccharolipids [SL]', 'Polyketides [PK]'
                    ]

DATASET_METADATA = {
            "SC": {
                    0: ["X170808_sample2_t0_A","X170808_sample2_t0_B","X170808_sample4_t0_A","X170808_sample4_t0_B"],
                    10: ["X170808_sample2_t1_A","X170808_sample2_t1_B","X170808_sample4_t1_A","X170808_sample4_t1_B"],
                    20: ["X170808_sample2_t2_A", "X170808_sample2_t2_B","X170808_sample4_t2_A","X170808_sample4_t2_B"],
                    30: ["X170808_sample2_t3_A","X170808_sample2_t3_B","X170808_sample4_t3_A","X170808_sample4_t3_B"]               
                    
                    },
            "SCD": {
                    0: ["X170808_sample1_t0_A","X170808_sample1_t0_B","X170808_sample3_t0_A", "X170808_sample3_t0_B"],
                    10: ["X170808_sample1_t1_A","X170808_sample1_t1_B","X170808_sample3_t1_A","X170808_sample3_t1_B"],
                    20: ["X170808_sample1_t2_A","X170808_sample1_t2_B","X170808_sample3_t2_A","X170808_sample3_t2_B"],
                    30: ["X170808_sample1_t3_A","X170808_sample1_t3_B","X170808_sample3_t3_A","X170808_sample3_t3_B"]
                    },
            "dpot1_SCD": {
                    0: ["X171018_sample2_t0_r1","X171018_sample2_t0_r2"],
                    10: ["X171018_sample2_t1_r1","X171018_sample2_t1_r2"],
                    20: ["X171018_sample2_t2_r1","X171018_sample2_t2_r2"],
                    30: ["X171018_sample2_t3_r1","X171018_sample2_t3_r2"]
                    },
            "dpot1_SC": {
                    0: ["X171018_sample3_t0_r1","X171018_sample3_t0_r2"],
                    10: ["X171018_sample3_t1_r1","X171018_sample3_t1_r2"],
                    20: ["X171018_sample3_t2_r1","X171018_sample3_t2_r2"],
                    30: ["X171018_sample3_t3_r1","X171018_sample3_t3_r2"]
                    },
            "Dpex6_SCD": {
                    0: ["X171018_sample4_t0_r1","X171018_sample4_t0_r2"],
                    10: ["X171018_sample4_t1_r1","X171018_sample4_t1_r2"],
                    20: ["X171018_sample4_t2_r1","X171018_sample4_t2_r2"],
                    30: ["X171018_sample4_t3_r1","X171018_sample4_t3_r2"]
                    },
            "Dpex6_SC": {
                    0: ["X171018_sample5_t0_r1","X171018_sample5_t0_r2"],
                    10: ["X171018_sample5_t1_r1","X171018_sample5_t1_r2"],
                    20: ["X171018_sample5_t2_r1","X171018_sample5_t2_r2"],
                    30: ["X171018_sample5_t3_r1","X171018_sample5_t3_r2"]
                    }
                                
        }
        

HIGH_LEVEL_CATEGORIES = ['FA', 'GL', 'GP', 'SP', 'ST', 'PR', 'SL', 'PK']

#function to load all of the data into the global variables
def load_data():
    
    global DATASET_EM_DF
    DATASET_EM_DF = pd.read_csv(DATASET_EM)
    
    global DATASET_PL_DF
    DATASET_PL_DF = pd.read_csv(DATASET_PEAKLIST)
    
    #generate a dictionary for the lipid classes
    global CLASS_DICT
    CLASS_DICT= {}
    
    f = open(LIPID_CLASS, 'r')

    for line in f:
        try:
            key = line.split('[')[1][:-2]
            CLASS_DICT[key] = line.split('[')[0].strip()
        except IndexError:
            continue
    
    global CLASS_LIST
    CLASS_LIST = list(CLASS_DICT.keys())
    
    global SUBCLASS_LIST
    SUBCLASS_LIST = [CLASS_LIST[i:i+8] for i in range(0,len(CLASS_LIST),8)]




#function to compare specific lipid classes between different conditions
def generate_plots_for_ds_compare(interested_list = []):
    with open(DATASET_DICTFILE_REV, 'r') as infile:
        ds_rev_dict = json.load(infile)
    
    shouldSave = False 
    
    if not len(interested_list):
        list_of_keys = list(ds_rev_dict.keys())
        list_of_keys.sort()
    else:
        list_of_keys = interested_list
        shouldSave = True
    
    media_list = ["SC", "dpot1_SC", "Dpex6_SC"] 
    media_list = ["SC", "SCD"]
        
    #for each page, we want 8 lipid classes at most
    groups_by_eight = [list_of_keys[i:i+8] for i in range(0,len(list_of_keys),8)]
    
    for group in groups_by_eight:
        plt.figure(figsize=(20,10))
        for category in group:
            plot_index = group.index(category)
            plt.subplot(int(str(24)+str(plot_index+1)))
            indices_of_int = ds_rev_dict[category]
            for media in media_list:
                
                timepoints = [time for time in DATASET_METADATA[media].keys()]
                condition = DATASET_METADATA[media]
                all_for_this_media = []
                for timepoint in timepoints:
                    all_for_this_media += condition[timepoint] 
                
                total_for_this_media = np.mean(np.sum(DATASET_PL_DF[all_for_this_media]))
                
                
                    
                values = [np.mean(np.mean(DATASET_PL_DF[DATASET_METADATA[media][time]].iloc[indices_of_int,:]))/total_for_this_media for time in timepoints]
                errors = [np.std(np.mean(DATASET_PL_DF[DATASET_METADATA[media][time]].iloc[indices_of_int,:]))/total_for_this_media for time in timepoints]
                plt.errorbar(timepoints,values, yerr = errors, label=media)
            plt.legend()
            plt.title(category)
        
        if shouldSave:
            plt.savefig('custom' + group[0] + 'to' + group[-1] + '.svg')
        plt.show()

   
def view_ds_data(interested_list = HIGH_LEVEL_CATEGORIES):
        
    #load the dictionary
    with open(DATASET_DICTFILE_REV, 'r') as infile:
        ds_rev_dict = json.load(infile)
    
    ds_rev_keys = list(ds_rev_dict.keys())
    ds_rev_keys.sort()
       
    
    media_list = ["SC", "SCD"] 
 
    media_dict={}
    for medium in media_list:
        media_dict[medium]=[]
        
    total_dict = {}
    
    #need a dictionary to compare across conditions
    compare_dict={}
    for item in interested_list:
        compare_dict[item]=copy.deepcopy(media_dict)
    
    
    
    plt.figure() 
    for media in media_list:   
        bar_dict = {}
        condition = DATASET_METADATA[media]
        timepoints = [time for time in condition.keys()]
        timepoints.sort()
        
        all_for_this_media = []
        for timepoint in timepoints:
            all_for_this_media += condition[timepoint] 
            
        
        total_for_this_media = np.sum(DATASET_PL_DF[all_for_this_media])
        
        for item in interested_list:
            current_entry =  np.sum(DATASET_PL_DF[all_for_this_media].iloc[ds_rev_dict[item],:])
            compare_dict[item][media].append(current_entry)
            bar_dict[item] = current_entry
        
        bar_y = [np.sum(bar_dict[item]) / np.sum(total_for_this_media) for item in interested_list]
        total_dict[media]=bar_dict[item]/np.sum(total_for_this_media)
        error_y = [np.std(bar_dict[item]) / np.sum(total_for_this_media) for item in interested_list]
        
        xoffset = media_list.index(media)
        
        
        xcoordinates = np.array(range(0,len(bar_y)))-0.25+xoffset*0.5
          
        plt.bar(xcoordinates, bar_y, width = 0.5, label=media)
        plt.errorbar(xcoordinates, bar_y, yerr = error_y, fmt='.', marker='None', ecolor = 'black')
    plt.xticks(range(0,len(interested_list)),FULL_LIPID_NAMES, rotation=90)
    
    for category in interested_list:
        a = np.array(compare_dict[category][media_list[0]])[0]
        b = np.array(compare_dict[category][media_list[1]])[0]
        t2, p2 = stats.ttest_ind(a,b)
    

    plt.legend()
    plt.title('Major groups')
    plt.savefig('bar_summary.svg')
    plt.show()
   

        
if __name__ == "__main__":
    load_data()    
    view_ds_data()    
    generate_plots_for_ds_compare(interested_list = ['PK01', 'PR0307'])
