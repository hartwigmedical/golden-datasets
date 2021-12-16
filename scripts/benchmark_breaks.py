#!/usr/bin/env python3
# -*- coding: utf-8 -*-



"""
Script to obtain the metrics of the benchmarking not considering
the type of SV for the evaluation. 

Rationale: 
    Check presence of a test variant in the truth file using only chr, start, 
    and end. Use different windows of bp to define it. 
    
    If there are several occurences that can be used to match the variant, 
    choose the closest one, and not allow for its use for a later comparison. 
    

"""


import pandas as pd
import numpy as np
import os
import argparse


def breaks_comparison(test, truth, window):
    
    """
    TP: All test variants that are within windows in the truth file.
        -NO check for already used truth variants in the evaluation
    FP: All test variants that are not within windows in the truth file.
    FN: All truth variants that have not been used after checking all test variants
        -Need to store in variables the truth variants that have been used
        -Wondering if it would be the same other way around and implement windows
        in the test variants. I think so, so we will go this way.
        
    """
    
    used_truth_TP = []
    used_test_TP = []
    used_test_FP = []
    used_truth_FN = []
    TP = 0
    FP = 0
    FN = 0
    
    # for index, row in test.iterrows():
    #     truth_matches = ""
    #     truth_matches = truth.loc[(truth['start_chrom'] == row['start_chrom']) & (truth['end_chrom'] == row['end_chrom']) \
    #                           & (abs(truth['start'] - row['start']) <  window) & (abs(truth['end'] - row['end']) <  window)]
    #     if truth_matches.empty is False:
    #         #Meaning if there is a truth variant that has been matched:
    #         pass
    #         #Now len(used_truth_TP) == len(used_test_TP) should match     
            
    #         
    #     else:
    #         FP+=1
        
    #         used_test_FP.append(row.values.tolist())


    for index, row in truth.iterrows():
        test_matches = ""
        test_matches = test.loc[(test['start_chrom'] == row['start_chrom']) & (test['end_chrom'] == row['end_chrom']) \
                              & (abs(test['start'] - row['start']) <  window) & (abs(test['end'] - row['end']) <  window)]
        if test_matches.empty is False:
            TP+=1
            used_truth_TP.append(row.values.tolist())
            
            #For the cases where there are several truth matches
            
            if len(test_matches) > 1:

            #Grab the closest one by start and end site. 
                
                values = []
                for idx, variant in test_matches.iterrows():
                    
                    values.append(abs(variant['start'] - row['start']) + abs(variant['end'] - row['end']))
                    #TODO: Add caution message in case they have the same min value
                #This way of doing it, only selects one variant.   
                selected_test_df = test_matches.iloc[values.index(min(values))]
                selected_test = selected_test_df.values.tolist()

                     
            else:
                selected_test = test_matches.values.tolist()[0]

            used_test_TP.append(selected_test)
            
        else:
          FN+=1
          used_truth_FN.append(row.values.tolist())

    used_test_FP = [item for item in test.values.tolist() if item not in used_test_TP]
    return(used_truth_TP,used_test_FP,used_truth_FN, used_test_TP)


    
def compute_metrics(used_truth_TP,used_test_FP,used_truth_FN, window):
    
    metrics = [["Window", "TP", "FP", "FN", "Recall", "Precision", "F1-score" ]]
    
    TP_df = pd.DataFrame(used_truth_TP, columns=['start_chrom','start','end_chrom', 'end', 'ref', 'alt', 'length', 'type'])
    FP_df = pd.DataFrame(used_test_FP, columns=['start_chrom','start','end_chrom', 'end', 'ref', 'alt', 'length', 'type'])
    FN_df = pd.DataFrame(used_truth_FN, columns=['start_chrom','start','end_chrom', 'end', 'ref', 'alt', 'length', 'type'])
    
    
    if TP_df.empty:
        TP = 0
    else:
        TP = TP_df.shape[0]
    
    if FP_df.empty:
        FP = 0
    else:
        FP = FP_df.shape[0]

    if FN_df.empty:
        FN = 0
    else:
        FN = FN_df.shape[0]

    if (TP + FN)==0:
        recall = 0
    else:
        recall = TP / (TP + FN)
    
    if (TP + FP)==0:
        precision= 0        
    else:
        precision = TP / (TP + FP)
    
    if (recall + precision) == 0:
        F1 = 0
    else:
        F1 = 2 * (recall * precision) / (recall + precision)

    metrics = [["Window", "TP", "FP", "FN", "Recall", "Precision", "F1-score" ]]
    
    metrics.append([window, TP, FP, FN, round(recall,2), round(precision,2), round(F1,2)])
    
    return(metrics)

    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("df_node", help="Path to SV dataframe test file")
    parser.add_argument("df_truth", help="Path to SV dataframe truth file")
    parser.add_argument('-w', '--window_list', help='delimited list input', type=str)
    parser.add_argument("-metrics", help="Save metrics in a file")
    #parser.add_argument("-metrics", help="If dataframe should be saved to CSV, give output file path here")

    args = parser.parse_args()

    test = pd.read_csv(args.df_node, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': int, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    truth = pd.read_csv(args.df_truth, index_col=0, dtype={'start_chrom': str, 'start': int, 'end_chrom': str, 'end': int, 'ref': str, 'alt': str, 'length': 'Int32', 'type': str})
    windows_list = [int(item) for item in args.window_list.split(',')]

    #This loop will be run as many times as windows have been introduced by the user
    agregate_metrics = []
    for window in windows_list:
        used_truth_TP,used_test_FP,used_truth_FN,used_test_TP = breaks_comparison(test, truth, window)
        metrics = compute_metrics(used_truth_TP,used_test_FP,used_truth_FN, window)
        agregate_metrics.append(metrics[1])
    
    metrics_df = pd.DataFrame(agregate_metrics, columns=["Window", "TP", "FP", "FN", "Recall", "Precision", "F1-score" ])


    #Save to file
    if args.metrics:
        metrics_df.to_csv(args.metrics, index=False, sep='\t')
        
  
        
if __name__ == "__main__":
    main()

