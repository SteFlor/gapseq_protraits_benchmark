import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
import os

def convert_and_plot(df1, df2, plot_name, flux_thr, outdir):
    df1 = df1.copy()
    df1 = df1.sort_index(ascending = False)
    df2 = df2.copy()
    # Convert df2 values to 1 where nonzero, and 0 otherwise
    df2_converted = df2.map(lambda x: '1' if abs(x) >= flux_thr else '0')

    # Ensure df2_converted has all the metabolites from df1
    for metabolite in df1.index:
        if metabolite not in df2_converted.index:
            df2_converted.loc[metabolite] = '0'  # Add metabolite with all 0s if not present

    # Align df2_converted to have the same order of metabolites as df1
    df2_converted = df2_converted.reindex(df1.index, fill_value='0')
    
    # Plotting
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 16))    

    real_positives = df1.T.apply(pd.Series.value_counts).loc['1'].fillna(0)
    real_negatives = df1.T.apply(pd.Series.value_counts).loc['0'].fillna(0)
    
    true_positives = (df1 == '1') & (df2_converted == '1')
    false_positives = (df1 == '0') & (df2_converted == '1')
    true_positives = true_positives.sum(axis=1)
    false_positives = false_positives.sum(axis=1)

    # Plot for 1s
    ax[0].barh(df1.index, real_positives, color='red', label='False negatives')
    ax[0].barh(df1.index, true_positives, color='green', label='True positives')
    ax[0].set_title('Having phenotype')
    ax[0].legend()

    # Plot for 0s
    ax[1].barh(df1.index, real_negatives, color='green', label='True negatives')
    ax[1].barh(df1.index, false_positives, color='red', label='False positives')
    ax[1].set_title('Not having phenotype')
    ax[1].legend()
       
        
    real_positives = real_positives.sum()
    true_positives = true_positives.sum()
    real_negatives = real_negatives.sum()
    false_positives = false_positives.sum()
    
    false_negatives = real_positives - true_positives
    true_negatives = real_negatives - false_positives
    
    precision = true_positives / (true_positives+false_positives)
    recall = true_positives / (true_positives+false_negatives)
    f1_score = 2*(precision*recall)/(precision+recall)
    accuracy = (true_positives + true_negatives) / (real_positives + real_negatives)
    specificity = true_negatives / (true_negatives+false_positives)
    
    precision = precision.round(2)
    recall = recall.round(2)
    f1_score = f1_score.round(2)
    accuracy = accuracy.round(2)
    specificity = specificity.round(2)
    
    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(outdir+"/"+plot_name+"_thr"+str(flux_thr)+"_p"+str(precision)+
                "_r"+str(recall)+"_f1"+str(f1_score)+"_a"+str(accuracy)+"_s"+str(specificity)+".png")

    return precision, recall, f1_score, accuracy, specificity

