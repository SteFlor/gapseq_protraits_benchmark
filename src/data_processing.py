
import pandas as pd

def binarize_results(df1, df2, flux_thr):
    df1 = df1.copy()
    df1 = df1.sort_index(ascending=False)
    df2 = df2.copy()
    # Convert df2 values to 1 where nonzero, and 0 otherwise
    df2_binarized = df2.map(lambda x: '1' if abs(x) >= flux_thr else '0')

    # Ensure df2_binarized has all the metabolites from df1
    for metabolite in df1.index:
        if metabolite not in df2_binarized.index:
            df2_binarized.loc[metabolite] = '0'  # Add metabolite with all 0s if not present

    # Align df2_binarized to have the same order of metabolites as df1
    df2_binarized = df2_binarized.reindex(df1.index, fill_value='0')

    real_positives_plot = df1.T.apply(pd.Series.value_counts).loc['1'].fillna(0)
    real_negatives_plot = df1.T.apply(pd.Series.value_counts).loc['0'].fillna(0)

    true_positives_plot = (df1 == '1') & (df2_binarized == '1')
    false_positives_plot = (df1 == '0') & (df2_binarized == '1')
    true_positives_plot = true_positives_plot.sum(axis=1)
    false_positives_plot = false_positives_plot.sum(axis=1)

    #Collect values for plots in a tuple
    plot_values = (real_positives_plot,true_positives_plot,real_negatives_plot,false_positives_plot)

    real_positives = real_positives_plot.sum()
    true_positives = true_positives_plot.sum()
    real_negatives = real_negatives_plot.sum()
    false_positives = false_positives_plot.sum()
    
    false_negatives = real_positives - true_positives
    true_negatives = real_negatives - false_positives
    
    precision = true_positives / (true_positives+false_positives)
    recall = true_positives / (true_positives+false_negatives)
    f1_score = 2*(precision*recall)/(precision+recall)
    accuracy = (true_positives + true_negatives) / (real_positives + real_negatives)
    specificity = true_negatives / (true_negatives+false_positives)

    # Rounding metrics and collect them in a tuple
    precision = precision.round(2)
    recall = recall.round(2)
    f1_score = f1_score.round(2)
    accuracy = accuracy.round(2)
    specificity = specificity.round(2)

    performance_metrics = (precision,recall,f1_score,accuracy,specificity)


    return df1, df2_binarized, plot_values, performance_metrics
