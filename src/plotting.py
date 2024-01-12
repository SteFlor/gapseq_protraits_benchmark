import seaborn as sns
from matplotlib import pyplot as plt
import os


def create_plot(df1, df2_converted, plot_values, performance_metrics, plot_name, outdir, flux_thr):
    """
    Creates and saves a plot based on the processed data.

    This function plots the comparison of true positives/negatives and false positives/negatives.
    It saves the plot in the specified directory with a detailed filename including metrics.

    """

    sns.set(style="whitegrid")
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 16))

    real_positives_plot,true_positives_plot,real_negatives_plot,false_positives_plot = plot_values
    precision,recall,f1_score,accuracy,specificity = performance_metrics


    # Plot for 1s
    ax[0].barh(df1.index, real_positives_plot, color='red', label='False negatives')
    ax[0].barh(df1.index, true_positives_plot, color='green', label='True positives')
    ax[0].set_title('Having phenotype')
    ax[0].legend()

    # Plot for 0s
    ax[1].barh(df1.index, real_negatives_plot, color='green', label='True negatives')
    ax[1].barh(df1.index, false_positives_plot, color='red', label='False positives')
    ax[1].set_title('Not having phenotype')
    ax[1].legend()

    plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(outdir + "/" + plot_name + "_thr" + str(flux_thr) + "_p" + str(precision) +
                "_r" + str(recall) + "_f1" + str(f1_score) + "_a" + str(accuracy) + "_s" + str(specificity) + ".png")

