import pandas as pd
import numpy as np
from src.accessory import clean_duplicated_uptake_exchangesFVA, read_and_merge_dataframes
from src.plotting import create_plot
from src.data_processing import binarize_results
import os

# Directory where the input dataframes (gapseq's output) are stored
input_dfs = "FVA_1_fullmedium"
directory = "DataFramesFromCluster/" + input_dfs

# Merge dataframes
uptakes_df = read_and_merge_dataframes(directory, 'uptakes_df')
mediums_df = read_and_merge_dataframes(directory, 'mediums_df')
growth_products_df = read_and_merge_dataframes(directory, 'growth_products_df')

# Save the complete medium, uptakes and products dataframes
outdir_intermediate = "Intermediate/" + input_dfs
os.makedirs(outdir_intermediate, exist_ok=True)
uptakes_df.to_csv(outdir_intermediate +"/merged_uptakes_df.csv")
mediums_df.to_csv(outdir_intermediate +"/merged_mediums_df.csv")
growth_products_df.to_csv(outdir_intermediate +"/merged_growth_products_df.csv")

outdir = "Results/" + input_dfs
os.makedirs(outdir, exist_ok=True)

# Add suffix to distinguish between production and uptake
uptakes_renamed = uptakes_df.copy()
uptakes_renamed.index = [reac+"_uptake" for reac in uptakes_renamed.index]
growth_products_renamed = growth_products_df.copy()
growth_products_renamed.index = [reac+"_production" for reac in growth_products_renamed.index]

# Create summary of all uptakes and productions
# Uptake originally had negative value, an Na indicates that the model cannot handle the molecule (hence equivalent to 0). The minus sign makes all uptake values positive
summary = (-uptakes_renamed).transpose().merge(growth_products_renamed.transpose(), left_index=True, right_index=True)
summary = summary.fillna(0.0)

# Filter out zero columns
summary = summary.loc[:, summary.sum() != 0]

# Load the genome table
genome_table = pd.read_csv("genome_table_protraits", sep="\t", header=None, index_col=1)
genome_table.columns = ["Tax_ID"]
genome_table.index.name = "id"

# Merge genome table with the summary of uptakes and productions
summary = summary.merge(genome_table, how="inner", left_index=True, right_index=True)

# Load protraits and merge it with the summary, to get a species level summary
protraits = pd.read_csv("ProTraits_binaryIntegratedPr0.90.txt", sep="\t",index_col=1)
summary_named = summary.merge(protraits["Organism_name"], how="inner", left_index=False, left_on="Tax_ID", right_index=True)
summary_named = summary_named.set_index("Organism_name").drop("Tax_ID",axis=1)
cleaned = clean_duplicated_uptake_exchangesFVA(summary_named) # remove reaction with redundant naming ("Exchange...") while keeping the flux inormation

# Create a dataframe with the species level aggreagate (mean, median and max would all make sense)
species_level = cleaned.groupby('Organism_name').max() #TODO:replace all occurrences of .max() with an argument (e.g. via groupby(...).agg(["max"]) or smth similar)

# Load dictionary that translates from Gapseq to Protraits molecule names
gapseq_to_protraits = pd.read_csv("gapseq_to_protraits.csv", index_col=0).to_dict()["Protraits_id"]

protraits = protraits.set_index("Organism_name")
rows = list(species_level.index)
columns = list(set(gapseq_to_protraits.values())) #all the protraits metabolites that can be linked to gapseq metabs
protraits = protraits.loc[rows]
protraits = protraits[columns]

# Select all the uptakes from the species-level aggregated dataframe and translate metabolite names from Gapseq to Protraits names
all_uptakes = species_level.copy()
all_uptakes = all_uptakes[[el for el in all_uptakes.columns if el.endswith("_uptake")]]
all_uptakes.columns = [el.replace("_uptake","") for el in all_uptakes.columns]
#TODO: move the " Exchange" postfix handling in the "clean_dupli..." accessory function
all_uptakes.columns = [el.replace(" Exchange","") for el in all_uptakes.columns]
all_uptakes.columns = [gapseq_to_protraits[el] if el in gapseq_to_protraits.keys() else el for el in all_uptakes.columns]
all_uptakes = all_uptakes.T.groupby(all_uptakes.columns).max() # after translating the gapseq names to protraits, duplicate columns may appear (many-to-one relationship b/w gapseq met. names and protraits met.names)

#take all uptakes and append indole_production
all_phenos = all_uptakes.copy()
all_phenos.loc["indole_production"] = summary_named["indol_production"].groupby('Organism_name').max()

# List of metabolites that are tested in Jakob's pipeline
jakobs_list = ["d-arabinose","d-galacturonate","d-xylose","fructose","g-amino-butyric_acid","glutamate","glycerol","indole_production","lactate","lactose","melibiose","sucrose","trehalose","alanine","l-arabinose","l-rhamnose","acetate","maltose","aspartate","glycine","histidine","l-arginine","leucine","mannose","serine","galactose","l-fucose","ribose"]

flux_thr = 2.0 # arbitrary threshold, is related (i.e. must be proportional to) the flux-molecule weight scaling used in the simulation script. When the growth will be assessed on each carbon source separately, flux-molecule weight scaling will not be necessary as a simulation input (so this "dependency" between these 2 parameter will not be there

df1_processed, df2_binarized, plot_values, performance_metrics = binarize_results(protraits.T, all_phenos, flux_thr)
create_plot(df1_processed, df2_binarized, plot_values, performance_metrics, input_dfs, outdir, flux_thr)
#precision, recall, f1_score, accuracy, specificity = convert_and_plot(protraits.T, all_phenos, input_dfs, flux_thr=flux_thr, outdir="Results/"+input_dfs)
precision, recall, f1_score, accuracy, specificity = performance_metrics
performance_metrics_string = "Precision: " + str(precision) + "\n" + "Recall: " + str(recall) + "\n" + "F1 Score: " + str(f1_score) + "\n" + "Accuracy: " + str(accuracy) + "\n" + "Specificity: " + str(specificity)


df1_processed, df2_binarized, plot_values, performance_metrics = binarize_results(protraits.T.loc[jakobs_list], all_phenos, flux_thr)
create_plot(df1_processed, df2_binarized, plot_values, performance_metrics, input_dfs+"_Jakobs_list", outdir, flux_thr)
#precision, recall, f1_score, accuracy, specificity = convert_and_plot(protraits.T.loc[jakobs_list], all_phenos, input_dfs + "_jakobs_list", flux_thr=flux_thr, outdir="Results/"+input_dfs)
precision, recall, f1_score, accuracy, specificity = performance_metrics
performance_metrics_string_jakobs_list = "Precision: " + str(precision) + "\n" + "Recall: " + str(recall) + "\n" + "F1 Score: " + str(f1_score) + "\n" + "Accuracy: " + str(accuracy) + "\n" + "Specificity: " + str(specificity)

df2_binarized.T.to_csv(outdir + "/binarized_predictions.csv")	



print("\nResults including all modelled metabolites matching between Gapseq and Protraits\n") 
print(performance_metrics_string+"\n")
print("_____________________________________\n")
print("Results including all metabolites matching between Gapseq and Protraits and that are tested in Jakob's pipeline\n") 
print(performance_metrics_string_jakobs_list+"\n")
print("\nPlot and binarized predicitons have been saved under the 'Results' folder\n")



