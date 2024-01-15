import cobra
import seaborn as sns
import pandas as pd 
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse


# Argument parser setup
parser = argparse.ArgumentParser(description='Process batch number and size.')
parser.add_argument('batch_num', type=int, help='Batch number for processing')
parser.add_argument('batch_size', type=int, help='Size of each batch')
parser.add_argument('models_dir', type=str, help='Directory where the input models are')
parser.add_argument('frac_opt', type=float, help='Fraction of optimum for FVA simulation')
parser.add_argument('out_dir', type=str, help='Directory where the output dataframes will be stored')
args = parser.parse_args()


def update_medium(models, reac_map_path, inflow_value):

    reac_map = pd.read_csv(reac_map_path, index_col = 1)
    medium_toadd = dict()
    for protraits_id in set(reac_map.index):
        medium_toadd[protraits_id] = [vals[0] if len(vals)==1 else vals for vals in reac_map.loc[protraits_id].values]
    for model in models:
        for metabs in medium_toadd.values():
            print(metabs)
            for metab in metabs:
                current_metab_and_model_ex_reacs = [reac for reac in models[model].boundary if reac.name.replace("-e0 Exchange","").replace("-e0","") == metab]
                for reac in current_metab_and_model_ex_reacs:
                    reac.lower_bound = -inflow_value/list(reac.metabolites.keys())[0].formula_weight # a rough way to keep the "caloric input" of each uptake comparable among each other
                    print(reac.id)
                    print(reac.lower_bound)


def fill_mediums_dfs(models, exchange_reac_name_postfix="-e0 Exchange"): #takes a dict of models, returns a dataframe containing the media of all models
    mediums_df = dict()
    
    for model_key in models:
        
         # select model
        model = models[model_key]
        
        # save medium as pd.Series in a dictionary
        medium = {reac.name.replace(exchange_reac_name_postfix,""): -reac.lower_bound for reac in model.boundary if reac.lower_bound < 0}
        medium = pd.Series(medium)
        mediums_df[model_key] = medium

    #transform mediums_df dictionary into a dataframe
    mediums_df = pd.DataFrame(mediums_df)
    return mediums_df


def get_uptakes_and_products_FVA(models, exchange_reac_name_postfix="-e0 Exchange", normalize_by_growth=True, fraction_of_optimum = 0.01):
    boundaries_dict = dict()
    uptakes_df = dict()
    growth_products_df = dict()
    for model_key in models:
        print(model_key)
    
        # select model
        model = models[model_key]
        
        # save list of boundary reactions
        boundary_names = {reac.name.replace(exchange_reac_name_postfix,"") for reac in model.boundary}
        boundaries_dict[model_key] = boundary_names        
        
        # select boundary ids for later selection of growth products and uptakes
        boundary_ids = [reac.id for reac in model.boundary]
        
        # calculate FVA of the required reactions, with the required fraction of optimum
        fva_result = cobra.flux_analysis.flux_variability_analysis(model, reaction_list=boundary_ids, fraction_of_optimum=fraction_of_optimum)

        # remove reactions postfix to simplify reactions naming
        fva_result.index = [model.reactions.get_by_id(el).name.replace(exchange_reac_name_postfix,"") for el in fva_result.index]
                
        ## Save fva optimal (minimum) uptakes as pd.Series in a dictionary
        # Modify the 'minimum' column: set values >= -1e-06 to 0
        fva_result.loc[fva_result['minimum'] >= -1e-06, 'minimum'] = 0.0
        # Modify the 'maximum' column: set values <= 1e-06 to 0
        fva_result.loc[fva_result['maximum'] <= 1e-06, 'maximum'] = 0.0
        # Extract the 'minimum' column as a Series
        uptakes = fva_result['minimum']
        # Extract the 'maximum' column as a Series
        growth_products = fva_result['maximum']
        #in case of duplicated compound name (same compund but different exchange reaction, group by compound name (id) and pick the maximum absolute value (min for uptakes, max for excretion)
        uptakes_df[model_key] = uptakes.groupby(level=0).min() 
        growth_products_df[model_key] = growth_products.groupby(level=0).max()
        
    uptakes_df=pd.DataFrame(uptakes_df)
    growth_products_df = pd.DataFrame(growth_products_df)
        
    return uptakes_df, growth_products_df, boundaries_dict
            



batch_num = args.batch_num
batch_size = args.batch_size
models_dir = args.models_dir
frac_opt = args.frac_opt
frac_opt_perc = str(int(frac_opt*100))
last_index = len(os.listdir(models_dir + "/"))-1
start_index = batch_size*batch_num
end_index= min(start_index + batch_size, last_index)
outdir_prefix = args.out_dir
outdir = outdir_prefix + "/dataframes_FVA_" + frac_opt_perc + "_fullmedium"
os.makedirs(outdir, exist_ok=True)

models = dict()
model_files = [model_file for model_file in os.listdir(models_dir + "/") if model_file.endswith("genomic.xml")]
for model_file in model_files[start_index:end_index]:
    models[model_file.replace(".xml","")] = cobra.io.read_sbml_model(models_dir + "/" + model_file)

update_medium(models, "dat/gapseq_to_protraits.csv", 400) #arbitrary value that fits the threshold of 2.0 on resulting uptakes for detection of growth
uptakes_df, growth_products_df, boundaries_dict = get_uptakes_and_products_FVA(models, fraction_of_optimum = frac_opt)

uptakes_df.to_csv(outdir + "/uptakes_df_"+str(batch_num)+".csv")

growth_products_df.to_csv(outdir + "/growth_products_df_"+str(batch_num)+".csv")

mediums_df = fill_mediums_dfs(models)

mediums_df.to_csv(outdir + "/mediums_df_"+str(batch_num)+".csv")






