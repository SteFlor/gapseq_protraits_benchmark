import pandas as pd
import numpy as np
import os

def clean_duplicated_uptake_exchangesFVA(df):
    df = df.copy()
    # List of all column names
    all_columns = df.columns

    for col in all_columns:
        # Check if it's a 'columnName Exchange' without a corresponding 'columnName'
        if (' Exchange_uptake' in col and col.replace(' Exchange_uptake', '_uptake') not in all_columns) or (' Exchange_production' in col and col.replace(' Exchange_production', '_production') not in all_columns):
            # Rename the standalone 'Exchange' column
            df.rename(columns={col: col.replace(' Exchange_uptake', '_uptake')}, inplace=True)
            df.rename(columns={col: col.replace(' Exchange_production', '_production')}, inplace=True)
        elif ' Exchange_uptake' not in col and ' Exchange_production' not in col:
            if not col.endswith("_uptake") and not col.endswith("_production"):
                raise ValueError("Colname " + col + " must either end in '_uptake' or '_production'")
            if col.endswith("_uptake"):
                exchange_col = col.replace("_uptake"," Exchange_uptake")
            if col.endswith("_production"):
                exchange_col = col.replace("_production"," Exchange_production")
            # if both "normal" 'metabolite_name_{uptake_or_production}' and 'metabolite_name_Exchange_{uptake_or_production}' versions of the reaction exists, only one will be kept
            if exchange_col in all_columns:
                # Keep the value that has the maximum value (I take the max (and not the sum) because the results are from FVA)
                df[col] = np.maximum(df[col],df[exchange_col])
                # Remove the Exchange "version" of the column
                df = df.drop(exchange_col, axis=1)  

    return df


def read_and_merge_dataframes(directory, file_identifier):
    # List all files in the directory
    all_files = os.listdir(directory)

    # Filter files that match the identifier
    filtered_files = [file for file in all_files if file_identifier in file]

    # Read and merge dataframes
    merged_df = pd.DataFrame()
    for file in filtered_files:
        df = pd.read_csv(os.path.join(directory, file),index_col=0)
        merged_df = pd.merge(merged_df, df, left_index=True, right_index=True, how="outer")

    return merged_df
