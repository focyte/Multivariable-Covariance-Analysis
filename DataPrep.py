import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
from functools import reduce

# Load the RNAseq data
rna_data = pd.read_csv("CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt", sep='\t')

# Load and rename columns for proteomics dataframes
proteins_files = {
    "PHD2": "EGLN1 (Q9GZT9) Relative Protein Expression Proteomics.csv",
    "HIF1A": "HIF1A (Q16665) Relative Protein Expression Proteomics.csv",
    "LIMD1": "LIMD1 (Q9UGP4) Relative Protein Expression Proteomics.csv"
    }

protein_dataframes = {}
for protein, file in proteins_files.items():
    df = pd.read_csv(file)
    df.rename(columns={df.columns[1]: protein}, inplace=True)
    protein_dataframes[protein] = df

# Merge all protein dataframes on shared columns, drop rows with one or more protein missing
merge_columns = ["Depmap ID", "Primary Disease", "Cell Line Name", "Lineage"]
combined_df = reduce(lambda left, right: pd.merge(left, right, on=merge_columns, how="outer"), protein_dataframes.values())
combined_df = combined_df.dropna(subset=['HIF1A'])
print(combined_df)

# Save a tsv file of the Protein Values for modelling
x = combined_df[["PHD2", "HIF1A", "LIMD1"]]

x.to_csv("x.tsv", sep='\t', index=False)
print(x)

# Save combined dataframe to CSV
combined_df.to_csv("merged_df.csv", index=False)

# Select cell lines from combined_df to filter RNAseq data before merging
cell_lines = set(combined_df["Cell Line Name"])
cell_lines.add("gene_id")

# Filter RNAseq data to include only selected columns
columns_to_select = [col for col in rna_data.columns if any(keyword in col for keyword in cell_lines)]
filtered_data = rna_data[columns_to_select]
filtered_data.columns = filtered_data.columns.str.split('_').str[0]
filtered_data['gene'] = filtered_data['gene'].str.split('.').str[0]

# Load mart data for gene annotations and merge
mart = pd.read_csv("mart.csv")
filtered_data = pd.merge(filtered_data, mart, on='gene').drop(columns=['gene'])

# Reshape the filtered data
filtered_data_melted = filtered_data.melt(id_vars=['Gene name'], var_name='Cell Line Name')
filtered_data_pivoted = filtered_data_melted.pivot_table(index='Cell Line Name', columns='Gene name', values='value', aggfunc='mean').reset_index()
print(filtered_data_pivoted)

# Save filtered data to CSV
filtered_data_pivoted.to_csv('filtered.csv', index=False)

# Merge proteomics and RNAseq data
final_df = pd.merge(combined_df, filtered_data_pivoted, on='Cell Line Name')
final_df.to_csv("final.csv", index=False)

'''
To generate a hypoxia score using Buffa et al. gene list 
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2816644/
Genes were filtered by the Buffa list, normalised using sklearn and mean calculated
Top 25% Buffa score samples were selected as Case and the rest as REST
'''

Buffa = ["VEGFA", "SLC2A1", "PGAM1", "ENO1", "LDHA", "TPI1", "P4HA1", "MRPS17",
        "CDKN3", "ADM", "NDRG1", "TUBB6", "ALDOA", "MIF", "ACOT7", "MCTS1",
        "PSRC1", "PSMA7", "ANLN", "TUBA1B", "MAD2L2", "GPI", "TUBA1C", "MAP7D1",
        "DDIT4", "BNIP3", "C20orf20", "GAPDH", "MRPL13", "CHCHD2", "YKT6",
        "CORO1C", "SEC61G", "ANKRD37", "ESRP1", "PFKP", "SHCBP1", "CTSL2",
        "KIF20A", "SLC25A32", "UTP11L", "SLC16A1", "MRPL15", "KIF4A", "LRRC42",
        "PGK1", "HK2", "AK3L1", "CA9", "HILPDA", "PNP", "MRGBP", "CTSV", "UTP11", "AK4"]

missing_genes = [gene for gene in Buffa if gene not in final_df.columns]
present_genes = [gene for gene in Buffa if gene in final_df.columns]

with open("Buffa-Score-notes.txt", "w") as file:
    if missing_genes:
        file.write("The following genes from the Buffa list are not in the data and will be ignored:\n")
        for gene in missing_genes:
            file.write(f"{gene}\n")

# Filter the final_df to just Buffa hypoxia genes and calculate the score
columns_to_select = ["Cell Line Name"] + present_genes

final_df_buffa = final_df[columns_to_select].copy()
final_df_buffa[present_genes] = final_df_buffa[present_genes].fillna(0)
final_df_buffa[present_genes] = final_df_buffa[present_genes] + 1
final_df_buffa[present_genes] = final_df_buffa[present_genes].apply(np.log2, axis=1)
final_df_buffa["BuffaScore"] = final_df_buffa[present_genes].apply(np.median, axis=1)

# Use the Buffa Score to classify the cell lines into CASE (high hypoxia) and CNTR
quantile_75 = final_df_buffa["BuffaScore"].quantile(q=0.75)
final_df_buffa['Hypoxia'] = np.where(final_df_buffa["BuffaScore"] >= quantile_75, 'CASE', 'CNTR')
print(final_df_buffa)
final_df_buffa.to_csv("Buffa.csv", index=False)

# Plot a violin plot of the Buffa Score across cell lines
plt.figure(figsize=(10, 6))
sns.violinplot(data=final_df_buffa, x="BuffaScore")
plt.title("Hypoxia Score")
plt.ylabel("Buffa Score")
plt.savefig("BuffaScore_violin_plot.pdf")
plt.show()

# Save a tsv file of the BuffaScore column for modelling
y = final_df_buffa["BuffaScore"]
y.to_csv("y.tsv", sep='\t', index=False)

