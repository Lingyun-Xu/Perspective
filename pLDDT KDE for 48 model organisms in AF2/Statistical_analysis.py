#This file contains code for more statistics analysis of pLDDT ditributions.
#Local file path: /home/lingyuncoding23/AF2_fasta/redo_PDB/extracted_AF2_model_organisms
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import scikit_posthocs as sp
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde

# Ensure the sequence lengths have been calculated
AF2_merged_shortened_superkindom_again['sequence_length'] = AF2_merged_shortened_superkindom_again['AAseq_one_letter_code'].str.len()

# Get the unique superkingdom groups
superkingdom_groups = AF2_merged_shortened_superkindom_again['superkindom'].unique()

# Initialize a dictionary to store results
results = []

# Loop through each superkingdom and compute correlations
for group in superkingdom_groups:
    # Subset data for the current superkingdom
    group_data = AF2_merged_shortened_superkindom_again[AF2_merged_shortened_superkindom_again['superkindom'] == group]
    
    # Compute Pearson correlation
    pearson_corr, pearson_p = stats.pearsonr(group_data['sequence_length'], group_data['plddt'])
    
    # Compute Spearman correlation
    spearman_corr, spearman_p = stats.spearmanr(group_data['sequence_length'], group_data['plddt'])
    
    # Store the results
    results.append({
        'superkingdom': group,
        'pearson_corr': pearson_corr,
        'pearson_p': pearson_p,
        'spearman_corr': spearman_corr,
        'spearman_p': spearman_p
    })

# Convert the results to a DataFrame for easy viewing
results_df = pd.DataFrame(results)

print(results_df)
'''    results_df output:
     superkingdom  pearson_corr     pearson_p  spearman_corr    spearman_p
0        Bacteria      0.188714  0.000000e+00       0.236299  0.000000e+00
1        Animalia     -0.186181  0.000000e+00      -0.113990  0.000000e+00
2        Protist      -0.243323  0.000000e+00      -0.181773  0.000000e+00
3           Fungi     -0.055514  1.225737e-53      -0.051769  7.213875e-47
4  Viridiplantae       0.147177  0.000000e+00       0.242971  0.000000e+00
5         Archaea      0.050844  3.229297e-02       0.142592  1.632023e-09  '''
mean_plddt_per_superkingdom = AF2_merged_shortened_superkindom_again.groupby('superkindom')['plddt'].mean()
''' mean_plddt_per_superkingdom output:
Mean plddt for each superkingdom:
superkindom
Animalia          75.155447
Archaea           89.796136
Bacteria          87.383879
Fungi             76.043622
Protist           71.441903
Viridiplantae     72.888478
Name: plddt, dtype: float64  '''

'''The Kruskal-Wallis test is robust for large datasets and non-normal distributions.
The Dunn’s post-hoc test allows you to explore pairwise differences if the Kruskal-Wallis test is significant.'''
superkingdom_groups = AF2_merged_shortened_superkindom_again['superkindom'].unique()

# Visualize distributions with boxplots
plt.figure(figsize=(10, 6))
sns.boxplot(data=AF2_merged_shortened_superkindom_again, x='superkindom', y='plddt')
plt.title('Distribution of plddt across Superkingdoms')
plt.xticks(rotation=45)
plt.show()
# Perform Kruskal-Wallis test for non-parametric comparison of distributions
kruskal_stat, kruskal_p = stats.kruskal(
    *[AF2_merged_shortened_superkindom_again[AF2_merged_shortened_superkindom_again['superkindom'] == group]['plddt']
      for group in superkingdom_groups]
)
print(f"Kruskal-Wallis Test: p-value = {kruskal_p}")
# Perform pairwise comparisons using Dunn's test for post-hoc analysis
dunn_test = sp.posthoc_dunn(
    AF2_merged_shortened_superkindom_again, 
    val_col='plddt', 
    group_col='superkindom', 
    p_adjust='bonferroni'
)
print("Pairwise post-hoc Dunn's test results (p-values):")
print(dunn_test)
'''Kruskal-Wallis Test: p-value = 0.0
Pairwise post-hoc Dunn's test results (p-values):
                    Animalia       Archaea      Bacteria         Fungi  
Animalia        1.000000e+00  0.000000e+00  0.000000e+00  2.344512e-95   
Archaea         0.000000e+00  1.000000e+00  2.545244e-16  0.000000e+00   
Bacteria        0.000000e+00  2.545244e-16  1.000000e+00  0.000000e+00   
Fungi           2.344512e-95  0.000000e+00  0.000000e+00  1.000000e+00   
Protist         0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00   
Viridiplantae   0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00   

                    Protist   Viridiplantae   
Animalia        0.000000e+00    0.000000e+00  
Archaea         0.000000e+00    0.000000e+00  
Bacteria        0.000000e+00    0.000000e+00  
Fungi           0.000000e+00    0.000000e+00  
Protist         1.000000e+00    3.339619e-51  
Viridiplantae   3.339619e-51    1.000000e+00 ''' 

#Perform the same code for all superkingdoms.
Archaea_df = AF2_merged_shortened_superkindom_again[AF2_merged_shortened_superkindom_again['superkindom'].str.contains('Archaea', case=False)]
# Plot the KDE and get density values
kde_plot = sns.kdeplot(data=Archaea_df, x="plddt", hue="superkindom", common_norm=True, common_grid=True)
density_values = kde_plot.get_lines()[0].get_ydata()
# Find peaks in the KDE
peak_indices = find_peaks(density_values)[0]
x_values = kde_plot.get_lines()[0].get_xdata()
peak_x_coords = x_values[peak_indices] 
print("Peak locations in Archaea_df:", peak_x_coords)
# To find the range where 95% of the data is located
plddt_values = Archaea_df['plddt'].dropna()
# Calculate the 2.5th and 97.5th percentiles (for the range covering 95% of the data)
percentile_2_5 = np.percentile(plddt_values, 2.5)
percentile_97_5 = np.percentile(plddt_values, 97.5)
print("Archaea_df 95% of data is located between:", percentile_2_5, "and", percentile_97_5)
# Filter the dataframe to get rows where plddt >= 70
plddt_gte_70 = Archaea_df[Archaea_df['plddt'] >= 70]
# Calculate the percentage of rows with plddt >= 70
percentage_gte_70 = (len(plddt_gte_70) / len(Archaea_df)) * 100
print(f"Archaea_df entries number is:{len(Archaea_df)}")
print(f"Archaea_df Percentage of data with plddt >= 70: {percentage_gte_70:.2f}%")
''' Peak locations in Archaea_df: [38.99735839 53.97425381 69.3351722  94.68068753]
Archaea_df 95% of data is located between: 62.336 and 97.66
Archaea_df entries number is:1773
Archaea_df Percentage of data with plddt >= 70: 95.04%

Peak locations in Bacteria_df: [27.57445694 95.01133439]
Bacteria_df 95% of data is located between: 57.37 and 97.64
Bacteria_df entries number is:61727
Bacteria_df Percentage of data with plddt >= 70: 91.94%

Peak locations in Fungi_df: [82.86129937 92.93241194]
Fungi_df 95% of data is located between: 39.31 and 96.26
Fungi_df entries number is:77026
Fungi_df Percentage of data with plddt >= 70: 68.99%

Peak locations in Animalia_df: [87.29089468]
Animalia_df 95% of data is located between: 42.919500000000006 and 95.21
Animalia_df entries number is:204479
Animalia_df Percentage of data with plddt >= 70: 67.31%

Peak locations in Viridiplantae_df: [64.41501301 82.38740874]
Viridiplantae_df 95% of data is located between: 40.99 and 94.54
Viridiplantae_df entries number is:166181
Viridiplantae_df Percentage of data with plddt >= 70: 59.87%

Peak locations in Protist_df: [59.63981651 80.43959227]
Protist_df 95% of data is located between: 38.11475 and 95.10525000000001
Protist_df Percentage of data with plddt >= 70: 57.34%
