#This file contains code for KDE of pLDDT distributions versus superkingdoms
#Local file path: /home/lingyuncoding23/AF2_fasta/redo_PDB/extracted_AF2_model_organisms
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 8),dpi=700)
AF2_merged_shortened_superkindom_again['superkindom'] = AF2_merged_shortened_superkindom_again['superkindom'].str.strip()
custom_palettes = {
    "Protist": "#FEE127",
    "Viridiplantae": "#0eb519",
     "Fungi": "#F68A21",
    "Animalia": "#E41E26",  
    "Bacteria": "#3D5BA9",
    "Archaea": "#787979"
}
sns.kdeplot(data=AF2_merged_shortened_superkindom_again[AF2_merged_shortened_superkindom_again['superkindom'] != 'Viridiplantae'],
            x="plddt", hue="superkindom", palette=custom_palettes, common_norm=False, common_grid=False, linewidth=4, alpha=1)
sns.kdeplot(data=AF2_merged_shortened_superkindom_again[AF2_merged_shortened_superkindom_again['superkindom'] == 'Viridiplantae'],
            x="plddt", hue="superkindom", palette=custom_palettes, common_norm=False, common_grid=False, linewidth=4, alpha=1.0)
# Add labels and a title
plt.xlabel('')
plt.ylabel('')
# plt.title("Normalized KDE Plot of plddt by Superkindom")
plt.xlim(0, 110)
plt.ylim(0, 0.1)
# Set custom x-axis and y-axis ticks
plt.xticks([0,20, 40, 60, 80, 100],[])
plt.yticks([0,0.02, 0.04, 0.06, 0.08, 0.1],[])
plt.tick_params(axis='both', width=3,length=6)
sns.despine()
ax = plt.gca()  # Get the current axes
ax.spines['left'].set_linewidth(3)   # Left spine thickness
ax.spines['right'].set_linewidth(3)  # Right spine thickness
ax.spines['top'].set_linewidth(3)    # Top spine thickness
ax.spines['bottom'].set_linewidth(3) # Bottom spine thickness
plt.legend([],[], frameon=False)
# Remove x-axis and y-axis numbers (ticks)
plt.savefig('new_plddt_density_plot_no_legend_no_tranparency_setframethickness_0_110tick4_111424.png', dpi=700)
# Show the plot
plt.show()
