#This file contains code to draw bar plot of pLDDT distributions versus superkingdoms.
#Local file path: /home/lingyuncoding23/AF2_fasta/redo_PDB/extracted_AF2_model_organisms
import seaborn as sns
import matplotlib.pyplot as plt

# Set the DPI and figure size
plt.figure(figsize=(14, 8),dpi=700)

# Clean up the 'superkindom' column by stripping whitespace and non-breaking spaces
AF2_merged_shortened_superkindom_again['superkindom'] = AF2_merged_shortened_superkindom_again['superkindom'].str.strip()

# Define the color palette as a dictionary
custom_palettes = {
    "Protist": "#FEE127",
    "Viridiplantae": "#0eb519",
    "Fungi": "#F68A21",
    "Animalia": "#E41E26",  
    "Bacteria": "#3D5BA9",
    "Archaea": "#787979"
}
custom_order = ["Archaea","Bacteria",  "Fungi", "Animalia", "Viridiplantae","Protist"  ]
flierprops = dict(marker='o', markerfacecolor='black', markersize=0.3, linestyle='none')
# Plot a horizontal box plot
ax=sns.boxplot(data=AF2_merged_shortened_superkindom_again, x="plddt", y="superkindom", 
            palette=custom_palettes, orient='h',linewidth=4,width=0.7,order=custom_order,whis=(5, 95),
           flierprops=flierprops)
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
plt.yticks([], []) 
# Add labels and a title
plt.xlabel("")
plt.ylabel("")
# plt.title("Box Plot of plddt by Superkindom")

# Customize the plot further if needed
plt.xlim(20, 100)  # Match the x-axis limits with your data range
# Set custom x-axis and y-axis ticks
plt.xticks([20, 40, 60, 80, 100],[])
plt.tick_params(axis='both', width=6,length=12)
plt.legend([],[], frameon=False)
# sns.despine()
sns.despine(left=True, right=False, top=True)

ax.spines['right'].set_linewidth(6)  # Set right spine thickness
ax.spines['top'].set_linewidth(6)    # Set top spine thickness
ax.spines['bottom'].set_linewidth(6) # Set bottom spine thickness
ax.spines['left'].set_linewidth(6)     # Make sure left spine is invisible
plt.savefig('new_plddt_boxplot_no_legend_yaxisswitch_setframethickness_14_8_thick4_111424.png', dpi=700)

# Show the plot
plt.show()
