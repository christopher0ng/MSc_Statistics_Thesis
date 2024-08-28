import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd

df = pd.read_csv("C:/Users/chris/OneDrive - Imperial College London/Desktop/Imperial Masters 2023/MSC Statistics/Summer Research Project/total_cases_quantiles.csv")

# Create a 2x2 subplot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Years to plot
years = [2020, 2021, 2022, 2023]

# Colors for each year
colors = ['blue', 'green', 'orange', 'red']

# Loop through each year and create the plot
for i, year in enumerate(years):
    ax = axes[i//2, i%2]  # Determine the subplot position
    start_idx = i  # Start index for this year
    mean_values = []
    lower_bound = []
    upper_bound = []
    
    # Extract data for all 8 variants for this specific year
    for variant in range(8):
        mean_values.append(df.iloc[start_idx + variant * 4, 1])
        lower_bound.append(df.iloc[start_idx + variant * 4, 0])
        upper_bound.append(df.iloc[start_idx + variant * 4, 2])
    
    variants = np.arange(1, 9)
    
    # Plot the mean line with a specific color
    ax.plot(variants, mean_values, marker='o', color=colors[i], label='Mean')
    
    # Plot the shaded area for the 95% confidence interval with the same color
    ax.fill_between(variants, lower_bound, upper_bound, color=colors[i], alpha=0.2, label='95% CI')
    
    # Set the title and labels
    ax.set_title(f'{year} Dengue Cases')
    ax.set_xlabel('Variant')
    ax.set_ylabel('Number of Cases')
    ax.set_xticks(variants)
    
    # Position the legend in the top-left corner
    ax.legend(loc="upper left")

# Adjust layout
plt.tight_layout()

# Save and show the plot
plt.savefig('variants_cases_by_year.png')
plt.show()

variant1 = pd.read_csv('Variant_1.csv')
v1_2020 = variant1.iloc[0::4].reset_index(drop=True)  # Rows 1, 5, 9, etc.
v1_2021 = variant1.iloc[1::4].reset_index(drop=True)  # Rows 2, 6, 10, etc.
v1_2022 = variant1.iloc[2::4].reset_index(drop=True)  # Rows 3, 7, 11, etc.
v1_2023 = variant1.iloc[3::4].reset_index(drop=True)  # Rows 4, 8, 12, etc.

gdf['2020_lower'] = v1_2020['2.5%']
gdf['2020_mean'] = v1_2020['50%']
gdf['2020_upper'] = v1_2020['97.5%']
gdf['2021_lower'] = v1_2021['2.5%']
gdf['2021_mean'] = v1_2021['50%']
gdf['2021_upper'] = v1_2021['97.5%']
gdf['2022_lower'] = v1_2022['2.5%']
gdf['2022_mean'] = v1_2022['50%']
gdf['2022_upper'] = v1_2022['97.5%']
gdf['2023_lower'] = v1_2023['2.5%']
gdf['2023_mean'] = v1_2023['50%']
gdf['2023_upper'] = v1_2023['97.5%']

#Variant 1 bins

bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2020_lower_bin'] = pd.cut(gdf['2020_lower'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2020_mean_bin'] = pd.cut(gdf['2020_mean'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2020_upper_bin'] = pd.cut(gdf['2020_upper'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2021_lower_bin'] = pd.cut(gdf['2021_lower'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2021_mean_bin'] = pd.cut(gdf['2021_mean'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2021_upper_bin'] = pd.cut(gdf['2021_upper'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2022_lower_bin'] = pd.cut(gdf['2022_lower'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2022_mean_bin'] = pd.cut(gdf['2022_mean'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2022_upper_bin'] = pd.cut(gdf['2022_upper'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2023_lower_bin'] = pd.cut(gdf['2023_lower'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2023_mean_bin'] = pd.cut(gdf['2023_mean'], bins=bins, labels=labels)
bins = [-float('inf'), 10, 100, 1000, float('inf')]
labels = ['less than 10','10-100', '100-1000','1000+']
gdf['2023_upper_bin'] = pd.cut(gdf['2023_upper'], bins=bins, labels=labels)

legends_list =['2020 2.5% quantile','2020 mean','2020 97.5% quantile','2021 2.5% quantile','2021 mean','2021 97.5% quantile','2022 2.5% quantile','2022 mean','2022 97.5% quantile','2023 2.5% quantile','2023 mean','2023 97.5% quantile']
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(15, 12))  # Adjust figsize as needed
axes = axes.flatten()  # Flatten the 2D array of axes to make iteration easier
columns = ['2020_lower_bin','2020_mean_bin','2020_upper_bin','2021_lower_bin','2021_mean_bin','2021_upper_bin','2022_lower_bin','2022_mean_bin','2022_upper_bin','2023_lower_bin','2023_mean_bin','2023_upper_bin']

for i, ax in enumerate(axes):
    gdf.plot(column=columns[i], ax=ax, legend=True, categorical=True,
                legend_kwds={'loc': "lower left",'fontsize': 'small'},
                cmap='OrRd')
    ax.set_title(legends_list[i])

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig('Variant 1 results')

# Show the plot
plt.show()

#Similar code for variants 2 to 8