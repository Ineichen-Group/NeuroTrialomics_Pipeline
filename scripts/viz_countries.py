import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import contextily as ctx
import geopandas as gpd
from matplotlib.colors import Normalize
from matplotlib.colorbar import ColorbarBase
from matplotlib.cm import ScalarMappable

def viz_countries_world_map(trial_metadata, output_file):
    unique_pairs_facility_countries = trial_metadata[['nct_id', 'country']].drop_duplicates()
    unique_pairs_facility_countries['country'] = unique_pairs_facility_countries['country'].replace({'United States': 'United States of America'})
    unique_pairs_facility_countries['country'] = unique_pairs_facility_countries['country'].replace({'Russian Federation': 'Russia'})
    unique_pairs_facility_countries['country'] = unique_pairs_facility_countries['country'].replace({'Korea, Republic of': 'South Korea'})

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    df_filtered = unique_pairs_facility_countries[unique_pairs_facility_countries['country'] != 'not reported']

    country_counts = df_filtered['country'].value_counts().reset_index()
    country_counts.columns = ['country', 'count']

    # Load world map
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

    # Merge the world GeoDataFrame with the country counts DataFrame
    world = world.merge(country_counts, how="left", left_on="name", right_on="country")

    # Apply a logarithmic transformation to the 'Frequency' column to deal with wide ranges in data
    world['log_count'] = np.log1p(world['count'])

    # Ensure that areas with zero (log1p(0) = 0) are left white by treating them as NaN
    world['log_count'] = world['log_count'].replace(0, np.nan)

    # Plot the data
    fig, ax = plt.subplots(1, 1, figsize=(15, 10))
    base = world.plot(ax=ax, column='log_count', cmap='Blues', legend=False, legend_kwds={'label': "Number of Clinical Trials by Country", 'orientation': "horizontal"})
    #ctx.add_basemap(ax, crs=world.crs.to_string(), source=ctx.providers.Stamen.Terrain)
    ax.set_axis_off()
    plt.title('World Map with of Clinical Trial Frequency', fontsize=15)

    # Create a custom colorbar
    norm = Normalize(vmin=world['log_count'].min(), vmax=world['log_count'].max())
    sm = ScalarMappable(cmap='Blues', norm=norm)
    sm._A = []  # Fake up the array of the scalar mappable.
    cb = plt.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.046, pad=0.04)
    cb.set_label('Number of Clinical Trials by Country', fontsize=14)

    # Format the ticks to show the actual counts
    tick_locs = np.linspace(world['log_count'].min(), world['log_count'].max(), num=5)
    cb.set_ticks(tick_locs)
    cb.set_ticklabels((np.exp(tick_locs) - 1).round().astype(int))  # Convert log count back to count

    # Save the figure as a PDF
    plt.savefig(output_file)

def main(metadata_file, output_file):
    # Load the input files
    trial_metadata = pd.read_csv(metadata_file)
    viz_countries_world_map(trial_metadata, output_file)

if __name__ == "__main__":
    if "snakemake" in globals():
        metadata_file = snakemake.input[0]
        output_file = snakemake.output[0]
        main(metadata_file, output_file)
    else:
        raise RuntimeError("This script should be run using Snakemake.")
