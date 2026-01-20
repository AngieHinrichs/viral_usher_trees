"""
Use BTE to read in an UShER protobuf file and calculate root-to-tip distances.
Read in dates from metadata, plot root-to-tip distances vs. dates, and calculate R^2.
"""

import argparse
import bte
import gzip
import re
import sys

import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def rtt_calc(node: bte.MATNode, node_distance: int, distances: dict):
    """Recursively descend from node, accumulating distance from root, and store leaf distances"""
    node_distance += len(node.mutations)
    if len(node.children) == 0:
        # Leaf node, add to distances
        distances['name'].append(node.id)
        distances['distance'].append(node_distance)
    else:
        for child in node.children:
            rtt_calc(child, node_distance, distances)


def get_rtt(input_tree_file: str):
    """Read in the tree with BTE, traverse it to get root-to-tip distance for each leaf,
    and return a dataframe mapping leaf names to distances."""
    tree = bte.MATree(input_tree_file)
    distances = {'name': [], 'distance': []}
    rtt_calc(tree.root, 0, distances)
    return pd.DataFrame(distances)


def get_dates(metadata_file: str, key_column: str = None):
    """Read in names and dates from a TSV metadata file.  If key_column is given, find names there,
    else use the first column.  Ignore dates that don't parse.
    Return a data frame mapping names to pandas datetimes."""
    dates = {'name': [], 'date': []}
    with gzip.open(metadata_file, "rt") if metadata_file.endswith(".gz") else open(metadata_file) as f:
        line = f.readline().rstrip()
        columns = line.split("\t")
        if key_column is None:
            key_column = columns[0]
        elif key_column not in columns:
            print(f"Error: key_column {key_column} is not in columns in first line of metadata file {metadata_file} ({columns})",
                  file=sys.stderr)
            sys.exit(1)
        name_ix = columns.index(key_column)
        try:
            date_ix = columns.index('date')
        except ValueError as e:
            print(f"Error: metadata file {metadata_file} does not have a 'date' column", file=sys.stderr)
            sys.exit(1)
        for line in f:
            columns = [x.strip() for x in line.split("\t")]
            date = columns[date_ix]
            m = re.match(r"^[0-9]{4}(-[0-9]{2}(-[0-9]{2})?)?$", date)
            if m:
                got_month, got_day = m.groups()
                if not got_month:
                    date += "-07-01"
                elif not got_day:
                    date += "-15"
                dates['name'].append(columns[name_ix])
                dates['date'].append(pd.to_datetime(date))
    return pd.DataFrame(dates)


def plot_distance_vs_date(df_distance, df_date, title="Distance vs Date"):
    """Plot distance vs date with regression line and r-squared."""
    # Merge the dataframes
    df = pd.merge(df_distance, df_date, on='name', how='inner')
   
    # Remove rows where either distance or date is missing
    df = df.dropna(subset=['distance', 'date'])

    df = df.sort_values('date')
    
    # Convert to numeric for regression
    date_numeric = (df['date'] - df['date'].min()).dt.days
    
    # Calculate regression
    print("Running stats.linregress on merged dataframe...")
    slope, intercept, r_value, p_value, std_err = stats.linregress(date_numeric, df['distance'])
    regression_line = slope * date_numeric + intercept
    print(f"slope={slope:.6f}, R^2={r_value**2:.6f}, intercept={intercept:.3f}, std_err={std_err:.6f}")
    
    # Create plot
    print("Plotting results...")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(df['date'], df['distance'], alpha=0.6, s=50, label='Distance')
    ax.plot(df['date'], regression_line, 'r-', linewidth=2, 
            label=f'Regression line (R² = {r_value**2:.3f})')
    
    ax.set_xlabel('Date')
    ax.set_ylabel('Distance')
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.xticks(rotation=45)
    plt.tight_layout()
    
    return fig, ax


def main():
    parser = argparse.ArgumentParser(description="Plot root-to-tip distances vs. date for an UShER tree with metadata")
    parser.add_argument("-i", "--input_tree", required=True,
                        help="Input UShER protobuf tree file")
    parser.add_argument("-m", "--metadata", required=True,
                        help="TSV metadata file for tree")
    parser.add_argument("-k", "--key_column",
                        help="Column of metadata that contains names found in tree")
    parser.add_argument("-t", "--title",
                        help="Title for the plot")
    parser.add_argument("-o", "--output_base", required=True,
                        help="Output file base name for .pdf and .png plots")
    args = parser.parse_args()
    print(f"Calculating root-to-tip distances from {args.input_tree}...")
    rtt = get_rtt(args.input_tree)
    print(f"Reading dates from {args.metadata}...")
    dates = get_dates(args.metadata, args.key_column)
    print("Plotting distance vs date...")
    fig, ax = plot_distance_vs_date(rtt, dates, args.title)
    
    print("Saving output and displaying plot")
    fig.savefig(f'{args.output_base}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{args.output_base}.png', dpi=300, bbox_inches='tight')
    

if __name__ == "__main__":
    main()
