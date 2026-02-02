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
        except ValueError:
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


def merge_dataframes(df_distance, df_date):
    """Merge df_distance and df_date on name, drop items where distance or date is missing, sort by date."""
    # Merge the dataframes
    df = pd.merge(df_distance, df_date, on='name', how='inner')
    df = df.dropna(subset=['distance', 'date'])
    df = df.sort_values('date')
    return df


def plot_distance_vs_date(tree_data: list, overall_title: str):
    """Make a stack of three plots: distance vs date with regression line and r-squared.
    tree_data is [ {df:..., title:...}, ... ] but if df is None leave a blank space (no data for that tree)"""
    fig, axes = plt.subplots(3, 1, figsize=(10, 12))
    summary = []
    for idx, tree_info in enumerate(tree_data):
        df = tree_info['df']
        title = tree_info['title']
        ax = axes[idx]
        if df is None:
            ax.text(0.5, 0.5, 'No tree',
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue
        if len(df) < 2:
            ax.text(0.5, 0.5, 'Insufficient data for regression',
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue
        # Convert to numeric for regression
        date_numeric = (df['date'] - df['date'].min()).dt.days
        # Calculate regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(date_numeric, df['distance'])
        summary.append([title, f"{r_value**2:.6f}", f"{slope:.6f}", f"{intercept:.3f}", f"{std_err:.6f}"])
        regression_line = slope * date_numeric + intercept
        # Create plot
        ax.scatter(df['date'], df['distance'], alpha=0.6, s=50, label='Sample')
        ax.plot(df['date'], regression_line, 'r-', linewidth=2,
                label=f'Regression line (R² = {r_value**2:.3f})')
        ax.set_xlabel('Date')
        ax.set_ylabel('Distance')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.xticks(rotation=45)
    fig.suptitle(overall_title, fontsize=14, y=0.995)
    plt.tight_layout()
    return fig, axes, summary


def main():
    parser = argparse.ArgumentParser(description="Plot root-to-tip distances vs. date for up to three UShER trees with metadata")
    parser.add_argument("-m", "--metadata", required=True,
                        help="TSV metadata file for all three trees")
    parser.add_argument("-a", "--input_tree_a", required=True,
                        help="Input UShER protobuf tree file A (first of three)")
    parser.add_argument("-A", "--key_column_a",
                        help="Column of metadata that contains names found in tree A")
    parser.add_argument("--title_a",
                        help="Title for plot A (first of three)")
    parser.add_argument("-b", "--input_tree_b",
                        help="Input UShER protobuf tree file B (second of three)")
    parser.add_argument("-B", "--key_column_b",
                        help="Column of metadata that contains names found in tree B")
    parser.add_argument("--title_b",
                        help="Title for plot B (second of three)")
    parser.add_argument("-c", "--input_tree_c",
                        help="Input UShER protobuf tree file C (third of three)")
    parser.add_argument("-C", "--key_column_c",
                        help="Column of metadata that contains names found in tree C")
    parser.add_argument("--title_c",
                        help="Title for plot C (third of three)")
    parser.add_argument("-t", "--title",
                        help="Overall title for stacked plots")
    parser.add_argument("-o", "--output_base", required=True,
                        help="Output file base name for .pdf and .png plots and .tsv regression stats")
    args = parser.parse_args()
    rtt_a = get_rtt(args.input_tree_a)
    dates_a = get_dates(args.metadata, args.key_column_a)
    df_a = merge_dataframes(rtt_a, dates_a)
    if args.input_tree_b:
        rtt_b = get_rtt(args.input_tree_b)
        dates_b = get_dates(args.metadata, args.key_column_b)
        df_b = merge_dataframes(rtt_b, dates_b)
    else:
        df_b = None
    if args.input_tree_c:
        rtt_c = get_rtt(args.input_tree_c)
        dates_c = get_dates(args.metadata, args.key_column_c)
        df_c = merge_dataframes(rtt_c, dates_c)
    else:
        df_c = None
    tree_data = [{'df': df_a, 'title': args.title_a if args.title_a else "A"},
                 {'df': df_b, 'title': args.title_b if args.title_b else "B"},
                 {'df': df_c, 'title': args.title_c if args.title_c else "C"}]
    fig, axes, summary = plot_distance_vs_date(tree_data, args.title)

    fig.savefig(f'{args.output_base}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{args.output_base}.png', dpi=300, bbox_inches='tight')

    with open(f"{args.output_base}.tsv", 'w') as f:
        f.write("\t".join(["method", "r_squared", "slope", "intercept", "std_err"]) + "\n")
        for row in summary:
            f.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
