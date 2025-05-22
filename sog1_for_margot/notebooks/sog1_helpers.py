import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.ndimage import gaussian_filter1d
import warnings
warnings.filterwarnings("ignore")
# Set plot font and PDF font type
plt.rcParams["font.family"] = "Helvetica"
plt.rcParams['pdf.fonttype'] = 42

### Choose which activity column to use
activity_col = "Activity_S3_1"

def return_activities(description, activity_file_path="../output/SK_recalc_scores.csv", use_regex=False):
    """
        Extracts and merges activity data with a sequence library based on a description.
        Args:
            description (str): A string or regex pattern used to filter rows in the sequence library based on the "Description" column.
            activity_file_path (str, optional): Path to the CSV file containing activity data. Defaults to "../output/SK_recalc_scores.csv".
            use_regex (bool, optional): If True, treats the description as a regex pattern. Defaults to False.

        Returns:
            pandas.DataFrame: A merged DataFrame containing filtered sequence library rows and corresponding activity data.

        Raises:
            KeyError: If required columns ("Start", "End", "ProteinSeq", or "Error") are missing in the input data files.
    """
    # Load sequence library and filter rows by description
    library = pd.read_csv("../data/visit2_seq_lib.csv", index_col=0)
    if use_regex:
        library_rows = library[library["Description"].str.contains(description, regex=True)]
    else:
        library_rows = library[library["Description"].str.contains(description)]
    library_rows["Start"] = library_rows["Start"].astype(int)
    library_rows["End"] = library_rows["End"].astype(int)

    # Process "ProteinSeq" column to create "tile" column
    library_rows["tile"] = library_rows["ProteinSeq"].str.strip().str.upper()
    library_rows = library_rows.drop(columns="ProteinSeq")

    # Load activity data and process "tile" and "Error" columns
    activities = pd.read_csv(activity_file_path)
    activities = activities.rename(columns={"ProteinSeq": "tile"})
    activities["tile"] = activities["tile"].str.strip().str.upper()
    
    if "Error" in activities.columns:
        activities["activ_err_start"] = activities["Activity_S3_1"] - activities["Error"]
        activities["activ_err_end"] = activities["Activity_S3_1"] + activities["Error"]
    else:
        raise KeyError("Missing 'Error' column in activities.")

    # Merge library rows with activity data
    return pd.merge(library_rows, activities, on="tile", how="left")

def find_difference_indices(str1, str2, adjust=0):
    """
        Finds the indices where two strings differ.

        Parameters:
            str1 (str): The first string to compare.
            str2 (str): The second string to compare.
            adjust (int, optional): A value to adjust the indices by. Defaults to 0.

        Returns:
            list or None: A list of indices where the two strings differ, adjusted by the 
            `adjust` parameter. If the strings are of different lengths, the index of the 
            first unmatched character is included. Returns None if the strings are identical 
            or if the inputs are not strings.
    """
    if not isinstance(str1, str) or not isinstance(str2, str):
        return None
    diff_indices = [i + adjust for i in range(min(len(str1), len(str2))) if str1[i] != str2[i]]
    if len(str1) != len(str2):
        diff_indices.append(min(len(str1), len(str2)))
    return diff_indices if diff_indices else None

def significant(a, b):
    """
    Determines if two intervals do not overlap.

    This function checks whether two intervals, represented as tuples or lists
    of two numeric values (start, end), do not overlap. It returns True if the
    intervals do not overlap and False otherwise.

    Parameters:
        a (tuple or list): The first interval, defined as (start, end).
        b (tuple or list): The second interval, defined as (start, end).

    Returns:
        bool: True if the intervals do not overlap, False otherwise.

    Example:
        >>> significant((1, 5), (6, 10))
        True
        >>> significant((1, 5), (4, 8))
        False
    """
    return not bool(max(0, min(a[1], b[1]) - max(a[0], b[0])))

def add_var_positions(var_df, ref_df, activity_col, add_AAs=False):
    """
    Adds variant positions, activity differences, and optional amino acid differences 
    to a DataFrame by merging variant and reference data.

    Args:
        var_df (pd.DataFrame): DataFrame containing variant data with columns 
            "Start", "mid", "End", and "tile_var".
        ref_df (pd.DataFrame): DataFrame containing reference data with columns 
            "Start", "mid", "End", and "tile_wt".
        activity_col (str): Column name representing activity values in both 
            variant and reference DataFrames.
        add_AAs (bool, optional): If True, adds amino acid differences to the 
            resulting DataFrame. Defaults to False.

    Returns:
        pd.DataFrame: A merged DataFrame with the following additional columns:
            - "var": Variant positions relative to the reference.
            - "activ_diff": Difference in activity values between variant and reference.
            - "activ_fold_change": Fold change in activity values between variant and reference.
            - "signif_diff": Boolean indicating statistical significance of activity differences.
            - "AAs" (optional): Amino acid differences, if `add_AAs` is True.

    Notes:
        - The function assumes that the "Start", "mid", and "End" columns in both 
          `var_df` and `ref_df` are integers.
        - The `significant` function is used to determine statistical significance 
          based on error ranges.
        - The `find_difference_indices` function is used to identify sequence 
          differences between variant and reference tiles.    """
    
    # Ensure "Start", "mid", and "End" columns are integers
    for col in ["Start", "mid", "End"]:
        var_df[col] = var_df[col].astype(int)
        ref_df[col] = ref_df[col].astype(int)

    # Merge variant and reference DataFrames
    merged = pd.merge(var_df, ref_df, on=["Start", "mid", "End"], how="left", suffixes=("_var", "_wt"))

    # Find sequence differences
    merged["var"] = [find_difference_indices(row['tile_var'], row['tile_wt']) for _, row in merged.iterrows()]
    merged["var"] = merged["var"].apply(lambda x: x if x is not None else 0) + merged["Start"]

    # Calculate activity differences and fold changes
    merged["activ_diff"] = merged[activity_col + "_var"] - merged[activity_col + "_wt"]
    merged["activ_fold_change"] = merged[activity_col + "_var"] / merged[activity_col + "_wt"]

    # Check for statistical significance
    merged["signif_diff"] = merged.apply(
        lambda row: significant(
            (row["activ_err_start_var"], row["activ_err_end_var"]),
            (row["activ_err_start_wt"], row["activ_err_end_wt"])
        ), axis=1
    )

    # Optionally add amino acid differences
    if add_AAs:
        merged["AAs"] = [
            row["tile_wt"][var - row["Start"]] if var is not None and 0 <= var - row["Start"] < len(row["tile_wt"]) else None
            for _, row in merged.iterrows()
        ]
    return merged

def add_all_var_positions(var_df, ref_df, activity_col, add_AAs=False):
    """
    Identifies sequence differences between variant and reference dataframes 
    and optionally adds amino acid information.

    Args:
        var_df (pd.DataFrame): DataFrame containing variant sequences with 
            columns ["Start", "mid", "End", "tile_var"].
        ref_df (pd.DataFrame): DataFrame containing reference sequences with 
            columns ["Start", "mid", "End", "tile_wt"].
        activity_col (str): Column name representing activity data (not used 
            directly in the function but may be relevant for downstream processing).
        add_AAs (bool, optional): If True, adds a column with amino acid 
            differences at the identified variant positions. Defaults to False.

    Returns:
        pd.DataFrame: A merged DataFrame with the following additional columns:
            - "vars": List of indices where sequence differences occur.
            - "AAs" (optional): List of amino acids at the variant positions 
              (only if `add_AAs` is True).
    """
    merged = pd.merge(var_df, ref_df, on=["Start", "mid", "End"], how="left", suffixes=("_var", "_wt"))
    merged["vars"] = [
        find_difference_indices(row['tile_var'], row['tile_wt'], row['Start'])
        for _, row in merged.iterrows()
    ]
    if add_AAs:
        merged["AAs"] = [
            [row["tile_wt"][var - row["Start"]] for var in row["vars"]]
            for _, row in merged.iterrows()
        ]
    return merged

def plot_tile(start, end, activity, ax, color, center=True, label=None, alpha=1):
    """
    Plots a horizontal line (tile) on the given axis and optionally centers a scatter point on it.

    Parameters:
        start (float): The starting x-coordinate of the horizontal line.
        end (float): The ending x-coordinate of the horizontal line.
        activity (float): The y-coordinate where the horizontal line is plotted.
        ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the tile.
        color (str or tuple): The color of the horizontal line and scatter point.
        center (bool, optional): If True, a scatter point is plotted at the center of the line. Defaults to True.
        label (str, optional): The label for the scatter point, used for legend. Defaults to None.
        alpha (float, optional): The transparency level of the scatter point. Defaults to 1.

    Returns:
        None
    """
    ax.hlines(y=activity, xmin=start, xmax=end, color=color, lw=1, alpha=0.5, zorder=0)
    if center:
        sns.scatterplot(x=[(start + end) / 2], y=[activity], color=color, ax=ax, s=15, zorder=1, label=label, edgecolor='none', alpha=alpha)

def plot_all_tiles(merged_df, y_col, ax, color='red', center=True, label=None, alpha=1):
    """
    Parameters:
    -----------
    merged_df : pandas.DataFrame
        A DataFrame containing the data to be plotted. It must include columns 
        "Start", "End", and the column specified by `y_col`.
    y_col : str
        The name of the column in `merged_df` to be used for the y-axis values.
    ax : matplotlib.axes.Axes
        The matplotlib axis on which to plot the tiles.
    color : str, optional
        The color of the tiles. Default is 'red'.
    center : bool, optional
        If True, centers the tiles on the x-axis. Default is True.
    label : str or None, optional
        The label for the tiles. Only the first tile will use this label, 
        subsequent tiles will not have a label. Default is None.
    alpha : float, optional
        The transparency level of the tiles. Default is 1 (opaque).

    Returns:
    --------
    None
        This function modifies the given axis in place and does not return anything.
    """
    for _, row in merged_df.iterrows():
        plot_tile(row["Start"], row["End"], row[y_col], ax, color=color, center=center, label=label, alpha=alpha)
        label = None

def create_combo_df(start, df_c, df_s):
    """
    Creates a DataFrame containing combinations of variables and their associated activities.

    This function processes input DataFrames to generate a new DataFrame that includes 
    combinations of variables, their activities, and associated errors. It also calculates 
    the count of variables in each combination and ensures that the wild-type (WT) is included.

    Args:
        start (str or int): The starting point or identifier used to filter the DataFrames.
        df_c (pd.DataFrame): A DataFrame containing variable combinations and their activities. 
            It must include the following columns:
            - "vars": A column containing variable combinations.
            - "Start": A column used for filtering based on the `start` parameter.
            - Columns for activity and error values (e.g., `activity_col + "_var"`, `activity_col + "_wt"`, "Error_var", "Error_wt").
        df_s (pd.DataFrame): A DataFrame containing individual variables and their activities. 
            It must include the following columns:
            - "var": A column containing individual variable names.
            - "Start": A column used for filtering based on the `start` parameter.
            - Columns for activity and error values (e.g., `activity_col + "_var"`, "Error_var").

    Returns:
        pd.DataFrame: A DataFrame with the following columns:
            - "combo": The variable combinations (including "WT").
            - `activity_col`: The activity values for each combination.
            - "Error": The error values for each combination.
            - "count": The number of variables in each combination (0 for "WT").
            - "combo_str": String representation of the variable combinations.

    Notes:
        - The function assumes that `activity_col` is a global variable defined elsewhere in the code.
        - The input DataFrames must have the required columns; otherwise, the function may raise errors.
        - The function sorts the resulting DataFrame by the "count" column in ascending order.
    """
    df_c["vars_str"] = df_c["vars"].astype(str)
    df_c["var_count"] = df_c["vars_str"].str.count(",") + 1
    single_tile_df = df_c[df_c["Start"] == start]
    all_row = single_tile_df[single_tile_df["var_count"] == 3]
    tile_vars = all_row["vars"].iloc[0]

    var_combo = list(single_tile_df["vars"])
    activities = list(single_tile_df[activity_col + "_var"])
    errors = list(single_tile_df["Error_var"])

    for var in tile_vars:
        row = df_s[(df_s["var"] == var) & (df_s["Start"] == start)]
        var_combo.append(row["var"].iloc[0])
        activities.append(row[activity_col + "_var"].iloc[0])
        errors.append(row["Error_var"].iloc[0])

    var_combo.append("WT")
    activities.append(single_tile_df[activity_col + "_wt"].iloc[0])
    errors.append(single_tile_df["Error_wt"].iloc[0])

    activities_df = pd.DataFrame({
        "combo": var_combo,
        activity_col: activities,
        "Error": errors
    })
    activities_df["count"] = activities_df["combo"].astype(str).str.count(",") + 1
    activities_df["combo_str"] = activities_df["combo"].astype(str)
    activities_df.loc[activities_df['combo'] == "WT", 'count'] = 0
    return activities_df.sort_values(by="count").reset_index(drop=True)

def plot_labels(start, palette_dict, var_str, activities_df, ax, suffix="", x=-25):
    """
        Parameters:
        -----------
        start : int
            The starting index for slicing the amino acid features.
        palette_dict : dict
            A dictionary mapping counts to colors for labeling.
        var_str : str
            The string to use for variable positions in the label.
        activities_df : pandas.DataFrame
            A DataFrame containing activity data, including mutation combinations and counts.
        ax : matplotlib.axes.Axes
            The matplotlib Axes object where the labels will be added.
        suffix : str, optional
            A suffix to append to column names in `activities_df` (default is "").
        x : int, optional
            The x-coordinate for placing the labels (default is -25).

        Notes:
        ------
        - The function reads amino acid features from a CSV file located at "../data/Sog1_AA_features.csv".
        - Labels are constructed based on the presence of mutations in the `activities_df` DataFrame.
        - Colors for the labels are determined using the `palette_dict` based on mutation counts.
        - Labels are added to the plot using the `ax.text` method.
    """
    sog1_aa_features = pd.read_csv("../data/Sog1_AA_features.csv")
    for j, row in activities_df.iterrows():
        if j == 0:
            label = "".join(sog1_aa_features.iloc[start:start + 40]["aa"])
            color = palette_dict[0]
        else:
            label = ""
            tile_aas = list(sog1_aa_features.iloc[start:start + 40]["aa"])
            mut_aas = row["combo" + suffix]
            if isinstance(mut_aas, np.int64):
                mut_aas = np.array([mut_aas])
            for i in range(len(tile_aas)):
                label += var_str if i + start + 1 in mut_aas else "-"
            color = palette_dict[row["count" + suffix]]
        ax.text(x, j, label, ha="right", va="center", font="monospace", color=color, fontsize="x-small")

def plot_combination_activities(start, ax, df_c, df_s, var_str="A"): 
    """
        Parameters:
        -----------
        start : int
            The starting index or value used for creating the combination dataframe.
        ax : matplotlib.axes.Axes
            The matplotlib Axes object where the plot will be drawn.
        df_c : pandas.DataFrame
            The first dataframe containing data for combination activities.
        df_s : pandas.DataFrame
            The second dataframe containing supplementary data for combination activities.
        var_str : str, optional
            A string representing the variable name to be used in plot labels. Default is "A".

        Returns:
        --------
        None
            The function modifies the provided Axes object `ax` to include the plot.

        Notes:
        ------
        - The function uses a seaborn barplot to visualize activity levels.
        - A color palette ('rocket_r') is used to differentiate activity levels.
        - A vertical dotted line is added to indicate a reference activity level.
        - The y-axis labels are hidden for a cleaner visualization.
    """
    activities_df = create_combo_df(start, df_c, df_s)
    palette = sns.color_palette('rocket_r', n_colors=4)
    palette_dict = {i: palette[i] for i in range(4)}
    sns.barplot(data=activities_df, x=activity_col, y="combo_str", hue="count", palette=palette_dict, ax=ax, edgecolor="none", zorder=2)
    plot_labels(start, palette_dict, var_str, activities_df, ax)
    ax.axvline(activities_df[activity_col].iloc[0], color=palette_dict[0], linestyle="dotted", alpha=0.5, zorder=1, lw=3)
    ax.tick_params(axis='y', labelleft=False)

def return_avg_activities(name, pos_regex):
    """
        Computes smoothed average activity values for a given name and position regex.

    Args:
        name (str): Identifier for the dataset or entity.
        pos_regex (str): Regular expression to filter positions.

    Returns:
        pandas.DataFrame: A DataFrame containing:
            - 'range': The range of positions.
            - 'Activity_S3_1': Average activity values for each position.
            - 'smoothed_avg_activity': Smoothed activity values using a Gaussian filter.
    """
    BasicArTh = return_activities(name, pos_regex=pos_regex)
    BasicArTh["range"] = [np.arange(s, e + 1) for s, e in zip(BasicArTh["Start"], BasicArTh["End"])]
    BasicArTh = BasicArTh.explode("range").reset_index(drop=True)
    avg_activities = BasicArTh[["range", "Activity_S3_1"]].groupby(["range"]).mean().reset_index()
    avg_activities["smoothed_avg_activity"] = gaussian_filter1d(avg_activities["Activity_S3_1"], sigma=3)
    return avg_activities
