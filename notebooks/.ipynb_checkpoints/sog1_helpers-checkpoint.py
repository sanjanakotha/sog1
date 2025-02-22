import pandas as pd
#from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#import warnings

# Ignore all warnings
#warnings.filterwarnings('ignore')

#plt.rcParams['text.usetex'] = True

plt.rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Helvetica" #somethings this one doesnt work
plt.rcParams['pdf.fonttype'] = 42


activity_col = "Activity_S3_1"

# Returns EC spike activities using the description column
# Extracts position using regex passed in
def return_activities(description, pos_regex = r'(\d+)_.*'):
    library = pd.read_csv("../data/visit2_seq_lib.csv", index_col = 0)    

    library_rows = library[library["Description"].str.contains(description)]
    library_rows["Start"] = library_rows["Description"].str.extract(pos_regex).astype(int)
    library_rows["Start"] = 10 * library_rows["Start"] - 9
    library_rows.loc[library_rows['Start'] ==411, 'Start'] = 410
    library_rows["mid"] = library_rows["Start"] + 20
    library_rows["End"] = library_rows["Start"] + 40
    library_rows["tile"] = library_rows["ProteinSeq"].astype(str).str.strip()
    library_rows = library_rows.drop(columns = {"ProteinSeq"})
    
    activities = pd.read_csv("../data/Sog1_library2_activities_with_reads.csv")
    #activities = pd.read_csv("../data/Sog1_library2_activities_with_reads_ECspike.csv")
    #activities = pd.read_csv("../data/Sog1_library2_activities_with_reads_EC.csv")
    activities = activities.rename(columns = {"AAseq" : "tile"})
    activities["tile"] = activities["tile"].astype(str).str.strip()
      
    return pd.merge(library_rows, activities[["tile", "Activity_S3_1", "Activity_S3_2", "lib2_avg"]], on = "tile", how = "left")

# Returns index of first difference between two strings
def find_difference_index(str1, str2):
    if type(str1) != str or type(str2) != str:
        return None
    # Iterate over the characters of the strings and find the index where they differ
    for i in range(min(len(str1), len(str2))):
        if str1[i] != str2[i]:
            return i
    # If the strings are identical until the shortest length, return the length of the shortest string
    if len(str1) != len(str2):
        return min(len(str1), len(str2))
    return None  # Return None if the strings are identical

# Returns all differing indices between two strings
def find_difference_indices(str1, str2, adjust = 0):
    if type(str1) != str or type(str2) != str:
        return None
    # Initialize a list to store indices where characters differ
    diff_indices = []
    # Iterate over the characters of the strings and find all indices where they differ
    for i in range(min(len(str1), len(str2))):
        if str1[i] != str2[i]:
            diff_indices.append(i + adjust)
    
    # If the strings are of different lengths, add the remaining indices
    if len(str1) != len(str2):
        diff_indices.append(min(len(str1), len(str2)))
        print("warning: Diff lens")
        
    return diff_indices if diff_indices else None  # Return None if no differences are found

# Adds index of first position that varies between var_df and corresponding sequence in ref_df
# var_df and ref_df must share Start, mid, end columns
def add_var_positions(var_df, ref_df, activity_col):
    merged = pd.merge(var_df, ref_df, on = ["Start", "mid", "End"], how = "left", suffixes = ("_var", "_wt"))
    diffs = []
    for i in merged.index:
        diffs.append(find_difference_index(merged['tile_var'].iloc[i], 
                                        merged['tile_wt'].iloc[i]))

    merged["var"] = diffs
    merged["var"] = merged["var"] + merged["Start"]
    merged["activ_diff"] = merged[activity_col + "_var"] - merged[activity_col + "_wt"]
    merged["activ_fold_change"] = merged[activity_col + "_var"] / merged[activity_col + "_wt"]
    return merged

# Adds all positions that vary between var_df and corresponding sequence in ref_df
# var_df and ref_df must share Start, mid, end columns
def add_all_var_positions(var_df, ref_df, activity_col):
    merged = pd.merge(var_df, ref_df, on = ["Start", "mid", "End"], how = "left", suffixes = ("_var", "_wt"))
    diffs = []
    for i in merged.index:
        diffs.append(find_difference_indices(merged['tile_var'].iloc[i], 
                                        merged['tile_wt'].iloc[i],
                                            merged['Start'].iloc[i]))

    merged["vars"] = diffs
    return merged

# Plots individual tile
def plot_tile(start, end, activity, ax, color, center = True):
    ax.hlines(y=activity, xmin=start, xmax=end, color=color, lw = 1, alpha = 0.5, zorder = 0)
    if center:
        sns.scatterplot(x = [(start + end) / 2], y = [activity], color = color, alpha = 1, ax = ax, s = 15, zorder = 1)

# Plots all tiles
def plot_all_tiles(merged_df, y_col, ax, color = 'red', center = True):
    for i in merged_df.index:
        row = merged_df.loc[i]
        plot_tile(row["Start"],
                  row["End"],
                  row[y_col],
                  ax,
                 color = color)  

# Input df where you have added all var positions to combinatorial and single
def plot_combination_activities(start, ax, df_c, df_s):
    sog1_aa_features = pd.read_csv("../data/Sog1_AA_features.csv")

    single_tile_df = df_c[df_c["Start"] == start]    
    
    all_row = single_tile_df[single_tile_df["var_count"] == 3]
    tile_vars = all_row["vars"].iloc[0]

    # Trio + pairs
    var_combo = list(single_tile_df["vars"])
    activities = list(single_tile_df[activity_col + "_var"])

    # Singles 
    for var in tile_vars:
        row = df_s[(df_s["var"] == var) & (df_s["Start"] == start)]
        var_combo.append(row["var"].iloc[0])
        activities.append(row[activity_col + "_var"].iloc[0])
    
    var_combo.append("WT")
    activities.append(single_tile_df[activity_col + "_wt"].iloc[0])
        
    activities_df = pd.DataFrame({"combo" : var_combo,
                  activity_col : activities})
    activities_df["count"] = activities_df["combo"].astype(str).str.count(",") + 1
    activities_df["combo_str"] = activities_df["combo"].astype(str)
    activities_df.loc[activities_df['combo'] == "WT", 'count'] = 0
    activities_df = activities_df.sort_values(by = "count", ascending = True)
    activities_df = activities_df.reset_index(drop = True)
    
    sns.set_context('talk')

    #display(activities_df)
    palette = sns.color_palette('rocket_r', n_colors = 4)
    palette_dict = {0: palette[0],
                    1: palette[1],
                    2: palette[2], 
                    3: palette[3]}
    sns.set_style("dark")

    sns.barplot(data = activities_df, x = activity_col, y = "combo_str", hue = "count", legend = False, palette = palette_dict, edgecolor = "none", zorder = 2, ax = ax)
    

    ax.axhline(0.5, color = "white", linestyle = "-", zorder = 0)
    ax.axhline(3.5, color = "white", linestyle = "-", zorder = 0)
    ax.axhline(6.5, color = "white", linestyle = "-", zorder = 0)
    sns.despine(ax = ax)
    #plt.ylabel("PS$\\rightarrow$A", rotation = 0, labelpad = 600)
    ax.set_ylabel(str(start) + "-\n" + str(start+40), labelpad = 325, rotation = 0, ha = "center")
    #ax.set_ylabel("")
    ax.set_xlabel("Activity")
    #ax.set_title("Tile " + str(start) + "-" + str(start+40), ha = "center")

    #plt.rcParams['font.family'] = 'monospace'
    
    ax.axvline(activities_df[activity_col].iloc[0], color = palette_dict[0], linestyle = "dotted", alpha = 0.5, zorder = 1, lw = 3)

    for j in activities_df.index:
        if j == 0:
            label = ("".join(sog1_aa_features.iloc[start : start + 40]["aa"]))
            color = palette_dict[0]

        else:
            label = ""
            tile_aas = list(sog1_aa_features.iloc[start : start + 40]["aa"])
            mut_aas = activities_df["combo"].iloc[j]
            #print(mut_aas)

            if isinstance(mut_aas, np.int64):
                mut_aas = np.array([mut_aas])

    
            for i in range(len(tile_aas)):
                if i + start + 1 in mut_aas:
                    #label += "$\\rightarrow$A"
                    label += "A"
                else:
                    label += "-"

            color = palette_dict[activities_df["count"].loc[j]]
        
    
            # if isinstance(mut_aas, np.int64):
            #     mut_aas = np.array([mut_aas])
    
            
            # for i in range(len(tile_aas)):
            #     facecolor='white'
            #     highlight = 'white'
            #     AA = tile_aas[len(tile_aas) - 1 - i]
            #     if type(mut_aas) == str:
            #         color = "white"
            #     elif i + start + 1 in mut_aas:
            #         color = "black"
            #         facecolor = 'red'
            #         AA = "A"
            #     else:
            #         color = "white"
                
            #     plt.text(-1 * spacing * i - spacing * 1.5, j, AA, ha = "center", va = "center", color = color, fontsize = "small", bbox=dict(facecolor=facecolor, 
            #                                                                                                                         alpha=1, 
            #                                                                                                                         edgecolor = 'none', 
            #                                                                                                                        pad = 0.2))
        ax.text(-25, j, label, ha = "right", va = "center", font = "monospace", color = color, fontsize = "x-small")

    #ax.gca().tick_params(axis='y', labelleft=False)
    ax.tick_params(axis='y', labelleft=False)