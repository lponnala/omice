
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

A = pd.read_excel("ATHENA-PGs-ABC1Ks-forLalit.xlsx",sheet_name=None,skiprows=1)
show_plots = False

# normalize each protein
for key in A.keys():
    print("\n",key,"\n",sep="")
    df = A[key].drop(columns=['our interests','group','lab annotation']).set_index('Accession')

    # remove proteins that have all zero values (since they create problems for the correlation-based distance matrix)
    if any(df.isna().sum(axis=1) == df.shape[1]):
        print(f"dropping {sum(df.isna().sum(axis=1) == df.shape[1])} all-missing row(s)")
        df = df.loc[df.isna().sum(axis=1) < df.shape[1],:]

    # fill missing values with zero
    df = df.fillna(0)

    # scale each row by its max (i.e. by the max across all tissues)
    df_s = df.apply(lambda x: x/x.max(), axis=1)
    df_s.to_csv(f"{key}_Scaled_data.csv")

    # note: dividing by the max seems better suited (instead of min-max scaling) because:
    # - we don't want to replace the minimum (across tissues) in each protein with zero
    # - it will be roughly the same as min-max scaling since we're using correlation as the distance measure

    # # min-max scaling of each row
    # df_mm = df.apply(lambda x: (x-x.min())/(x.max()-x.min()), axis=1)

    # apply z-score transformation for each row
    df_z = df.apply(lambda x: (x-x.mean())/x.std(ddof=1), axis=1)
    df_z.to_csv(f"{key}_Zscore_data.csv")

    # # plot as bars
    # fig,axs = plt.subplots(figsize=(12,12))
    # df.T.plot(kind='bar',ax=axs)
    # axs.set_title(key)
    # if show_plots:
    #     plt.show()
    # else:
    #     plt.savefig(f"{key}_bars.png")

    # # plot as lines
    # fig,axs = plt.subplots(figsize=(12,12))
    # df.T.plot(kind='line',ax=axs)
    # axs.set_title(key)
    # if show_plots:
    #     plt.show()
    # else:
    #     plt.savefig(f"{key}_lines.png")

    # plot as heatmap
    fig,axs = plt.subplots(figsize=(12,12))
    sns.heatmap(df, cmap="icefire", annot=False, ax=axs)
    axs.set_title(key)
    if show_plots:
        plt.show()
    else:
        plt.savefig(f"{key}_heatmap.png")

