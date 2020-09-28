
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

A = pd.read_excel("ATHENA-PGs-ABC1Ks-forLalit.xlsx",sheet_name=None,skiprows=1)
show_plots = False

# # -- prepare the data by normalizing each protein --
# for key in A.keys():
#     print("\n",key,"\n",sep="")
#     df = A[key].drop(columns=['our interests','group','lab annotation']).set_index('Accession')
# 
#     # remove proteins that have all zero values (since they create problems for the correlation-based distance matrix)
#     if any(df.isna().sum(axis=1) == df.shape[1]):
#         print(f"dropping {sum(df.isna().sum(axis=1) == df.shape[1])} all-missing row(s)")
#         df = df.loc[df.isna().sum(axis=1) < df.shape[1],:]
# 
#     # fill missing values with zero
#     df = df.fillna(0)
# 
#     # scale each row by its max (i.e. by the max across all tissues)
#     df_s = df.apply(lambda x: x/x.max(), axis=1)
#     df_s.to_csv(f"{key}_Scaled_data.csv")
# 
#     # note: dividing by the max seems better suited (instead of min-max scaling) because:
#     # - we don't want to replace the minimum (across tissues) in each protein with zero
#     # - it will be roughly the same as min-max scaling since we're using correlation as the distance measure
# 
#     # # min-max scaling of each row
#     # df_mm = df.apply(lambda x: (x-x.min())/(x.max()-x.min()), axis=1)
# 
#     # apply z-score transformation for each row
#     df_z = df.apply(lambda x: (x-x.mean())/x.std(ddof=1), axis=1)
#     df_z.to_csv(f"{key}_Zscore_data.csv")


# -- plot the abundance without normalization --
for key in A.keys():
    print("\n",key,"\n",sep="")
    df = A[key].drop(columns=['our interests','group','lab annotation']).set_index('Accession')

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
    # arrange in order of tissue-based clusters: see "cluster the tissues" in clust.R
    if key == 'PGs':
        # get protein_order that matches dendrogram: run pvclust() as in clust.R section "cluster the proteins (show p-values)" and then cat("'",paste0(pv$hclust$labels[pv$hclust$order],collapse="','"),"'",sep="")
        protein_order = ['AT5G41120.1','AT3G24190.1','AT1G17420.1','AT1G73750.1','AT4G32770.1','AT5G13800.1','AT1G28150.1','AT4G31390.1','AT1G54570.1','AT4G22240.1','AT3G07700.1','AT3G27110.1','AT1G52590.1','AT1G06690.1','AT2G42130.1','AT1G32220.1','AT4G38970.1','AT3G43540.1','AT2G46910.1','AT2G34460.1','AT1G71810.1','AT1G79600.1','AT5G05200.1','AT4G04020.1','AT2G35490.1','AT4G13200.1','AT3G23400.1','AT3G58010.1','AT2G41040.1','AT3G10130.1','AT4G19170.1','AT1G78140.1','AT5G08740.1']
        # get tissue_order that matches dendrogram: run clust.R "cluster the tissues (no p-values)" and copy-paste the tissue names in order
        tissue_order = ['root','callus','egg like callus','cell culture early','cell culture late','carpel','silique','flower pedicle','root tip','root upper zone','sepal','stamen','petal','flower','cotelydons','leaf petiole','cauline leaf','shoot tip','leaf distal','leaf proximal','hypocotyl','node','internode','embryo','seed','seed imbibed','pollen','senescent leaf','silique septum','silique valves']
        df = df.loc[protein_order,tissue_order]
    elif key == '17-ABC1Ks':
        protein_order = ['AT5G24970.1','AT1G61640.1','AT2G40090.1','AT5G24810.1','AT5G50330.1','AT1G11390.1','AT4G01660.1','AT1G65950.1','AT5G64940.1','AT3G24190.1','AT2G39190.1','AT4G31390.1','AT4G24810.1','AT1G71810.1','AT3G07700.1','AT1G79600.1','AT5G05200.1']
        tissue_order = ['flower','carpel','petal','stamen','egg like callus','cell culture early','callus','cell culture late','pollen','root upper zone','root','root tip','seed','seed imbibed','embryo','cauline leaf','sepal','senescent leaf','silique septum','silique valves','cotelydons','internode','hypocotyl','node','leaf distal','silique','flower pedicle','leaf petiole','leaf proximal','shoot tip']
        df = df.loc[protein_order,tissue_order]

    fig,axs = plt.subplots(figsize=(12,12))
    sns.heatmap(df, cmap="icefire", annot=False, ax=axs)
    axs.set_title(key)
    if show_plots:
        plt.show()
    else:
        plt.savefig(f"{key}_heatmap.png")

