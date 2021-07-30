# Welcome to the HTA package!

The HTA package can statistically assess the level of both spatial, and global, heterogeneity within a spatial sample. HTA was specifically designed to handle multivariate spatial transcriptomics data, such as Visium samples (10x Genomics), but can be used in other domains (see our paper [1] for further details).

## Prerequisites
Python 3.8 or above.

## Support and bug reports
If you encounter any issues, please contact levyalona at gmail.

# Example 1 - Synthetic Data

### Input format:

Since generating a trait-combination matrix may be complicated, HTA generates it for you from a simpler form of input: a stacked set of matrices where each matrix represents one trait, and each entry indicates if the trait manifests or not (0/1) at the corresponding spatial position. 

For **example**, for 2 traits and a 2D space of 32x32, let's generate some random input to HTA:
	
	from hta.stats import HTA
	
    n_traits = 2
    t_shape = (32, 32, n_traits) 
    t = np.random.random(t_shape)   # random values between 0 and 1
    t = (t > 0.5)*1   # binarising to 0/1

Now `t` contains two stacked matrices, each of shape 32x32. The first represents trait no. 1, the second represents trait no. 2 and each entry has either 0/1 (1 indicates that the trait manifests in that entry).

### HTA:

Once `t` is ready, we can run HTA. 

For **example**, using `t` from above we can decide on a region size of 8 (i.e., each cell in the grid is 8x8) and run the following:

    region_size = 8   # you can also use [8,8] or [8,16] etc., corresp. to [x,y].
    hta = HTA(t, region_size)
    hta_stat, hta_pval = hta.calc()
    print('HTA {:.2f},  p-val: {:.2e}'.format(hta_stat, hta_pval))

### Heterogeneity maps:

To correctly produce the heterogeneity maps, you must specify the name of the traits as they are ordered in `t` (i.e., the trait at index 0 should be the first (index 0) in the stacked matrices in `t`, and the trait at index n_traits-1 should be last in `t`). You can then call the function `hta.plot_heterogeneity_map(..)`. This function will return a pyplot object which can then be modified to your liking (e.g. adding a legend, title etc.). 

For **example**:

    trait_names = ['Trait {}'.format(i+1) for i in range(n_traits)]
    
    hm = hta.plot_heterogeneity_map(trait_names, dot_size=5)
    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    title = 'HTA {:.2f} (p-val: {:.2e}), region_size={}' \
        .format(hta_stat, hta_pval, region_size)
    font_dict = {'family': 'arial', 'size': 11}
    hm.title(title, fontdict=font_dict)
    hm.savefig('../out/result.jpeg', dpi=350)
    hm.close()

### Region report:

To produce the region report mentioned in the paper (which provides additional information on each region) you can use the following code:

    rr = hta.region_report(trait_names)
    rr.to_csv('../out/region_report.csv')

# Example 2 - Visium data

If you are analysing Visium spatial gene expression data (e.g. any of those listed [here](https://support.10xgenomics.com/spatial-gene-expression/datasets), or [this one](https://support.10xgenomics.com/spatial-gene-expression/datasets/1.1.0/V1_Breast_Cancer_Block_A_Section_1) used in our examples below) you can use our `Visium` built-in class to generate the input format for HTA (i.e., the `t` as in Example 1). 

Make sure you have the following folders from your Visium data:

> filtered_feature_bc_matrix
> spatial

and place them in a folder hierarchy as shown below (all shown files are required):
 ```
data_folder
└───filtered_feature_bc_matrix
│   │   barcodes.tsv.gz
│   │   features.tsv.gz
│   │   matrix.mtx.gz
└───spatial
	│   tissue_positions_list.csv
	│   ...
```

You are now ready to load your Visium data and use HTA:

    from hta.stats import HTA  
    from hta.utils import Visium    
    
    path = "../res/data_folder"  # path to 'data folder' in above hierarchy
    trait_names = ['ERBB2', 'CD8A']   # names of features to use in features.tsv.gz  
      
    # load and prepare visium data for HTA  
    visium = Visium(path)  
    visium.load()  
    t, t_mask = visium.prep(trait_names)  
      
    # compute HTA and HTA p-val  
    region_size = 15   # modify region_size as needed
    hta = HTA(t, region_size=region_size, tissue_mask=t_mask) 
    hta_stat, hta_pval = hta.calc()

>What is `t_mask`? It identifies, using barcodes.tsv.gz, which barcodes are under
> the tissue, and is used to discard barcodeds that are not.

Now we can proceed to produce the heterogeneity map and region report (we've left the p-val title formatting in the example below for your convenience, but you can replace it with your own title format):

  
    import math
    
    # generate heterogeneity map and legend
    hm = hta.plot_heterogeneity_map(trait_names, dot_size=8)  
    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=9)  
      
    # set format for p-val and title
    if hta_pval <= 10**-10000:  
        title = 'HTA {:.2f} (p-val ~ 0), region size: {}'.format(hta_stat, region_size)  
    elif hta_pval < 0.00001:  
        p_power_with_base_10 = math.log10(hta_pval)  
        hta_pval = p_power_with_base_10  
        title = 'HTA {:.2f} (p-val: 10^{:.0f}), region size: {}'.format(hta_stat, hta_pval, region_size)  
    else:  
        title = 'HTA {:.2f} (p-val: {:.2f}), region size: {}'.format(hta_stat, hta_pval, region_size)  
    
    # save heterogeneity map 
    font_dict = {'family': 'normal', 'size': 9}  
    hm.title(title, fontdict=font_dict)  
    hm.savefig('../out/{}_hetero_map.jpeg'.format('_'.join(trait_names)), dpi=350)  
    hm.close()  
      
    # save region report  
    rr = hta.region_report(trait_names)  
    rr.to_csv('../out/{}_region_report.csv'.format('_'.join(trait_names)))




# Example 3 - Visium with cluster id per barcode

You can also use HTA with cluster IDs generated per barcode. The best example is using the cluster IDs provided in Visium's analysis folder, but you can use your own cluster IDs, provided they have the same format. 

As an example, you can place Visium's 'analysis' folder (see links above) under your 'data_folder'. The 'analysis' folder contains many k-means clustering results where each barcode has a cluster ID.  

The main differences in the code below compared to Example 2 above are in the lines of code marked with (***):


    from hta.utils import Visium  
    from hta.stats import HTA  
      
    path = "../res/data_folder"  
    k = 10
    clusters_path = '{}/analysis/clustering/kmeans_{}_clusters/clusters.csv'.format(path, k)
    trait_names = [str(i+1) for i in range(k)]  
      
    visium = Visium(path)  
    visium.load()  
    t, t_mask = visium.prep_clusters(clusters_path)  # (***)
    
    # copute HTA
    region_size = 15  
    hta = HTA(t, region_size=region_size, tissue_mask=t_mask)  
    hta_stat, hta_pval = hta.calc()  
      
	# plot heretogeneity map
    hm = hta.plot_heterogeneity_map(trait_names, dot_size=8, is_clusters=True)  # (***)
    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=9)  
    
    # save heterogeneity map 
    title = 'HTA {:.2f} (p-val: {:.2f}), region size: {}'.format(hta_stat, hta_pval, region_size)
    font_dict = {'family': 'normal', 'size': 9}  
    hm.title(title, fontdict=font_dict)  
    hm.savefig('../out/{}_hetero_map.jpeg'.format('_'.join(trait_names)), dpi=350)  
    hm.close()  
      
    # save region report  
    rr = hta.region_report(trait_names)  
    rr.to_csv('../out/{}_region_report.csv'.format('_'.join(trait_names)))

[1] Assessing heterogeneity in spatial data using the HTA index with applications to spatial transcriptomics and imaging. *Levy-Jurgenson et al.* 