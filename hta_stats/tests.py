# from hta_stats.hta import *


def _make_heterogeneous_regions(shape, region_size):
    '''
    Create a sample with perfectly heterogeneous regions
    '''
    M_1, M_2 = np.zeros(shape), np.zeros(shape)
    M_1[::2, :] = 1
    M_2[:, ::2] = 1
    M = np.concatenate([M_1, M_2], axis=-1)
    return M


def _make_homogeneous_regions(shape, region_size):
    '''
    Create a sample with perfectly homogeneous regions
    '''
    M_1, M_2 = np.zeros(shape), np.zeros(shape)
    for i in range(0,shape[0],region_size):
        if (i / region_size) % 2 == 0:
            M_1[i:i+region_size, :] = 1
            M_2[:, i:i+region_size] = 1
    M = np.concatenate([M_1, M_2], axis=-1)
    return M


def run_synthetic(region_size, tensor_gen_fn):
    print("Running synthetic example.")

    trait_names = ['Trait 1', 'Trait 2']
    shapes = [(32, 32)]

    for shape in shapes:

        t = tensor_gen_fn(shape=(shape[0], shape[1], 1), region_size=region_size)

        # calc HTA index and p-val
        hta = HTA(t, region_size)
        hta_stat, hta_pval = hta.calc()

        # plot heterogeneity map
        hm = hta.plot_heterogeneity_map(trait_names, dot_size=10, linecolor='darkslategray')
        hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
        title = 'HTA: {:.2f} (p-val: {:.2e}), region_size: {}' \
            .format(hta_stat, hta_pval, region_size)
        font_dict = {'family': 'arial', 'size': 11}
        hm.title(title, fontdict=font_dict)
        hm.savefig('../out/{}.jpeg'.format(tensor_gen_fn.__name__), dpi=350)
        hm.close()

        # make region report
        rr = hta.region_report(trait_names)
        rr.to_csv('../out/region_report_{}.csv'.format(tensor_gen_fn.__name__))

        print_results(tensor_gen_fn.__name__, hta_stat, hta_pval)


def run_visium_example(region_size=15):
    print("\nRunning visium example")

    trait_names = ['ESR1', 'GATA3']
    run_name = 'visium'

    # load data
    t, t_mask = np.load('../res/t_visium.npy'), np.load('../res/t_visium_tissue_mask.npy')

    # calc HTA index and p-val
    hta = HTA(t, region_size=region_size, tissue_mask=t_mask)
    hti_stat, _, _ = hta._HTI(t)
    hta_stat, hta_pval = hta.calc()

    # plot heterogeneity map
    hm = hta.plot_heterogeneity_map(trait_names, dot_size=5)
    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)
    title = 'HTA {:.2f} (p-val: {:.2e}), region_size={}' \
        .format(hta_stat, hta_pval, region_size)
    font_dict = {'family': 'arial', 'size': 11}
    hm.title(title, fontdict=font_dict)
    hm.savefig('../out/{}.jpeg'.format(run_name), dpi=350)
    hm.close()

    # make region report
    rr = hta.region_report(trait_names)
    rr.to_csv('../out/region_report_{}.csv'.format(run_name))

    print_results(run_name, hta_stat, hta_pval)


def run_visium(gene_exp=True, testing_rand=False, lower_half_of_tissue=False, is_clusters=False, curated_traits_biocarta=False):

    species = 'human'
    if species == 'mouse':
        matrix_dir = "../res/coronal_mouse_kidney"
        trait_names = ['Umod', 'Slc22a12']

    if species == 'human':
        path = "../res/human_breast_cancer/Block_A_sec_1"
        # path = '../res/human_brain_section_1'
        trait_names = ['ERBB2', 'CD8A']
        # trait_names = ['ESR1', 'GATA3']
        # trait_names = ['ESR1', 'GATA3', 'FOXA1']
        # trait_names = ['CD274','ERBB2', 'CD8A']
        # trait_names = ['MYC', 'ESR1',  'ERBB2',  'GATA3',  'FOXA1',  'TP53', 'CDK4']

        if is_clusters:
            trait_names = [str(i+1) for i in range(10)]
            print(trait_names)
            is_clusters=True

        if curated_traits_biocarta:
            df = pd.read_csv('../res/signalling_pathways/expression_sp.txt', sep='\t')
            cols = [b.replace('.', '-') for b in df.columns.to_list()]
            df.columns = cols
            print(df.head(n=10))
            df = df[df.index.isin(['BIOCARTA_CYTOKINE_PATHWAY', 'BIOCARTA_BCR_PATHWAY','BIOCARTA_TCR_PATHWAY'])]#, 'BIOCARTA_DEATH_PATHWAY'])]
            trait_names = df.index.to_list()
            trait_names = [t.split('_')[1] for t in trait_names]
            trait_names = [t if 'BREAST' != t else 'ESR1' for t in trait_names ]
            print(trait_names)

    region_size = 15

    visium = Visium(path)
    visium.load()

    if gene_exp:
        t, t_mask = visium.prep(trait_names)#, combinations)
    elif is_clusters:
        t, t_mask = visium.prep_clusters('{}/analysis/clustering/kmeans_{}_clusters/clusters.csv'.format(path, len(trait_names)))
    elif curated_traits_biocarta:
        t, t_mask = visium.prep_curated_traits(df)

    # keep_combinations = [[0], [1]]
    # keep_combinations = [[i] for i in range(len(trait_names))]
    hta = HTA(t, region_size=region_size, tissue_mask=t_mask)#, keep_combinations=keep_combinations)

    # taking part of the tensor (here, lower half)
    if lower_half_of_tissue:
        half_tensor_mask_aligned = visium._raw_tensor_mask_aligned[visium._raw_tensor_mask_aligned.shape[0] // 2:, :, :]
        half_tensor_mask_aligned = round_tensor_shape(half_tensor_mask_aligned)
        t, t_mask = half_tensor_mask_aligned[:, :, :-1], half_tensor_mask_aligned[:, :, -1]
        visium.trait_tensor, visium.tissue_mask = t, t_mask
        visium._binarize(trait_names, np.median)
        hta = HTA(visium.trait_tensor, region_size=region_size, tissue_mask=t_mask)

    if testing_rand:
        name = 'test_rand'
        t = hta._get_rand_perm_tensor()  # accurate for CLT but computationally expensive
        # t = apply_tissue_mask(t_visium, visium.tissue_mask)
        hta = HTA(t, region_size=region_size, tissue_mask=t_mask)

    hti_stat, _, _ = hta._HTI(t)  # TODO: hti static and can be accessed as attribute
    hta_stat, hta_pval_clt = hta.calc()

    hm = hta.plot_heterogeneity_map(trait_names, dot_size=8, is_clusters=is_clusters)

    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=9)
    # hm.legend(loc='lower left', ncol=3)
    p = hta_pval_clt

    # format p-value in title
    if p <= 10**-10000:
        title = 'HTA {:.2f} (p-val ~ 0), region size: {}'.format(hta_stat, region_size)
    elif p < 0.00001:
        p_power_with_base_10 = math.log10(p)
        p = p_power_with_base_10
        title = 'HTA {:.2f} (p-val: 10^{:.0f}), region size: {}'.format(hta_stat, p, region_size)
    else:
        title = 'HTA {:.2f} (p-val: {:.2f}), region size: {}'.format(hta_stat, p, region_size)

    font_dict = {'family': 'normal', 'size': 9}
    hm.title(title, fontdict=font_dict)
    hm.savefig('../out/visium_{}_{}.jpeg'.format(title, '_'.join(trait_names)), dpi=350)
    hm.close()

    rr = hta.region_report(trait_names)
    rr.to_csv('../out/region_report_{}_{}.csv'.format('visium', '_'.join(trait_names)))

def run_visium_clusters():
    from hta_stats.utils import Visium
    from hta_stats.hta import HTA

    path = "../res/human_breast_cancer/Block_A_sec_1"
    k = 10
    clusters_path = '{}/analysis/clustering/kmeans_{}_clusters/clusters.csv'.format(path, k)
    trait_names = [str(i + 1) for i in range(k)]

    visium = Visium(path)
    visium.load()
    t, t_mask = visium.prep_clusters(clusters_path)

    # compute HTA
    region_size = 15
    hta = HTA(t, region_size=region_size, tissue_mask=t_mask)
    hta_stat, hta_pval = hta.calc()

    # plot heretogeneity map
    hm = hta.plot_heterogeneity_map(trait_names, dot_size=8, is_clusters=True)  # Note the 'is_clusters = True'
    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=9)

    # save heterogeneity map
    title = 'HTA {:.2f} (p-val: {:.2f}), region size: {}'.format(hta_stat, hta_pval, region_size)
    font_dict = {'family': 'normal', 'size': 9}
    hm.title(title, fontdict=font_dict)
    hm.savefig('../out/{}_hetero_map.jpeg'.format('_'.join(trait_names)), dpi=350)
    hm.close()


def run_visium_simple_for_readme():
    from hta_stats.hta import HTA
    from hta_stats.utils import Visium

    path = "../res/human_breast_cancer/Block_A_sec_1"  # path to 'data folder' in above hierarchy
    trait_names = ['ERBB2', 'CD8A']  # names of features to use in features.tsv.gz

    # load and prepare visium data for HTA
    visium = Visium(path)
    visium.load()
    t, t_mask = visium.prep(trait_names)

    # compute HTA and HTA p-val
    region_size = 15  # modify region_size as needed
    hta = HTA(t, region_size=region_size, tissue_mask=t_mask)
    hta_stat, hta_pval = hta.calc()

    import math

    # generate heterogeneity map and legend
    hm = hta.plot_heterogeneity_map(trait_names, dot_size=8)
    hm.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=9)

    # set format for p-val and title
    if hta_pval <= 10 ** -10000:
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

    hti_stat, _, _ = hta._HTI(t)

if __name__ == '__main__':

    # run_synthetic(region_size=8, tensor_gen_fn=_make_homogeneous_regions)
    # run_synthetic(region_size=8, tensor_gen_fn=_make_heterogeneous_regions)
    #
    # run_visium(gene_exp=True)
    # run_visium(gene_exp=True, lower_half_of_tissue=True)
    # run_visium(gene_exp=True, testing_rand=True)
    # run_visium(gene_exp=False, is_clusters=True)
    run_visium_clusters()
    # run_visium(gene_exp=False, curated_traits_biocarta=True)
    # run_visium_simple_for_readme()

    # run_visium_example()