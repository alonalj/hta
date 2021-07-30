
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from collections import Counter

from hta.utils import apply_tissue_mask

CLT_REPEAT = 1000
EPSILON = 1e-100

'''
This code is the implementation of the method described in the paper:

"Assessing heterogeneity in spatial data using the HTA index with applications to spatial transcriptomics and imaging"
by Levy-Jurgenson et al.
'''


def _powerset(iterable):
    '''
    Create the set of all trait combinations.
    '''
    from itertools import chain, combinations
    xs = list(iterable)
    return list(chain.from_iterable(combinations(xs, n) for n in range(len(xs) + 1)))


def _np_coords_to_xyz(coords, shape):
    '''
    coords: the numpy region position in the matrix (the row, column... indices)

    Translate numpy indices to x,y,z region bounds (due to the fact that, e.g., x axis is respresented by columns, which
    is the second np index in a numpy matrix)
    '''
    if len(coords) == 2:
        return coords[1] - 0.5, shape[0] - coords[0] - 0.5
    return coords[1] - 0.5, shape[0] - coords[0] - 0.5, coords[2] - 0.5


def _xy_region_size_for_np(coords):
    '''
    coords: the shape of the region size, e.g. for 2d a square region size of 8 would be (8,8).

    Translates x,y,z region sizes to numpy matrix representation (e.g. x axis is actually respresented by columns, which
    is the second np index in a numpy matrix)
    '''
    a = [coords[1], coords[0]]
    a.extend(coords[2:])
    return a


def _get_rand_tensor_for_shape(shape, n_traits, region_size, p):
    '''
    Generates a random sample for the given number of traits and shape, using probability p of observing the trait
    '''
    M = []
    for i in range(n_traits):
        _M = np.random.choice(np.array([1, 0]), size=shape, p=[p[i], 1 - p[i]])
        if len(M) == 0:
            M = _M
        else:
            M = np.concatenate([M, _M], axis=-1)
    return M


def _get_comb_ids_kept(keep_combinations):
    '''Handles user input if they choose a subset of trait-combinations to be kept'''
    combination_ids_kept = []
    for tc in keep_combinations:
        trait_comb_id = sum([2 ** i for i in tc])
        combination_ids_kept.append(trait_comb_id)
    return combination_ids_kept


class HTA:
    '''
    Class used to compute HTA.
    '''

    def __init__(self, t, region_size, n_repeat=CLT_REPEAT, keep_combinations=[], tissue_mask=[]):
        '''
        t: the trait tensor -- a numpy ndarray where the first two or three dimensions
        represent the spatial dimensions (i.e., the x, y, z axes). The remaining dimensions, each, represent a single trait,
        with binary values that depict whether that trait manifests in that x,y,z position (there could be more if e.g., a
        time dimension is added, or similar, but we did not demonstrate this in the paper), e.g. (x,y, n_traits) or
        (x,y, z, n_traits) for 3D like MRI.
        For example, in a 2d space, a trait_tensor of shape (16,16,4) describes a 16x16 space (x,y) and in each such
        (x,y) position, there are 4 traits. The element (x,y,0), would have 1 or 0, depending on whether the first trait out
        of four manifests in this (x,y) position or not, respectively.

        region_size: the size of the region in which HTIs will individually computed. HTA is the average of HTIs across
        all regions.

        n_repeat: the number of repeats with which to perform the random permutations when computing HTA p-val (see paper)

        keep_combinations: for future, not currently available

        tissue_mask: a binary tissue mask of the same spatial dimensions as the t (e.g. if t is (16,16,4) for
        (x,y, n_traits) then tissue_mask is (16,16)) where 0 is for disregarding that value (e.g. if is expected to be
        a 'noise' measurement, like a barcoded spot that has not tissue under it).

        '''
        self._t = t
        self._set_region_size(region_size)
        self._n_repeat = n_repeat
        self._combination_ids_kept = _get_comb_ids_kept(keep_combinations)
        self._tissue_mask = tissue_mask
        self._n_trait_combs_actual, self._trait_comb_to_linear_id = self._calc_trait_comb_info()
        self._t_regions, self._t_weights, self._t_regions_ijk_nonempty = None, None, None
        self._t_hti_per_region = None
        self._t_comb_counts_per_region = None
        self._is_counts = False
        self._is_3d = len(t.shape) > 3

    def _set_region_size(self, region_size):
        if isinstance(region_size, int) or isinstance(region_size, float):
            self._region_size = [region_size, region_size, region_size]
        else:
            assert len(region_size) > 1, "Expected a list of length > 1, e.g. [region_size_x, region_size_y], " \
                                               "but got {}".format(self._region_size)
            self._region_size = _xy_region_size_for_np(region_size)

    def _check_is_valid_input(self):
        assert self._t.shape[-1] >= 2, "Illegal number of traits. Input must contain at least 2 traits, but received " \
                                       "shape: {}, with number of traits {}.".format(self._t.shape, self._t.shape[-1])

    def _check_is_legal_region_size_for_shape(self):
        assert (len(self._t.shape) > 3 and len(self._region_size) >= 3) or (len(self._t.shape) <= 3), "Shape of input implies > 3D data, but you supplied less than 3 values of region sizes."

    def calc(self):
        ''' Main function for computing HTA. Handles both HTA index value and p-value.
        returns: hta and Lyapunov CLT p-value
        '''

        self._t_regions, self._t_weights, self._t_regions_ijk_nonempty = self._get_region_data_from_tensor(self._t)

        get_stat_clt, get_pval_clt = self._get_CLT_hta_p_val()

        hta, hti_per_region, trait_comb_ids_per_region, trait_comb_ids_to_counts_per_region = self._HTA(self._t_regions,
                                                                                                        self._t_weights)
        self._t_hti_per_region = hti_per_region
        self._t_comb_counts_per_region = trait_comb_ids_to_counts_per_region

        sample_mean = hta
        stat_clt = get_stat_clt(sample_mean)

        p_val_clt = get_pval_clt(stat_clt)
        return hta, p_val_clt

    def _get_region_data_from_tensor(self, t, region_indices_to_include=[]):
        '''
        Returns the data per region in t, and the weights of each region (based on the number of samples in the region)

        t: the trait tensor -- a numpy ndarray where the first two or three dimensions
        represent the spatial dimensions (i.e., the x, y, z axes). The remaining dimensions, each, represent a single trait,
        with binary values that depict whether that trait manifests in that x,y,z position (there could be more if e.g., a
        time dimension is added, or similar, but we did not demonstrate this in the paper), e.g. (x,y, n_traits) or
        (x,y, z, n_traits) for 3D like MRI.
        For example, in a 2d space, a trait_tensor of shape (16,16,4) describes a 16x16 space (x,y) and in each such
        (x,y) position, there are 4 traits. The element (x,y,0), would have 1 or 0, depending on whether the first trait out
        of four manifests in this (x,y) position or not, respectively.

        region_indices_to_include: will be supported in the futute -- will enable users to select specific trait combinations.

        '''
        def _process_region(region_data, weights, regions_i_j_with_data, region_id):
            n_any_positive_in_region = np.count_nonzero(np.sum(region, axis=-1))  # e.g. number of tissue barcodes
            if n_any_positive_in_region > 0:
                region_data.append(region)
                weights.append(n_any_positive_in_region)
                regions_i_j_with_data.append(region_id)
            elif region_id in region_indices_to_include:
                region_data.append(region)
            return region_data, weights, regions_i_j_with_data

        region_data, weights, regions_i_j_with_data = [], [], []
        n_any_positive = np.count_nonzero(np.sum(t, axis=-1))
        for i in range(0, t.shape[0], self._region_size[0]):
            for j in range(0, t.shape[1], self._region_size[1]):
                if self._is_3d:
                    for k in range(0, t.shape[2], self._region_size[2]):
                        region = t[i:i + self._region_size[0], j:j + self._region_size[1], k:k + self._region_size[2], ...]
                        region_data, weights, regions_i_j_with_data = _process_region(region_data, weights, regions_i_j_with_data, (i,j,k))
                else:
                    region = t[i:i + self._region_size[0], j:j + self._region_size[1], ...]
                    region_data, weights, regions_i_j_with_data = _process_region(region_data, weights,
                                                                                  regions_i_j_with_data, (i,j))

        assert sum(weights) == n_any_positive
        return region_data, np.array(weights) / n_any_positive, regions_i_j_with_data

    def _calc_trait_comb_info(self):
        '''
        Produces trait combination ids, also used for plotting purposes later (to assign each combination a different color)
        returns: n_trait_combs_actual -- the actual number of trait combinations present in the data
                trait_comb_to_linear_id -- a number representing the trait combination, running from 0 to n_trait_combs_actual
        '''
        n_traits = self._t.shape[-1]
        trait_comb_ids = list(set(sum([self._t[...,i] * 2 ** i for i in range(n_traits)]).flatten()))
        trait_comb_to_linear_id = {trait_comb_ids[i]:i for i in range(len(trait_comb_ids))}
        if len(self._combination_ids_kept) > 0:
            n_trait_combs_actual = len(self._combination_ids_kept)
        else:
            n_trait_combs_actual = len(trait_comb_ids)-1*(0 in trait_comb_to_linear_id.keys())  # -1 if empty combination in trait_comb_to_linear_id.keys()
        # print("actual N trait combs = {}".format(n_trait_combs_actual))
        return n_trait_combs_actual, trait_comb_to_linear_id

    def _HTI(self, x):
        '''
        Computes the heterogeneity index (HTI). 1 = heterogeneous, 0 = homogeneous.
        x: is a tensor molecular cartography, where each layer corresponds to a single trait. Expects a numpy array
        of shape: X x Y x N_TRAITS containing binary values.
        Example:
        The tensor molecular cartogrpahy will be reshaped to: N_LOCATIONS X N_TRAITS (i.e. X*Y X N_TRAITS), e.g.
        data = np.array([[0,1], [1,1], [1,0]]) represents three locations and two traits. The first location only has the
        second trait, the second location has both traits, and the third location has only the first trait. HTI will be 1.

        :return: HTI
        '''
        from scipy.stats import entropy
        n_traits = x.shape[-1]
        x = x.astype(int)
        n_trait_combinations_actual = self._n_trait_combs_actual  # 2**n_traits-1  # TODO: move  out
        if self._is_counts:
            # TODO: future
            # already in format: 1 X N_TYPES where the i-th element holds # locations with the i-th type
            n_positive_loc = np.count_nonzero(x)
        else:
            # turn into format: N_LOCATIONS X N_TYPES the i-th element indicates 0/1 if the location has the i-th type
            x = x.reshape(np.prod(x.shape[:-1]), x.shape[-1])
            trait_comb_ids = sum([x[:, i] * 2 ** i for i in range(n_traits)])
            if len(self._combination_ids_kept) > 0:
                trait_comb_ids = [i for i in trait_comb_ids if i in self._combination_ids_kept]
            n_positive_loc = np.count_nonzero(trait_comb_ids)
            trait_comb_ids = [t for t in trait_comb_ids if t != 0]
            trait_comb_to_count = Counter(trait_comb_ids)
            trait_comb_counts = np.array(list(trait_comb_to_count.values()))
            _x = trait_comb_counts

        if n_positive_loc == 0:
            # all positions are empty - trivially homogeneous
            return 0, trait_comb_ids, trait_comb_to_count

        probs = _x / n_positive_loc

        if probs[0] == 1:
            # all positions are identical - trivially homogeneous
            return 0, trait_comb_ids, trait_comb_to_count

        hti = entropy(probs, base=n_trait_combinations_actual)
        return hti, trait_comb_ids, trait_comb_to_count

    def _HTA(self, regions, weights):
        '''
        Given the regions' data and corresponding region weights, computes the HTA (weighted_mean_hti) and returns
        additional details including hti_per_region, trait_comb_ids_per_region, trait_comb_counts_per_region.

        regions: a list of the regions' data (flattened submatrix of t).
        weights: a list of the number of elements in the regions (elements = have some positive trait) divided by the total number of elements in t

        '''
        hti_per_region, trait_comb_ids_per_region, trait_comb_counts_per_region = [], [], []
        assert len(weights) == len(regions)
        for r in range(len(regions)):
            region = regions[r]
            region_hti, region_trait_comb_ids, region_trait_comb_count = self._HTI(region)
            hti_per_region.append(region_hti)
            trait_comb_ids_per_region.append(region_trait_comb_ids)
            trait_comb_counts_per_region.append(region_trait_comb_count)
        weighted_mean_hti = np.dot(hti_per_region, weights)
        return weighted_mean_hti, hti_per_region, trait_comb_ids_per_region, trait_comb_counts_per_region

    def _get_direct_hta_p_val(self):
        '''
        Compute HTA p-value by direct (monte-carlo) simulations. Less accurate than Lyapunov CLT unless
         the number of regions is small (e.g. 4 regions).
        '''
        hta_scores_all_runs = []
        for n in range(self._n_repeat):
            t_rand = self._get_rand_perm_tensor()
            if len(self._tissue_mask) > 0:
                t_rand = apply_tissue_mask(t_rand, self._tissue_mask)
            regions, _, _ = self._get_region_data_from_tensor(t_rand, self._t_regions_ijk_nonempty)
            hta, _, _, _ = self._HTA(regions, self._t_weights)
            hta_scores_all_runs.append(hta)

        p_direct_fn = lambda HTA: sum(np.array(hta_scores_all_runs) <= HTA) / len(hta_scores_all_runs)
        p_val_direct = p_direct_fn
        return p_val_direct

    def _get_CLT_hta_p_val(self):
        '''
        Compute HTA p-value based on Lyapunov CLT. Returns the statistic and the pvalue.
        '''
        whti_scores_all_runs_all_regions, weights_all_runs_all_regions = [], []
        for n in range(self._n_repeat):
            t_rand = self._get_rand_perm_tensor()
            if len(self._tissue_mask) > 0:
                t_rand = apply_tissue_mask(t_rand, self._tissue_mask)

            regions, _, _ = self._get_region_data_from_tensor(t_rand, self._t_regions_ijk_nonempty)
            estimated_mean_hti_per_region, hti_per_region, _, _ = self._HTA(regions, self._t_weights)
            weighted_hti_per_region = list(hti_per_region * self._t_weights)
            whti_scores_all_runs_all_regions.append(weighted_hti_per_region)  # rows = runs, cols = regions (1,..,k)

        # Weighting for Lyapunov
        whti_scores_all_runs_all_regions = np.array(whti_scores_all_runs_all_regions)
        estimated_mean_hti_per_region = np.mean(whti_scores_all_runs_all_regions, axis=0)
        estimated_var_hti_per_region = np.var(whti_scores_all_runs_all_regions, axis=0,ddof=1)
        sum_mu_regions = sum(estimated_mean_hti_per_region)  # sum_k=1^n mu_k (n=no. regions, mu_k is mean estimate for region k)
        sum_var_regions = sum(estimated_var_hti_per_region)
        sn = np.sqrt(sum_var_regions)
        get_stat_clt = lambda whti: (whti - sum_mu_regions) / (sn + EPSILON)
        get_pval_clt = lambda stat: stats.norm.cdf(stat)

        return get_stat_clt, get_pval_clt

    def _get_rand_perm_tensor(self):
        '''
        Generate a sample from the null-model (i.e. a random uniform permutation of the data, in-place, i.e. maintaining
        spatial coordinates that have data so that the region weights remain the same).

        Returns t_perm (a randomly permuted t).
        '''
        if self._is_3d:
            t_perm = np.zeros_like(self._t)
            rows, cols, depths = np.where(self._t.any(axis=-1))
            perm = np.random.permutation(len(rows))
            t_perm[rows[perm], cols[perm], depths[perm]] = self._t[rows, cols, depths]
        else:
            t_perm = np.zeros_like(self._t)
            rows, cols = np.where(self._t.any(axis=-1))
            perm = np.random.permutation(len(rows))
            t_perm[rows[perm], cols[perm]] = self._t[rows, cols]
        assert np.sum(self._t) == np.sum(t_perm), "_t and t_perm have different sums."

        return t_perm

    def plot_heterogeneity_map(self, trait_names, dot_size=3, t=None, with_grid=True,linecolor='gray', is_clusters=False):
        def scatter_2d(t):
            for i in range(0, t.shape[0]):
                for j in range(0, t.shape[1]):
                    x = j
                    y = t.shape[0] - i - 1
                    if is_clusters:
                        trait_comb_id = np.argmax(t[i, j, ..., :] == 1) + 1
                    else:
                        trait_comb_id = sum([2 ** k * t[i, j, ..., k] for k in range(t.shape[-1])])
                    if len(self._combination_ids_kept) > 0 and trait_comb_id not in self._combination_ids_kept:
                        continue
                    trait_comb = ' '.join(trait_names[k] for k in range(self._t.shape[-1]) if t[i, j, k] == 1)
                    if is_clusters:
                        c = linear_id_to_color[trait_comb_id]
                    else:
                        c = linear_id_to_color[trait_comb_id_to_linear_id[trait_comb_id]]
                    row_list.append([x, y, c, trait_comb])
            df = pd.DataFrame(row_list, columns=['x', 'y', 'color', 'trait_comb'])
            fig, ax = plt.subplots()
            for trait_comb in list(set(df['trait_comb'])):
                if trait_comb == '':
                    continue
                subset = df[df['trait_comb'] == trait_comb]
                plt.scatter(subset.x, subset.y, c=subset.color, s=dot_size,
                            label=trait_comb)  # , alpha=0.8)# cmap=plt.cm.get_cmap('Greens', M.shape[-1]))
                # plt.colorbar()
            grid_loci_y = [t.shape[0]- i * self._region_size[0] - 0.5 for i in range(t.shape[0] // self._region_size[0] + 1)]
            grid_loci_x = [i * self._region_size[1] - 0.5 for i in range(t.shape[1] // self._region_size[1] + 1)]
            ax.set_xticks(grid_loci_x, minor=False)
            ax.set_yticks(grid_loci_y, minor=False)
            if with_grid:
                ax.grid(True, which='major', axis='both', color=linecolor)
            plt.ylim([-0.5 - 0.05, t.shape[0] - 0.5 + 0.05])
            plt.xlim([-0.5 - 0.05, t.shape[1] - 0.5 + 0.05])
            for tick in ax.get_xticklabels():
                tick.set_rotation(45)

            # shrink box figure so that the legend fits
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.2,
                             box.width *0.8, box.height * 0.8])

            print("number of combinations = {}".format(n_trait_combs_nonempty)) #TODO: remove
            if n_trait_combs_nonempty <= 15:
                plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.15), ncol=n_trait_combs_nonempty//2)
            handles, labels = ax.get_legend_handles_labels()
            # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0])) # if you want to sort the legend
            # ax.legend(handles, labels)
            return plt

        def scatter_3d(t):
            from mpl_toolkits.mplot3d import Axes3D
            assert is_clusters is False, "Clusters not yet supported in 3D. See the 2D case to adapt."
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.dist = 13
            for l in reversed(range(0, t.shape[2])):
                for i in range(0, t.shape[0]):
                    for j in range(0, t.shape[1]):
                        x = j
                        y = t.shape[0] - i - 1
                        z = l
                        trait_comb_id = sum([2 ** k * t[i, j, l, k] for k in range(t.shape[-1])])
                        if len(self._combination_ids_kept) > 0 and trait_comb_id not in self._combination_ids_kept:
                            continue
                        trait_comb = ' '.join(trait_names[k] for k in range(t.shape[-1]) if t[i, j, l, k] == 1)
                        c = linear_id_to_color[trait_comb_id_to_linear_id[trait_comb_id]]
                        row_list.append([x, y, z, c, trait_comb])

            df = pd.DataFrame(row_list, columns=['x', 'y', 'z', 'color', 'trait_comb'])
            for trait_comb in list(set(df['trait_comb'])):
                if trait_comb == '':
                    continue
                subset = df[df['trait_comb'] == trait_comb]
                ax.scatter(subset.x, subset.y, subset.z, c=subset.color, s=dot_size,
                            label=trait_comb, alpha=0.8, depthshade=True)# cmap=plt.cm.get_cmap('Greens', M.shape[-1]))
                # # plt.colorbar()
                # grid_loci = [i * self._region_size - 0.5 for i in
                #              range(max(t.shape[1], t.shape[0]) // self._region_size + 1)]
                # ax.set_xticks(grid_loci, minor=False)
                # ax.set_yticks(grid_loci, minor=False)
                # if with_grid:
                #     ax.grid(True, which='major', axis='both')
            handles, labels = ax.get_legend_handles_labels()
            for h in handles:
                h.set_sizes([20])
            # # sort both labels and handles by labels
            # labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
            # ax.legend(handles, labels)
            plt.xlim([-0.5 - 0.05, t.shape[0] - 0.5 + 0.05])
            plt.ylim([-0.5 - 0.05, t.shape[0] - 0.5 + 0.05])
            for tick in ax.get_xticklabels():
                tick.set_rotation(45)
            plt.xlabel('x', labelpad=10)
            plt.ylabel('y', labelpad=10)
            # ax.yaxis.grid(True, which='minor')
            # plt.show()
            print("number of combinations = {}".format(n_trait_combs_nonempty))
            if n_trait_combs_nonempty <= 15:
                plt.legend(loc='best')#, ncol=1)
            return plt

        if t is None:
            t = self._t
        # # M = apply_empty_mask(M)
        # assert len(t.shape) == 3, "Heterogeneity map supported only for 2D spatial data: (x,y,n_traits). " \
        #                           "Got shape: {}".format(t.shape)
        row_list = []
        added_colors = []  # [(182, 0, 14),(235, 140, 1), (181, 191, 1), (1, 181, 183)]#, (182, 140, 184), (204, 90, 113), (128, 71, 94), (238, 123, 48)]
        # added_colors = [(240, 174, 205), (203, 235, 170)]
        # # adding mean between first two non-whites
        # added_colors.append(((added_colors[1][0]+added_colors[2][0]) / 2, (added_colors[1][1]+added_colors[2][1]) / 2, (added_colors[1][2]+added_colors[2][2]) / 2))
        # added_colors = [(c[0]/255., c[1]/255., c[2]/255.) for c in added_colors]
        # cmap = plt.get_cmap('tab20')
        cmap = plt.get_cmap('tab10')

        # cmap = plt.get_cmap('PiYG',11)
        cmap_list = [cmap(i) for i in range(cmap.N)]
        # np.random.shuffle(cmap_list)
        cmap_list.extend(added_colors)
        cmap_list.reverse()
        linear_id_to_color = {}
        n_trait_combs_nonempty, trait_comb_id_to_linear_id = self._n_trait_combs_actual, self._trait_comb_to_linear_id

        if len(cmap_list) < n_trait_combs_nonempty + 1:
            # cmap = plt.get_cmap('gist_rainbow')
            cmap = plt.get_cmap('rainbow')
            cmap_list = [cmap(i) for i in range(cmap.N)]
            cmap_list_extension = [tuple(np.add(cmap_list[i], cmap_list[i + 1]) / 2) for i in range(len(cmap_list) - 1)]
            cmap_list.extend(cmap_list_extension)
            while len(cmap_list) < n_trait_combs_nonempty+1:
                cmap_list_extension = [tuple(np.add(cmap_list_extension[i],cmap_list_extension[i+1])/2) for i in range(len(cmap_list_extension)-1)]
                cmap_list.extend(cmap_list_extension)
            assert len(cmap_list) >= n_trait_combs_nonempty + 1, "Number of trait combinations is {} which larger than" \
                                                                 " available colors. Apply additional colormap extentions.".format(len(cmap_list))
            n_per_slice = len(cmap_list) // (n_trait_combs_nonempty + 1)
            new_cmaps_list = []
            for i in range(n_trait_combs_nonempty + 1):
                current_slice = cmap_list[n_per_slice * i: n_per_slice * (i + 1)]
                new_cmaps_list.append(current_slice[n_per_slice // 2])
            cmap_list = new_cmaps_list

        for linear_id in range(n_trait_combs_nonempty + 1):
            linear_id_to_color[linear_id] = cmap_list.pop()

        if len(t.shape) == 3:
            return scatter_2d(t)
        else:
            return scatter_3d(t)

    def region_report(self, trait_names):
        '''
        Generates the region report, describing for each region (represented by its upper left (x.y) corner) which
        train combinations are present and how many, sorted from highest count to lowest.
        '''
        def map_comb_id_to_name(trait_comb_ids, trait_comb_names):
            for i in range(len(trait_comb_names)):
                s = 0
                for j in trait_comb_ids[i]:
                    s += 2 ** j
                comb_id_to_comb_name[s] = trait_comb_names[i]
            return comb_id_to_comb_name

        trait_ix = list(range(len(trait_names)))
        region_loc, hti, comb_most_to_least_frequent = [], [], []
        comb_id_to_comb_name = {}
        trait_comb_ids, trait_comb_names = _powerset(trait_ix), _powerset(trait_names)
        comb_id_to_comb_name = map_comb_id_to_name(trait_comb_ids, trait_comb_names)
        # fill in region info
        for i in range(len(self._t_regions_ijk_nonempty)):
            comb_ids_decr_frequency = []
            region_loc.append(_np_coords_to_xyz(self._t_regions_ijk_nonempty[i], self._t.shape))
            hti.append(self._t_hti_per_region[i])
            for comb in self._t_comb_counts_per_region[i].most_common():
                comb_name, comb_count = comb_id_to_comb_name[comb[0]], comb[1]
                comb_ids_decr_frequency.append((comb_name, comb_count))
            comb_most_to_least_frequent.append(
                comb_ids_decr_frequency)
        # save as dataframe
        data = {'region_upper_left_xy': region_loc, 'weights': self._t_weights, 'HTI': hti,
                'comb_most_to_least_frequent': comb_most_to_least_frequent}
        df = pd.DataFrame.from_dict(data)
        df = df.sort_values(['region_upper_left_xy'], ascending=True)
        return df


def print_results(name, hta, hta_pva_clt):
    print("\nExample {}:\n------------------\nHTA = {:.2f} (p-val = {:.2e})"
          .format(name, hta, hta_pva_clt))






