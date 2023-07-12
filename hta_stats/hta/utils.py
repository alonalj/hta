import os

import numpy as np
import csv
import gzip
import scipy.io
import pandas as pd
import matplotlib.pyplot as plt
import math


def apply_tissue_mask(trait_tensor, tissue_mask):
    '''
    trait_tensor: the trait-tensor -- a numpy ndarray where the first two or three dimensions
    represent the spatial dimensions (i.e., the x, y, z axes). The remaining dimensions, each, represent a single trait,
    with binary values that depict whether that trait manifests in that x,y,z position (there could be more if e.g., a
    time dimension is added, or similar, but we did not demonstrate this in the paper).
    For example, in a 2d space, a trait_tensor of shape (16,16,4) describes a 16x16 space (x,y) and in each such
    (x,y) position, there are 4 traits. The element (x,y,0), would have 1 or 0, depending on whether the first trait
    manifests or not, respectively.

    tissue_mask: binary numpy 2d array, of shape (trait_tensor.shape[0], trait_tensor.shape[1]) where 1 represents 'keep'
    and 0 represents 'disregard'.

    Multiplies the trait tensor element-wise by the binary tissue_mask, of shape (trait_tensor.shape[0], trait_tensor.shape[1])
    to enable disregarding positions where measurements are irrelevant (e.g. noise due to barcoded spot but not tissue).
    '''
    tissue_mask_tensor = [np.expand_dims(tissue_mask, -1) for i in range(trait_tensor.shape[-1])]
    tissue_mask_tensor = np.concatenate(tissue_mask_tensor, -1)
    t = np.multiply(trait_tensor, tissue_mask_tensor)
    return t


def apply_empty_mask(t):
    '''
    Removes entire rows/columns where there is zero across all traits for that row/column

    t: the trait-tensor -- a numpy ndarray where the first two or three dimensions
    represent the spatial dimensions (i.e., the x, y, z axes). The remaining dimensions, each, represent a single trait,
    with binary values that depict whether that trait manifests in that x,y,z position (there could be more if e.g., a
    time dimension is added, or similar, but we did not demonstrate this in the paper).
    For example, in a 2d space, a trait_tensor of shape (16,16,4) describes a 16x16 space (x,y) and in each such
    (x,y) position, there are 4 traits. The element (x,y,0), would have 1 or 0, depending on whether the first trait
    manifests or not, respectively.
    '''
    if len(t.shape) == 2:
        t = t[~np.all(t == 0, axis=1)]
    else:
        mask = ~(t == 0).all(axis=(1, 2))
        t = t[mask]
        mask = ~(t == 0).all(axis=(0, 2))
        t = t[:, mask, :]
    return t


def round_tensor_shape(t):
    '''
    Rounds the shape of the tensor to the nearest 10 by padding with zero. Can be used to make region sizes visually even
    (the number of elements in each region might still differ among regions, but when plotting, the grid would be
    evenly dispersed across the sample).

    t: the trait-tensor -- a numpy ndarray where the first two or three dimensions
    represent the spatial dimensions (i.e., the x, y, z axes). The remaining dimensions, each, represent a single trait,
    with binary values that depict whether that trait manifests in that x,y,z position (there could be more if e.g., a
    time dimension is added, or similar, but we did not demonstrate this in the paper).
    For example, in a 2d space, a trait_tensor of shape (16,16,4) describes a 16x16 space (x,y) and in each such
    (x,y) position, there are 4 traits. The element (x,y,0), would have 1 or 0, depending on whether the first trait
    manifests or not, respectively.
    '''
    def roundup(x):
        return int(math.ceil(x / 10.0)) * 10

    if t.shape[0] != t.shape[1] or t.shape[0] % 10 != 0:
        max_shape = max(t.shape[0], t.shape[1])
        max_shape = roundup(max_shape)
        new_shape = list(t.shape)
        new_shape[0], new_shape[1] = max_shape, max_shape
        new_t = np.zeros(new_shape)
        for i in range(t.shape[-1]):
            new_t[:t.shape[0], :t.shape[1], ...] = t
        t = new_t
    return t


def plot_heatmap(matrix, name, out_path='../out/'):
    matrix = apply_empty_mask(matrix)
    plt.imshow(matrix, cmap='magma',
               interpolation='nearest')  # TODO: remove interpolation? # viridis
    plt.colorbar()
    plt.savefig(out_path+'{}.jpeg'.format(name))  # TODO add p_trait to name to make sure it computes it correctly
    plt.close()


class Visium():

    def __init__(self, path):
        self.path = path
        self.feature_matrix, self.traits_available, self.barcode_ids = None, None, None
        self.barcode_positions = None
        self.barcode_to_row_col = {}
        self.tissue_positions = None
        self.matrix_shape = None
        self.trait_tensor = None
        self.tissue_mask = None
        self.shape = None

    def load(self):
        '''
        Loads the feature matrix into self.feature_matrix, barcode ids into self.barcode_ids, traits (e.g. genese)
        available into self.traits_available and barcode positions (row, col in tissue) into the dictionary
        self.barcode_to_row_col
        '''
        def load_feature_matrix():
            path = self.path + '/filtered_feature_bc_matrix/'
            mat = scipy.io.mmread(os.path.join(path, "matrix.mtx.gz"))

            # feature and barcode sequences corresponding to row and column indices respectively
            features_path = os.path.join(path, "features.tsv.gz")
            feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]  # rows
            genes_available = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
            feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
            barcodes_path = os.path.join(path, "barcodes.tsv.gz")
            barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]  # cols

            m = mat.tocsr()  # genes x barcodes (rows x cols)
            # b = m.toarray()
            # b = m[:1000,:1000].toarray()
            self.feature_matrix = m
            self.traits_available = genes_available
            self.barcode_ids = barcodes

        def load_barcode_positions_with_tissue():
            path = self.path + '/spatial/'
            path = os.path.join(path, "tissue_positions_list.csv")
            cols = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
            barcode_positions = pd.read_csv(path, names=cols)
            self.barcode_positions = barcode_positions # TODO: check if this is used anywhere and if not remove self.barcode_positions so doesn't take mem
            self.barcode_positions = self.barcode_positions.loc[self.barcode_positions['in_tissue'] == 1, :]
            assert len(self.barcode_positions['barcode']) == len(self.barcode_ids)

            for b in self.barcode_ids:
                d_barcode_positions = self.barcode_positions[self.barcode_positions['barcode'] == b]
                row, col = d_barcode_positions['array_row'].values[0], d_barcode_positions['array_col'].values[0]
                position = (row, col)
                self.barcode_to_row_col[b] = position

        load_feature_matrix()
        load_barcode_positions_with_tissue()

    def _set_shape_where_tissue(self):
        x_max, y_max = self.barcode_positions['array_row'].max(), self.barcode_positions['array_col'].max()
        shape = max(x_max, y_max)
        matrix_shape = (shape + 1, shape + 1)
        self.matrix_shape = matrix_shape

    def _make_tissue_mask(self):
        tissue_mask = np.zeros(self.matrix_shape)

        tissue_positions = []
        for b in self.barcode_ids:
            position = self.barcode_to_row_col[b]
            row, col = position
            tissue_positions.append(position)
            tissue_mask[row, col] = 1

        # plot_heatmap(tissue_mask, 'tissue_mask') # TODO
        self.tissue_positions = tissue_positions
        self.tissue_mask = tissue_mask

    def _create_matrix_for_trait(self, gene_row_ix):
        '''
        :param m: matrix with trait in rows and barcodes in cols
        :param gene_row_ix:
        :param positions: positions list for each barcode. The j'th entry (position tuple) is for the j'th column (barcode) in m
        :return:
        '''
        assert self.feature_matrix.shape[1] == len(self.tissue_positions)
        matrix_trait = np.zeros(self.matrix_shape)
        for j in range(len(self.tissue_positions)):
            row, col = self.tissue_positions[j]
            matrix_trait[row, col] = self.feature_matrix[gene_row_ix, j]

        # matrix_trait_tight = matrix_trait[~np.all(matrix_trait == 0, axis=1)]
        # matrix_trait_tight = matrix_trait_tight[:, ~np.all(matrix_trait_tight == 0, axis=0)]

        return matrix_trait

    def _make_trait_tensor(self, trait_names):
        matrices = []
        self._assert_trait_names_in_data(trait_names)
        for trait in trait_names:
            trait_row_ix = self.traits_available.index(trait)
            trait_matrix = self._create_matrix_for_trait(trait_row_ix)
            trait_matrix = np.expand_dims(trait_matrix, axis=-1)
            matrices.append(trait_matrix)
        trait_tensor = np.concatenate(matrices, axis=-1)
        trait_tensor = apply_tissue_mask(trait_tensor, self.tissue_mask)
        self.trait_tensor = trait_tensor

    def _make_cluster_tensor(self, clusters_filepath):
        cluster_df = pd.read_csv(clusters_filepath)
        clusters_sorted = sorted(list(set(cluster_df.iloc[:, 1])))  # assumes clusters are in second col and are numbers from 1,2,3...
        n_clusters = len(clusters_sorted)
        trait_tensor = np.zeros((self.matrix_shape[0], self.matrix_shape[1], n_clusters))
        for i in range(len(cluster_df)):
            row, col = self.barcode_to_row_col[cluster_df.iloc[i, 0]]
            cluster_id = cluster_df.iloc[i, 1]-1
            trait_tensor[row, col, cluster_id] = 1
        trait_tensor = apply_tissue_mask(trait_tensor, self.tissue_mask)
        self.trait_tensor = trait_tensor

    def _make_curated_tensor(self, df, threshold_fn=np.median):
        barcodes = df.columns
        df = df.values
        n_traits = df.shape[0]
        trait_tensor = np.zeros((self.matrix_shape[0], self.matrix_shape[1], n_traits))
        for trait_ix in range(n_traits):
            trait_threshold_value = threshold_fn(df[trait_ix, :])
            for barcode_ix in range(len(barcodes)):
                row, col = self.barcode_to_row_col[barcodes[barcode_ix]]
                value = (df[trait_ix, barcode_ix] >= trait_threshold_value)*1
                trait_tensor[row, col, trait_ix] = value
        trait_tensor = apply_tissue_mask(trait_tensor, self.tissue_mask)
        self.trait_tensor = trait_tensor

    def _preprocess_tensor_and_tissue_mask(self):
        trait_tensor_mask_aligned = [self.trait_tensor[:, :, i] for i in range(self.trait_tensor.shape[-1])]
        trait_tensor_mask_aligned.append(self.tissue_mask)
        trait_tensor_mask_aligned = [np.expand_dims(i, -1) for i in trait_tensor_mask_aligned]
        trait_tensor_mask_aligned = np.concatenate(trait_tensor_mask_aligned, -1)
        trait_tensor_mask_aligned = apply_empty_mask(trait_tensor_mask_aligned)
        self._raw_tensor_mask_aligned = trait_tensor_mask_aligned
        trait_tensor_mask_aligned = round_tensor_shape(trait_tensor_mask_aligned)

        trait_tensor, tissue_mask = trait_tensor_mask_aligned[:, :, :-1], trait_tensor_mask_aligned[:, :, -1]
        self.trait_tensor = trait_tensor
        self.tissue_mask = tissue_mask

    def _assert_trait_names_in_data(self, trait_names):
        for t in trait_names:
            assert t in self.traits_available, "Trait {} not in data. Use .traits_available to see which are available.".format(t)

    # def plot_heatmap_per_trait(self, out_path):
    #     for trait_ix in range(self.trait_tensor.shape[-1]):
    #         plot_heatmap(self.trait_tensor[:, :, trait_ix], 'visium_{}'.format(self.traits_available[trait_ix]), out_path)

    def _binarize(self, trait_names, threshold_fn=np.median): #TODO: separate binarization from submatrix with trait names
        for trait_ix in range(len(trait_names)):
            trait_row_ix = self.traits_available.index(trait_names[trait_ix])
            trait_values = self.feature_matrix[trait_row_ix, :].data
            trait_threshold_value = threshold_fn(trait_values)
            self.trait_tensor[:, :, trait_ix] = self.trait_tensor[:, :, trait_ix] >= trait_threshold_value

    def prep(self, trait_names=None, threshold_fn=np.median):
        self._assert_trait_names_in_data(trait_names)
        self._set_shape_where_tissue()
        self._make_tissue_mask()
        self._make_trait_tensor(trait_names)
        self._preprocess_tensor_and_tissue_mask()
        self._binarize(trait_names, threshold_fn)

        return self.trait_tensor, self.tissue_mask

    def prep_clusters(self, clusters_filepath):
        self._set_shape_where_tissue()
        self._make_tissue_mask()
        self._make_cluster_tensor(clusters_filepath)
        self._preprocess_tensor_and_tissue_mask()
        return self.trait_tensor, self.tissue_mask

    def prep_curated_traits(self, df, threshold_fn=np.median):
        self._set_shape_where_tissue()
        self._make_tissue_mask()
        self._make_curated_tensor(df, threshold_fn)
        self._preprocess_tensor_and_tissue_mask()
        return self.trait_tensor, self.tissue_mask


class Images():
    def __init__(self, path, trait_names, slice_names, img_format):
        self._path = path
        self._traits = trait_names
        self._slices = slice_names
        self._img_format = img_format
        self.trait_tensor = None

    def prep(self):
        import PIL

        def is_hot(a):
            ''' If the image is grayscale, it checks if it's above the middle threshold (more white)'''
            in_upper_hot_quadrant = (128 < a[:,:,0])
            return in_upper_hot_quadrant

        t_3d_for_trait_all = []
        for slice_name in self._slices:
            trait_images_per_slice = []
            for trait in self._traits:
                path_trait = self._path + '/{}'.format(trait)
                path_slice = path_trait + '/{}.{}'.format(slice_name, self._img_format)
                im = PIL.Image.open(path_slice)

                im = im.convert('RGB')
                im_arr = np.asarray(im)
                im_arr = is_hot(im_arr) * 1

                plot_heatmap(im_arr, 'images_bin_intensity_{}_{}'.format(trait,path_slice.split('/')[-1]))
                im_arr = np.expand_dims(im_arr, -1)
                im_arr = np.expand_dims(im_arr, -1)
                trait_images_per_slice.append(im_arr)
                im.close()
            t_3d_for_trait = np.concatenate(trait_images_per_slice, axis=-1)
            t_3d_for_trait_all.append(t_3d_for_trait)
        self.trait_tensor = np.concatenate(t_3d_for_trait_all, axis=2)

        return self.trait_tensor
