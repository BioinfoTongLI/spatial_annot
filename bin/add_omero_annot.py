#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

"""
import fire
from typing import Dict
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader
import anndata as ad
from glob import glob
import yaml
import dask.array as da
import json

from shapely import wkt
from shapely.geometry import Polygon

from skimage.draw import polygon
try:
    import cupy as xp
    print("Using cupy")
except ImportError:
    import numpy as xp
    print("Using numpy")

# import scanpy as sc
# from scipy.sparse import csr_matrix
# import json


def read_and_qc(sample_name:str,
                count_file='filtered_feature_bc_matrix.h5'):
    r""" This function reads the data for one 10X spatial experiment into the anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.

    :param sample_name: Name of the sample
    :param count_file: path to matrix
    """
    adata = sc.read_visium(str(sample_name),
                           count_file=count_file)
    
    fiducial_file = glob(f"{sample_name}/*_fiducial.json")[0]
    with open(fiducial_file) as f:
        # Load the JSON data
        align_data = json.load(f)
    adata.obs['sample'] = list(adata.uns['spatial'])[0]
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)
    # fix TypeError when read in obsm
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)
    # Calculate QC metrics
    adata.X = adata.X.toarray()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.X = csr_matrix(adata.X)
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    
    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    adata.uns["cytassist_align"] = align_data
    return adata 


# def main(omero_annot_yaml:str, spaceranger_folder:str):
    # adata = read_and_qc(spaceranger_folder)
    # print(adata.uns['spatial'][list(adata.uns['spatial'])[0]].keys())
    # print(adata.uns['spatial'][list(adata.uns['spatial'])[0]]['scalefactors'])
    # print(adata.uns["cytassist_align"]['transform'])
    # regions = {}
    # for y in glob(f"{omero_annot_yaml}/*.yaml"):
    #     with open(y, 'r') as f:
    #         data = yaml.safe_load(f)
    #     regions[data['roi_name']] = wkt.loads(data['poly'])
    # print(regions)
    # adata.write_h5ad(f"{spaceranger_folder}.h5ad")


def load_annotations_from_yaml(roi_folder: str) -> Dict[str, Polygon]:
    """
    Load regions from YAML files in a folder.

    :param roi_folder: Path to the folder containing the YAML files.
    :return: A dictionary with the regions loaded from the YAML files.
    """
    regions = {}
    for y in glob(f"{roi_folder}/*.yaml"):
        with open(y, 'r') as f:
            data = yaml.safe_load(f)
        regions[data['roi_name']] = Polygon(wkt.loads(data['poly']))
    return regions


def read_label_zarr(label_zarr_url: str) -> da.array:
    """
    Reads the label zarr file from the given URL using ome-zarr and returns the label image as a numpy array.

    :param label_zarr_url: The URL of the label zarr file.
    :type label_zarr_url: str

    :return: The label image as a dask array.
    :rtype: dask.array
    """
    # Read in the label zarr using ome-zarr
    reader = Reader(parse_url(label_zarr_url))
    lab_img = list(reader())[0].data[0]

    return lab_img


# Assign manual annotations to Visium spots
def main(label_zarr:str, roi_folder:str, out:str):
    rois = load_annotations_from_yaml(roi_folder)
    print(rois)

    # read in the label zarr using ome-zarr
    lab_img = read_label_zarr(f"{label_zarr}/0").squeeze()
    print(lab_img.shape)

    visium_spots_in_rois = {}
    for r in rois:
        # Get the coordiantes of pixels inside the polygon
        rr, cc = polygon(rois[r].exterior.coords.xy[1], rois[r].exterior.coords.xy[0], shape=lab_img.shape)
        try:
            visium_spots_in_rois[r] = [int(i) for i in xp.unique(xp.array(lab_img)[rr, cc]).get()]
        except:
            visium_spots_in_rois[r] = [int(i) for i in xp.unique(xp.array(lab_img)[rr, cc])]
    print(visium_spots_in_rois)

    with open(out, 'w') as f:
        json.dump(visium_spots_in_rois, f)

    
if __name__ == "__main__":
    fire.Fire(main)