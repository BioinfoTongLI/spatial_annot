#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

"""
import fire
import json
import anndata as ad


def main(annot_json:str, anndata_zarr:str, out:str):
    with open(annot_json, 'r') as f:
        config = json.load(f)
    for child in config["tree"][0]['children']:
        selected_ids = [i[0] for i in child['set']]

    adata = ad.read_zarr(anndata_zarr)
    adata.obs["Annotation"] = "NA"
    adata.obs.loc[selected_ids, "Annotation"] = config["tree"][0]['children'][0]['name']
    zarr_out = f"{out}_anndata_with_annotations.h5ad"
    adata.write_h5ad(zarr_out)


if __name__ == '__main__':
    fire.Fire(main)
