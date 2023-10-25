#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

"""
import fire
from typing import Dict
import json
import anndata as ad


VITESSCE_CONFIG_TEMPLATE = {
    "datatype": "obs",
    "version": "0.1.3",
    "tree": {
        0: {
            "name": "User Annotations",
            "children": {},
            "color":{
                0: 50,
                1: 50,
                2: 50
            }
        }
    }
}


def construct_set(object_name:str, object_ids:str):
    set = {i: {0:str(n), 1: None} for i, n in enumerate(object_ids) if n != 0}
    return {
        "name": object_name,
        "color": {
            0: 50,
            1: 50,
            2: 50
        },
        "set": set
    }


def main(annot_json:str, anndata_zarr:str, out:str):
    with open(annot_json, 'r') as f:
        rois_with_ids = json.load(f)
    children = {i: construct_set(k, rois_with_ids[k]) for i, k in enumerate(rois_with_ids)}

    config = VITESSCE_CONFIG_TEMPLATE.copy()
    config['tree'][0]['children'] = children
    # print(config['tree'][0]['children'])
    with open(f"{out}_vitessce.json", 'w') as f:
        json.dump(config, f)

    adata = ad.read_zarr(anndata_zarr)
    adata.obs["Annotation"] = "NA"

    reversed_dict = {}
    for key, values in rois_with_ids.items():
        for value in values:
            if value not in reversed_dict:
                reversed_dict[value] = key
            else:
                reversed_dict[value] += f" - {key}"
    del reversed_dict[0] # remove background label

    for t in reversed_dict:
        adata.obs["Annotation"].loc[str(t)] = reversed_dict[t]
    zarr_out = f"{out}_anndata_with_annotations.h5ad"
    adata.write_h5ad(zarr_out)


if __name__ == '__main__':
    fire.Fire(main)
