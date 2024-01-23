#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

"""
import fire
from typing import Dict
import json
import anndata as ad
import pandas as pd


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


def main(annot_csv:str, anndata_zarr:str, out:str):
    annot_count_df = pd.read_csv(annot_csv, index_col=0).astype(int)
    annot_count_df.index = annot_count_df.index.astype(str)

    annot_cols = annot_count_df.columns
    adata = ad.read_zarr(anndata_zarr)
    merged_df = adata.obs.merge(
        annot_count_df, left_index=True, right_index=True, how='outer'
    )
    for col in annot_cols:
        merged_df[col].fillna(0, inplace=True)

    print(merged_df)
    zarr_out = f"{out}_anndata_with_annotations.h5ad"
    adata.obs = merged_df  # Update adata.obs to merged_df

    adata.write_h5ad(zarr_out)

    # children = {i: construct_set(k, rois_with_ids[k]) for i, k in enumerate(rois_with_ids)}

    # config = VITESSCE_CONFIG_TEMPLATE.copy()
    # config['tree'][0]['children'] = children
    # # print(config['tree'][0]['children'])
    # with open(f"{out}_vitessce.json", 'w') as f:
    #     json.dump(config, f)


if __name__ == '__main__':
    fire.Fire(main)
