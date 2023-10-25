#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
"""

"""
import fire
import json
import anndata as ad
import s3fs


def read_remote_zarr(s3_url: str) -> ad.AnnData:
    # pip install s3fs zarr anndata

    # horribleness for Sanger
    s3 = s3fs.S3FileSystem(anon=True, client_kwargs = {
        'endpoint_url':"https://cog.sanger.ac.uk",
        'region_name':'us-east-1'
    })

    # path to anndata.zarr
    # s3_store = 'xenium/public/KR_XEN/KR_XEN-output-XETG00155__0004120__Region_2__20231011__125831-anndata.zarr'

    # create zarr store
    store = s3fs.S3Map(root=s3_url,s3=s3)
    return ad.read_zarr(store)


def main(annot_json:str, anndata_zarr:str, out:str):
    with open(annot_json, 'r') as f:
        config = json.load(f)
    for child in config["tree"][0]['children']:
        selected_ids = [i[0] for i in child['set']]

    try:
        adata = ad.read_zarr(anndata_zarr)
    except:
        adata = read_remote_zarr(anndata_zarr)
    adata.obs["Annotation"] = "NA"
    adata.obs.loc[selected_ids, "Annotation"] = config["tree"][0]['children'][0]['name']
    zarr_out = f"{out}_anndata_with_annotations.h5ad"
    adata.write_h5ad(zarr_out)


if __name__ == '__main__':
    fire.Fire(main)
