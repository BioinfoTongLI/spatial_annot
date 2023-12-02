#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2
params.out_dir = '.output'


process ids_to_rois {
    debug true
    cache true

    container 'bioinfotongli/spatial_annot:latest'
    /*containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs --nv':'-v /lustre:/lustre -v /nfs:/nfs --gpus all'}"*/
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    // publishDir params.out_dir, mode:"copy"
    storeDir params.out_dir

    input:
    tuple val(meta), path(label_zarr), path(roi_folder)

    output:
    tuple val(meta), path(out)

    script:
    out = "${meta['id']}.json"
    def args = task.ext.args ?: ''
    """
    add_omero_annot.py \
        -label_zarr ${label_zarr} \
        -roi_folder ${roi_folder} \
        -out ${out} \
        ${args}
    """
}


process to_vitessce_json {
    debug true
    cache true

    container 'bioinfotongli/spatial_annot:latest'
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    publishDir params.out_dir, mode:"copy"

    input:
    tuple val(meta), path(annot_json), path(anndata_zarr)

    output:
    tuple val(meta), path("${out}_vitessce.json"), path("${out}_anndata_with_annotations.h5ad")

    script:
    out = meta['id']

    def args = task.ext.args ?: ''
    """
    to_Vitessce_annot.py \
        -annot_json ${annot_json} \
        -anndata_zarr ${anndata_zarr} \
        -out ${out} \
        ${args}
    """
}


process vitessce_annot_to_anndata {
    debug true
    cache true

    container 'bioinfotongli/spatial_annot:latest'
    containerOptions "${workflow.containerEngine == 'singularity' ? '-B /lustre,/nfs':'-v /lustre:/lustre -v /nfs:/nfs'}"
    publishDir params.out_dir, mode:"copy"

    input:
    tuple val(meta), path(annot_json), path(anndata_zarr)

    output:
    tuple val(meta), path("${out}_anndata_with_annotations.h5ad")

    script:
    out = meta['id']
    def args = task.ext.args ?: ''
    """
    add_vitessce_annot_to_anndata.py \
        -annot_json ${annot_json} \
        -anndata_zarr ${anndata_zarr} \
        -out ${out} \
        ${args}
    """
}


workflow {
    zarrs = channel.fromPath("/nfs/team283_imaging/LP_GBM/playground_Tong/webatlas_convert/spaceranger_ID_map - visium_sample_map.csv")
        .splitCsv(header: true, sep: ',')
        .multiMap { row ->
            label_annot: [['id': row['SANGER_ID']],
                file("/nfs/team283_imaging/LP_GBM/playground_Tong/webatlas_convert/zarrs/0.3.2/LP_GBM-" + row['SANGER_ID']+ "-label.zarr", checkIfExists:true),
                file("/nfs/team283_imaging/LP_GBM/playground_Tong/ROIs/" + row['omero_id'] + "_download", checkIfExists:true)]
            anndata: [['id': row['SANGER_ID']], 
                file("/nfs/team283_imaging/LP_GBM/playground_Tong/webatlas_convert/zarrs/0.3.2/LP_GBM-" + row['SANGER_ID']+ "-anndata.zarr", checkIfExists:true)]
        }
    ids_to_rois(zarrs.label_annot)
    to_vitessce_json(ids_to_rois.out.join(zarrs.anndata))
}

params.vitessce_annot_transfer_input = [
    [['id': 'KR_XEN_XETG00155__0004120__Region_1__20231011__125831'], '/nfs/team283_imaging/KR_XEN/playground_Tong/20231025_mouse_brain_QC/My Selections_vitessce-obs-hierarchy KR_XEN_XETG00155__0004120__Region_1__20231011__125831.json', 'xenium/public/KR_XEN/KR_XEN-output-XETG00155__0004120__Region_1__20231011__125831-anndata.zarr'],
    [['id': 'KR_XEN_XETG00155__0004120__Region_2__20231011__125831'], '/nfs/team283_imaging/KR_XEN/playground_Tong/20231025_mouse_brain_QC/My Selections_vitessce-obs-hierarchy KR_XEN_XETG00155__0004120__Region_2__20231011__125831.json', 'xenium/public/KR_XEN/KR_XEN-output-XETG00155__0004120__Region_2__20231011__125831-anndata.zarr'],
    [['id': 'KR_XEN_XETG00155__0004142__Region_1__20231011__125831'], '/nfs/team283_imaging/KR_XEN/playground_Tong/20231025_mouse_brain_QC/My Selections_vitessce-obs-hierarchy KR_XEN_XETG00155__0004142__Region_1__20231011__125831.json', 'xenium/public/KR_XEN/KR_XEN-output-XETG00155__0004142__Region_1__20231011__125831-anndata.zarr'],
    [['id': 'KR_XEN_XETG00155__0004142__Region_2__20231011__125831'], '/nfs/team283_imaging/KR_XEN/playground_Tong/20231025_mouse_brain_QC/My Selections_vitessce-obs-hierarchy KR_XEN_XETG00155__0004142__Region_2__20231011__125831.json', 'xenium/public/KR_XEN/KR_XEN-output-XETG00155__0004142__Region_2__20231011__125831-anndata.zarr'],
]

workflow add_vitessce_annotations_to_anndata {
    vitessce_annot_to_anndata(channel.from(params.vitessce_annot_transfer_input))
}
