#!/usr/bin/env/ nextflow

nextflow.enable.dsl=2

params.zarrs = [
    // [['id': 'NJ_EMB-section_070_SOB26-WEM-FFPE-1'], '/home/ubuntu/Documents/spatial_annot/test_data/NJ_EMB-section_070_SOB26-WEM-FFPE-1-label.zarr', '/home/ubuntu/Documents/spatial_annot/test_data/44619_download', '/home/ubuntu/Documents/spatial_annot/test_data/NJ_EMB-section_070_SOB26-WEM-FFPE-1-anndata.zarr'],
    // [['id': 'test'], 's3://bayraktar/NJ_EMB/Visium/0.3.2/NJ_EMB-section_070_SOB26-WEM-FFPE-1-anndata.zarr', 's3://bayraktar/NJ_EMB/Visium/0.3.2/NJ_EMB-section_070_SOB26-WEM-FFPE-1-label.zarr'],

    /*[['id': 'Slide3_HandE_A1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_A1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52440_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_A1-anndata.zarr'],*/
    /*[['id': 'Slide3_HandE_B1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_B1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52441_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_B1-anndata.zarr'],*/
    /*[['id': 'Slide3_HandE_C1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_C1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52442_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_C1-anndata.zarr'],*/
    /*[['id': 'Slide3_HandE_D1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_D1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52443_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_D1-anndata.zarr'],*/

    [['id': 'Slide4_HandE_A1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_A1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52444_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide4_HandE_A1-anndata.zarr'],
    [['id': 'Slide4_HandE_B1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_B1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52445_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide4_HandE_B1-anndata.zarr'],
    [['id': 'Slide4_HandE_C1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_C1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52446_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide4_HandE_C1-anndata.zarr'],
    [['id': 'Slide4_HandE_D1'], '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide3_HandE_D1-label.zarr', '/lustre/scratch126/cellgen/team283/tl10/omero_roi_image_download/GBM/52447_download', '/lustre/scratch126/cellgen/team283/tl10/nf_runs/webatlas/JL_CVLNG/output/0.3.2/JL_CVLNG-Slide4_HandE_D1-anndata.zarr'],
]
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
    label_annot = channel.from(params.zarrs).map { meta, label_zarr, roi_folder, anndata_zarr ->
        tuple(meta, label_zarr, roi_folder)
    }
    anndata = channel.from(params.zarrs).map { meta, label_zarr, roi_folder, anndata_zarr ->
        tuple(meta, anndata_zarr)
    }
    ids_to_rois(label_annot)
    to_vitessce_json(ids_to_rois.out.join(anndata))
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
