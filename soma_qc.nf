#!/Users/tie_zhao/miniconda3/envs/dsl1 nextflow


params.naga_dir = '/Users/tie_zhao/Desktop/Soma_QC_re/DataSrc/NAGA'
params.phom_dir = '/Users/tie_zhao/Desktop/Soma_QC_re/DataSrc/PHOM'
params.result_dir = '/Users/tie_zhao/Desktop/Soma_QC_re/Results'
params.script = '/Users/tie_zhao/Desktop/Soma_QC_re/Scripts'

Channel
    .fromFilePairs("${params.naga_dir}/*.{adat,sample}", size: 2, flat: true)
    .set { naga_adat_ch }

Channel
    .fromFilePairs("${params.phom_dir}/*.{adat,sample}", size: 2, flat: true)
    .set { phom_adat_ch }


process naga_adat {

    tag "NAGA: ${adat_name}"

    conda '/Users/tie_zhao/miniconda3'
    publishDir "${params.result_dir}/01.aptamer_sample_qc/NAGA", mode: 'symlink'

    input:
    tuple val(adat_name), path(adat_file), path(sample_file) from naga_adat_ch

    output:
    tuple val(adat_name), file(qc_out_rds) into naga_qc_out_ch

    script:
    qc_out_rds = "${adat_name}.qc.rds"
    """
    Rscript ${params.script}/adat_qc.r -d ${adat_file} -s ${sample_file} -o ${adat_name}
    """
}

process phom_adat {

    tag "PHOM: ${adat_name}"

    conda '/Users/tie_zhao/miniconda3'
    publishDir "${params.result_dir}/01.aptamer_sample_qc/PHOM", mode: 'symlink'

    input:
    tuple val(adat_name), path(adat_file), path(sample_file) from phom_adat_ch

    output:
    tuple val(adat_name), file(qc_out_rds) into phom_qc_out_ch

    script:
    qc_out_rds = "${adat_name}.qc.rds"
    """
    Rscript ${params.script}/adat_qc.r -d ${adat_file} -s ${sample_file} -o ${adat_name}
    """
}
