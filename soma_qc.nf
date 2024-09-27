#!/Users/tie_zhao/miniconda3/envs/dsl1 nextflow


params.naga_dir = '/Users/tie_zhao/Desktop/Soma_QC_re/DataSrc/NAGA'
params.phom_dir = '/Users/tie_zhao/Desktop/Soma_QC_re/DataSrc/PHOM'
params.result_dir = '/Users/tie_zhao/Desktop/Soma_QC_re/Results'
params.script = '/Users/tie_zhao/Desktop/Soma_QC_re/Scripts'
params.adat_summary = '/Users/tie_zhao/Desktop/Soma_QC_re/DataSrc/adat_summary.xlsx'

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

naga_qc_out_ch
    .map { name, path -> path }
    .reduce { a, b -> [a, b] }
    .map { paths -> ["NAGA"].plus(paths) }
    .set { naga_merge_ch }

process merge_naga_qc {

    tag "${group_name}"

    conda '/Users/tie_zhao/miniconda3'
    publishDir "${params.result_dir}/02.merged_naga", mode: 'symlink'

    input:
    tuple val(group_name), path(naga_qc_rds_1), path(naga_qc_rds_2) from naga_merge_ch

    output:
    tuple val(group_name), file(merged_rds) into naga_out_ch

    script:
    merged_rds = "${group_name}.qc.merge.rds"
    """
    Rscript ${params.script}/merge_rds.r -r ${naga_qc_rds_1} -s ${naga_qc_rds_2} -o ${group_name}
    """
}

phom_qc_out_ch
    .filter { it[0].contains('SS-229363') }
    .map { item -> ["PHOM", item[1]]}
    .set {phom_day0_in_ch}


process phom_day0_selection {

    tag "${group_name}:Day0"

    conda '/Users/tie_zhao/miniconda3'
    publishDir "${params.result_dir}/03.phom_day0", mode: 'symlink'

    input:
    tuple val(group_name), path(rds_file) from phom_day0_in_ch

    output:
    tuple val(group_name), file(phom_day0_rds) into phom_day0_out_ch, phom_day0_uniprot_ch
    file(phom_ex_csv)
    file(phom_apt_csv)

    script:
    phom_day0_rds = "${group_name}.qc.day0.rds"
    phom_ex_csv = "${group_name}.ex.day0.csv"
    phom_apt_csv = "${group_name}.apt.day0.csv"
    """
    Rscript ${params.script}/rds_day0.r -r ${rds_file} -o ${group_name}

    Rscript -e '
    data <- readRDS("${phom_day0_rds}")
    data[[1]] <- data[[1]][!grepl("QC", rownames(data[[1]])), ]
    write.csv(data[[1]], "${group_name}.ex.day0.csv", row.names = TRUE, col.names = TRUE)
    write.csv(data[[2]], "${group_name}.apt.day0.csv", row.names = FALSE, col.names = TRUE)
    '
    """
}

naga_out_ch
    .combine(phom_day0_out_ch)
    .map { item -> ['PHOMday0_vs_NAGA', item[3], item[1]] }
    .set {phom_day0_vs_naga_in_ch}

process phom_day0_vs_naga {
    
    tag "${group_name}"
    
    conda '/Users/tie_zhao/miniconda3'
    publishDir "${params.result_dir}/04.phom_day0_vs_naga", mode: 'symlink'
    
    input:
    tuple val(group_name), path(phom_day0_rds), path(naga_rds) from phom_day0_vs_naga_in_ch
    path(adat_summary) from params.adat_summary

    output:
    tuple val(group_name), file(phom_day0_vs_naga_rds) into phom_day0_vs_naga_out_ch
    file(phom_day0_vs_naga_ex_csv)
    file(phom_day0_vs_naga_apt_csv)
    file(phom_day0_vs_naga_qc_plot)
    
    script:
    phom_day0_vs_naga_rds = "${group_name}.qc.merge.rds"
    phom_day0_vs_naga_ex_csv = "${group_name}.ex.csv"
    phom_day0_vs_naga_apt_csv = "${group_name}.apt.csv"
    phom_day0_vs_naga_qc_plot = "${group_name}.qc.pca_plot.html"
    """
    Rscript ${params.script}/merge_rds.r -r ${phom_day0_rds} -s ${naga_rds} -o ${group_name}

    Rscript -e '
    data <- readRDS("${phom_day0_vs_naga_rds}")
    data[[1]] <- data[[1]][!grepl("QC", rownames(data[[1]])), ]
    write.csv(data[[1]], "${group_name}.ex.csv", row.names = TRUE, col.names = TRUE)
    write.csv(data[[2]], "${group_name}.apt.csv", row.names = FALSE, col.names = TRUE)
    '

    Rscript ${params.script}/qc_pca_plot.r -r ${phom_day0_vs_naga_rds} -s ${adat_summary} -o ${group_name} 
    """
}


phom_day0_uniprot_ch
    .map { item -> ['PHOMday0', item[1]] }
    .mix(phom_day0_vs_naga_out_ch)
    .set { seqid2uniprot_in_ch }


process seqid2uniprot_in_ch {

    tag "${group_name}"

    conda '/Users/tie_zhao/miniconda3'
    publishDir "${params.result_dir}/05.seqid2uniprot", mode: 'symlink'

    input:
    tuple val(group_name), path(seqid_rds) from seqid2uniprot_in_ch

    output:
    tuple val(group_name), file(uniprot_rds) into seqid2uniprot_out_ch
    file(ex_uniprot_csv)
    file(apt_uniprot_csv)
    file(summary_uniprot_csv)

    script:
    uniprot_rds = "${group_name}.uniprot.rds"
    ex_uniprot_csv = "${group_name}.ex.uniprot.csv"
    apt_uniprot_csv = "${group_name}.apt.uniprot.csv"
    summary_uniprot_csv = "${group_name}.summary.uniprot.csv"
    """
    Rscript ${params.script}/seqid2uniprot.r -r ${seqid_rds} -o ${group_name}
    Rscript -e '
    data <-readRDS("${group_name}.uniprot.rds")
    data[[1]] <- data[[1]][!grepl("QC", rownames(data[[1]])), ]
    write.csv(data[[1]], "${group_name}.ex.uniprot.csv", row.names = TRUE, col.names = TRUE)
    write.csv(data[[2]], "${group_name}.apt.uniprot.csv", row.names = FALSE, col.names = TRUE)
    write.csv(data[[3]], "${group_name}.summary.uniprot.csv", row.names = FALSE, col.names = TRUE)
    '
    """
}

