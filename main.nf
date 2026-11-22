nextflow.enable.dsl=2

process FASTQC {

    tag "$sample"

    container 'staphb/fastqc:0.12.1'

    publishDir 'results/fastqc', mode: 'copy'

    cpus 2
    memory '2 GB'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_fastqc.html"), path("*_fastqc.zip")

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        --outdir . \
        ${reads[0]} \
        ${reads[1]}
    """
}

process MULTIQC {

    tag 'Aggregate FastQC reports'

    container 'multiqc/multiqc:v1.35'

    publishDir 'results/multiqc', mode: 'copy'

    cpus 1
    memory '2 GB'

    input:
    path(fastqc_zips)

    output:
    path 'multiqc_report.html'
    path 'multiqc_report_data'

    script:
    """
    multiqc . \
        --outdir . \
        --filename multiqc_report.html \
        --force
    """
}

process SALMON_QUANT {

    tag "$sample"

    container 'combinelab/salmon:1.10.3'

    publishDir 'results/salmon/quant', mode: 'copy'

    cpus 2
    memory '3 GB'

    input:
    tuple val(sample), path(reads)
    path index, stageAs: 'salmon_index'

    output:
    path "${sample}"

    script:
    """
    salmon quant \
        -i salmon_index \
        -l A \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        --validateMappings \
        --gcBias \
        --seqBias \
        --threads ${task.cpus} \
        -o ${sample}
    """
}

process DESEQ2 {

    tag 'Gene-level differential expression'

    container 'rnaseq-deseq2:1.0'

    publishDir 'results', mode: 'copy'

    cpus 2
    memory '4 GB'

    input:
    path quant_dirs
    path samplesheet
    path tx2gene
    path analysis_script

    output:
    path 'deseq2'

    script:
    """
    Rscript "${analysis_script}" \
        "${samplesheet}" \
        "${tx2gene}" \
        deseq2
    """
}

workflow {

    samples = Channel
        .fromPath('data/metadata/samplesheet.csv')
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.run_accession,
                [
                    file(row.fastq_1),
                    file(row.fastq_2)
                ]
            )
        }

    fastqc_results = FASTQC(samples)

    fastqc_zips = fastqc_results
        .map { sample, html, zip -> zip }
        .collect()

    MULTIQC(fastqc_zips)

    salmon_index = Channel.value(
        file('results/salmon/flybase_r6.69_index')
    )

    SALMON_QUANT(samples, salmon_index)

    salmon_quant_dirs = SALMON_QUANT.out.collect()

    samplesheet_ch = Channel.value(
        file('data/metadata/samplesheet.csv')
    )

    tx2gene_ch = Channel.value(
        file('data/reference/flybase_r6.69/tx2gene.tsv')
    )

    analysis_script_ch = Channel.value(
        file('scripts/deseq2.R')
    )

    DESEQ2(
        salmon_quant_dirs,
        samplesheet_ch,
        tx2gene_ch,
        analysis_script_ch
    )
}
