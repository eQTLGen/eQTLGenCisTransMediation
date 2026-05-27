#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

cis-trans mediation analysis v${workflow.manifest.version}"
===========================================================
Pipeline for running cis-eQTL - trans-eQTL mediation analysis for each fine-mapped trans-eQTL

Usage:

nextflow run main.nf 
--finemapped \
--vcf \
--expression_matrix \
--covariate_matrix \
--covariates \
--OutputDir

Mandatory arguments:
--finemapped                eQTLGen file with fine-mapped eQTL results.
--covariates                Comma-separated string with covariates to use.
--vcf                       Folder with imputed vcf files (Imputation pipeline output).
--expression_matrix         Preprocessed gene expression matrix (DataQC pipeline output).
--covariate_matrix          Covariate file  (DataQC pipeline output).

Optional arguments:
--OutputDir                 Output directory. Defaults to "results".
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
// params.covariates = 'GenPC1,GenPC2,GenPC3,GenPC4,ExpPC1,ExpPC2,ExpPC3,ExpPC4,ExpPC5,ExpPC6,ExpPC7,ExpPC8,ExpPC9,ExpPC10,ExpPC11,ExpPC12,ExpPC13,ExpPC14,ExpPC15,ExpPC16,ExpPC17,ExpPC18,ExpPC19,ExpPC20'
params.expression_matrix = ''

//Show parameter values
log.info """================================================================
eQTLGen cis-trans mediation analysis pipeline v${workflow.manifest.version}"
==========================================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['Finemapped results']                       = params.finemapped
summary['Genotype vcf folder']                      = params.vcf
summary['Expression matrix']                        = params.expression_matrix
summary['Covariate matrix']                         = params.covariate_matrix
summary['Covariates']                               = params.covariates

// import modules
include { PREPAREINPUTS; MEDIATIONANALYSIS; MediationAnalysis; PrepareInputs } from './modules/MediationAnalysis.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
finemapped_ch = Channel.fromPath(params.finemapped)

genotype_ch = Channel
    .fromPath("${params.vcf}/chr*.filtered.vcf.gz")
    .map { file ->
        def chr = (file.baseName =~ /chr(\d+)\.filtered\.vcf/)[0][1]
        tuple(chr, file)
    }

expression_matrix_ch = Channel.fromPath(params.expression_matrix)
covariate_matrix_ch = Channel.fromPath(params.covariate_matrix)

// Define parameter channels
// Currently no parameter channels
covariate_ch = Channel.value(params.covariates)

workflow {

        PREPAREINPUTS(finemapped_ch)

        triplets_ch = PREPAREINPUTS
        .out
        .flatten()
        .map { file ->
        def chr = (file.baseName =~ /^(\d+)_batch\d+$/)[0][1]
        tuple(chr, file)
        }

        genotype_ch.view()

        triplets_ch = triplets_ch.combine(genotype_ch, by: 0)

        combined_ch = triplets_ch
            .combine(covariate_ch)
            .combine(expression_matrix_ch)
            .combine(covariate_matrix_ch)

        MEDIATIONANALYSIS(combined_ch)

        MEDIATIONANALYSIS.out.collectFile(name: 'CisTransMediationResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")

}
workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
