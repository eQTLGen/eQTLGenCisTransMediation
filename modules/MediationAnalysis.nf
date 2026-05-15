#!/bin/bash nextflow

process PREPAREINPUTS {

    scratch '$TMPDIR'

    input:
        path(finemapped_results)

    output:
        path("*_batch*")

    shell:
        """
        PrepareInputs.R \
        --finemapped !{finemapped_results}
        """

}


process MEDIATIONANALYSIS {

    scratch '$TMPDIR'

    input:
        tuple val(chr), path(triplets), path(vcf_file), val(covariates), path(exp_file), path(cov_file)

    output:
        path("*_SobelResults.txt")

    shell:
        """
        MediationAnalysis.R \
        --triplets !{triplets} \
        --cov !{covariates} \
        --vcf_file !{vcf_file} \
        --exp_file !{exp_file} \
        --cov_file !{cov_file}
        """

}

workflow MediationAnalysis {
    take:
        data

    main:
        mediationanalysis_output_ch = MediationAnalysis(data)

    emit:
        mediationanalysis_output_ch
}

workflow PrepareInputs {
    take:
        data

    main:
        prepareinputs_output_ch = PrepareInputs(data)

    emit:
        prepareinputs_output_ch
}
