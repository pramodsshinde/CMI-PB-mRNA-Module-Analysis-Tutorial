##### Global variables
ANTIBODY_GEN_IDS = ['IgG-FHA', 'IgG-PT','IgG-PRN']
ANTIBODY_GEN_FILE_LABELS = ['igg_fha', 'igg_pt','igg_prn']
DAYS = [0, 1, 3, 7]

##### Dump the RNA-seq table from the CMI-PB website 
rule dump_cmi_pb_tables:
    output:
        ss_table = 'output/database_dump/subject_sample_days.tsv',
        ab_dump = 'output/database_dump/ab_titers_dump.tsv',
        cytof_dump = 'output/database_dump/cytof_dump.tsv',
        olink_dump = 'output/database_dump/olink_dump.tsv',
        rnaseq_dump = 'output/database_dump/rnaseq_dump.tsv'
    params:
        outdir = 'output/database_dump/'
    log:
        'logs/rule_dump_cmi_pb_tables.log'
    shell:
        """
            python scripts/dump_cmi_pb_tables.py {params.outdir}
        """

##### Clean RNA-seq data and split by day
# In addition I generate a full scale cleaned RNA-seq 
# data table and some characterizing analyses files. 
checkpoint generate_wgcna_input:
    input:
        rules.dump_cmi_pb_tables.output.rnaseq_dump
    output:
        outdir = directory('output/wgcna/input')
    log:
        'logs/checkpoint_generate_wgcna_input.log'
    shell:
        """
            mkdir -p {output}
            python scripts/generating_wgcna_input_data.py {input} {output}
        """

##### Dump the Ab titer table from the CMI-PB website
# Apply some cleaning as well
rule extract_ab_titers_wgcna:
    input:
        rules.generate_wgcna_input.output.outdir
    output:
        ab_pheno = 'output/wgcna/phenotypes/ab_titers_as_clinical_phenotypes.tsv'
    log:
        'logs/rule_extract_ab_titers_wgcna.log'
    shell: 
        """
            python scripts/extracting_ab_titer_for_WGCNA.py {output}
        """

##### Create a sample dendrogram for each split from generate_wgcna_input
rule wgcna_sample_dendrogram:
    input:
        daybased_rnaseq = 'output/wgcna/input/wgcna_input_rnaseq_day_{day}.tsv'
    output:
        dendrogram = 'output/wgcna/runs/wgcna_day_{day}.sampleClustering.png',
        dendroConst = 'output/wgcna/runs/wgcna_day_{day}.dendrogramConstructed.RData'
    log:
        'logs/rule_wgcna_sample_dendrogram_{day}.log'
    shell:
        """
            Rscript scripts/wgcna_draw_sample_dendrogram.R {input.daybased_rnaseq} {output}
        """

##### Filter samples 
rule wgcna_sample_dendrogram_filtering:
    input:
        dendroConst = rules.wgcna_sample_dendrogram.output.dendroConst,
        cutoff_txt = 'output/wgcna/manual_inputs/wgcna_day_{day}.sample_dendrogram_cutoffs.txt'
    output:
        fltExpr = 'output/wgcna/runs/wgcna_day_{day}.filteredExpression.RData'
    log:
        'logs/rule_wgcna_sample_dendrogram_filtering_{day}.log'
    shell:
        """
            Rscript scripts/wgcna_filter_samples_using_dendrogram.R {input} {output}
        """

##### Generate WGCNA module for each split from generate_wgcna_input
rule wgcna_create_rnaseq_modules:
    input:
        rules.wgcna_sample_dendrogram_filtering.output.fltExpr
    output:
        net = 'output/wgcna/runs/wgcna_day_{day}.networkConstruction-auto.RData',
        modules = 'output/wgcna/runs/wgcna_day_{day}.modules.tsv',
        eigen_expr = 'output/wgcna/runs/wgcna_day_{day}.eigengenes.expression.tsv',
        module_map = 'output/wgcna/runs/wgcna_day_{day}.modulemap.tsv'
    log:
        'logs/rule_wgcna_create_rnaseq_modules_{day}.log'
    shell:
        """
            Rscript scripts/wgcna_create_rnaseq_modules.R {input} {output}
        """

##### get eigengene files
def get_eigengene_files(wildcards):

    checkpoint_dir = checkpoints.generate_wgcna_input.get(**wildcards).output[0]

    # getting all days 
    frags_blanks = os.path.join(checkpoint_dir, 'wgcna_input_rnaseq_day_{day}.tsv')
    days = glob_wildcards(frags_blanks).day

    # returning all eigengene files 
    eg_fns = expand('output/wgcna/runs/wgcna_day_{day}.eigengenes.expression.tsv', day=days)
    return eg_fns

##### Master rule to run from checkpoint generate_wgcna_input to rule wgnca_create_rnaseq_moduls
rule wgcna_module_analyses_complete:
    input:
        get_eigengene_files
    output:
        'output/wgcna/runs/wgcna_module_analyses.complete.txt'
    log:
        'logs/rule_wgcna_module_analyses_complete.log'
    shell:
        """
            touch {output}
        """

##### Correlates mRNA modules with Ab titer information
rule wgcna_correlate_modules_and_traits:
    input:
        ss_table = rules.dump_cmi_pb_tables.output.ss_table,
        ge_rdata = rules.wgcna_sample_dendrogram_filtering.output.fltExpr,
        net_rdata = rules.wgcna_create_rnaseq_modules.output.net,
        datTraits = rules.extract_ab_titers_wgcna.output.ab_pheno
    output:
        png = 'output/wgcna/runs/wgcna_day_{day}.module_vs_igg_corr.png',
        correlations = 'output/wgcna/runs/wgcna_day_{day}.module_vs_igg_corr.tsv',
        pvals = 'output/wgcna/runs/wgcna_day_{day}.module_vs_igg_pvals.tsv'
    log:
        'logs/rule_wgcna_correlate_modules_and_traits_{day}.log'
    shell:
        """
            Rscript scripts/wgcna_correlate_modules_and_traits.R {input} {output}
        """

##### Correlates mRNA modules with Ab titer information
rule wgcna_gathering_ml_input_data:
    input:
        ge_rdata = expand('output/wgcna/runs/wgcna_day_{day}.filteredExpression.RData', day=DAYS),
        net_rdata = expand('output/wgcna/runs/wgcna_day_{day}.networkConstruction-auto.RData', day=DAYS),
        moduleTraitCor_fn = expand('output/wgcna/runs/wgcna_day_{day}.module_vs_igg_corr.tsv', day=DAYS),
        moduleTraitPvalue_fn = expand('output/wgcna/runs/wgcna_day_{day}.module_vs_igg_pvals.tsv', day=DAYS)
    output:
        expand("output/wgcna/runs/wgcna_gene_expression_members_of_sigmod_{igg_ag}.tsv", igg_ag=ANTIBODY_GEN_FILE_LABELS)
    log:
        'logs/rule_wgcna_gathering_ml_input_data.log'
    shell:
        """
            Rscript scripts/wgcna_gathering_ml_input_data.R {input} {output}
        """


##### Model ab titers using least squares, elastic net and lasso
def get_ab_name(wildcards):
    ab_label = wildcards.ab_titer
    ab_name = ab_label.replace('igg', 'IgG')
    ab_name = ab_name.split('_')
    ab_name[1] = ab_name[1].upper()
    ab_name = '-'.join(ab_name)
    return ab_name

rule model_ab_titers:
    input:
        top_genes = 'output/wgcna/runs/wgcna_gene_expression_members_of_sigmod_{ab_titer}.tsv',
        ab_data = rules.extract_ab_titers_wgcna.output.ab_pheno
    params:
        igg_assay = get_ab_name,
        trials = 300
    output:
        summary_png = 'output/wgcna/runs/sklearn_model_summary_{ab_titer}.png'
    log:
        'logs/rule_model_ab_titers_{ab_titer}.log'
    shell:
        """
            python scripts/modelling_antibody_titers_with_expression_from_top_modules.py {input} {output} {params}
        """
