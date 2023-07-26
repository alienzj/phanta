##### STEP ZERO - Import required modules

from os.path import join
import sys
import snakemake

##### STEP ONE - Get information from config file.

def get_sample_reads(sample_file):
  sample_reads = {}
  paired_end = ''
  with open(sample_file, 'r') as sf:
    for line in sf:
      #print(line)
      line = line.rstrip('\n').split('\t')
      if len(line) == 1 or line[0] == 'Sample' or line[0] == '#Sample' or line[0].startswith('#'):
        continue
      else: # not header
        sample = line[0]
        #print(sample)
      if (len(line) == 3): # paired end
        reads = [line[1], line[2]]
        if paired_end != '' and not paired_end:
          sys.exit('All samples must be paired or single ended.')
        paired_end = True

      elif (len(line) == 2): # single end specified
        reads = [line[1]]
        if paired_end != '' and paired_end:
          sys.exit('All samples must be paired or single ended.')
        paired_end = False

      if sample in sample_reads:
        raise ValueError("Non-unique sample encountered!")
      else:
        #print('here')
        sample_reads[sample] = reads

    return (sample_reads, paired_end)

# call function, determine whether reads are paired or single ended
sample_reads, paired_end = get_sample_reads(config['sample_file'])

#print(sample_reads, paired_end)

if paired_end:
  paired_string = '--paired'
else:
  paired_string = ''
sample_names = sample_reads.keys()

# also determine whether gzipped
# print(config['gzipped'] == True)
if config['gzipped'] == True:
   gzipped_string = '--gzip-compressed'
else:
   gzipped_string = ''

# also read in desired confidence threshold for Kraken
confidence_threshold = config['confidence_threshold']

# define output directory
outdir = config['outdir']

# get set up to run the single ended Kraken2 if that is the user's desire
single_end_files = []
if config["single_end_krak"]:
    single_end_files.append(expand(join(outdir, "classification/single_end/fwd_run/{samp}.krak"), samp=sample_names))
    single_end_files.append(expand(join(outdir, "classification/single_end/rev_run/{samp}.krak"), samp=sample_names))

##### STEP TWO - Define desired pipeline outputs.

rule all:
  input:
    join(outdir, "final_merged_outputs/total_reads.tsv"),
    join(outdir, "final_merged_outputs/counts.txt"),
    join(outdir, "final_merged_outputs/relative_read_abundance.txt"),
    join(outdir, "final_merged_outputs/relative_taxonomic_abundance.txt"),
    join(outdir, "pipeline_completed.txt"),
    single_end_files

##### STEP THREE - Run Kraken2, and filter report based on user-defined thresholds.

rule kraken:
  input:
    reads = lambda wildcards: sample_reads[wildcards.samp]
  output:
    krak = join(outdir, "classification/{samp}.krak"),
    krak_report = join(outdir, "classification/{samp}.krak.report")
  params:
    db = config['database'],
    paired_string = paired_string,
    gzipped_string = gzipped_string,
    confidence_threshold = confidence_threshold
  threads: config['class_threads']
  resources:
    mem_mb=config['class_mem_mb'],
    runtime=6*60,
  shell: """
    kraken2 --db {params.db} --threads {threads} --output {output.krak} \
    --report {output.krak_report} --report-minimizer-data {params.paired_string} \
    {params.gzipped_string} {input.reads} --confidence {params.confidence_threshold}
    """

rule filter_kraken:
  input:
    krak_report = join(outdir, "classification/{samp}.krak.report")
  output:
    krak_species = temp(join(outdir, "classification/{samp}.krak.report.species")),
    #krak_species_final = join(outdir, "classification/{samp}.krak.report.species.final"),
    krak_report_filtered = join(outdir, "classification/{samp}.krak.report.filtered"),
    filtering_decisions = join(outdir, "classification/{samp}.krak.report.filtering_decisions.txt")
  params:
    repo_dir = config['pipeline_directory'],
    db = config['database'],
    cov_thresh_bacterial = config['cov_thresh_bacterial'],
    cov_thresh_viral = config['cov_thresh_viral'],
    cov_thresh_arc = config['cov_thresh_arc'],
    cov_thresh_euk = config['cov_thresh_euk'],
    minimizer_thresh_bacterial = config['minimizer_thresh_bacterial'],
    minimizer_thresh_viral = config['minimizer_thresh_viral'],
    minimizer_thresh_arc = config['minimizer_thresh_arc'],
    minimizer_thresh_euk = config['minimizer_thresh_euk']
  shell: """
    python {params.repo_dir}/pipeline_scripts/filter_kraken_reports.py {input.krak_report} {params.db} \
    {params.cov_thresh_bacterial} {params.cov_thresh_viral} {params.minimizer_thresh_bacterial} \
    {params.minimizer_thresh_viral} \
    {params.cov_thresh_arc} {params.cov_thresh_euk} {params.minimizer_thresh_arc} \
    {params.minimizer_thresh_euk} 
  """

# need the below rule for protection against Bracken erroring if all species with
# # >= thresh reads have been removed from the report
rule process_filtered_kraken:
  input:
    krak_report_filtered = join(outdir, "classification/{samp}.krak.report.filtered")
  output:
    completed = temp(join(outdir, "processed_filtered_kraken/{samp}.txt")),
  params:
    repo_dir = config['pipeline_directory'],
    threshold = config['filter_thresh'],
    db = config['database']
  shell: """
    python {params.repo_dir}/pipeline_scripts/process_filtered_kraken.py {input.krak_report_filtered} \
    {params.threshold} {params.db}
    touch {output.completed}
    """

# single ended Kraken2 (optional - user will specify whether also to run this, in config file)
rule kraken_single_end:
  input:
    fwd_reads = lambda wildcards: sample_reads[wildcards.samp][0],
    rev_reads = lambda wildcards: sample_reads[wildcards.samp][1]
  output:
    fwd_krak = join(outdir, "classification/single_end/fwd_run/{samp}.krak"),
    rev_krak = join(outdir, "classification/single_end/rev_run/{samp}.krak")
  params:
    db = config['database'],
    gzipped_string = gzipped_string,
    confidence_threshold = confidence_threshold
  threads: config['class_threads']
  resources:
    mem_mb=config['class_mem_mb'],
    runtime=6*60,
  shell: """
    kraken2 --db {params.db} --threads {threads} --output {output.fwd_krak} \
    {params.gzipped_string} {input.fwd_reads} --confidence {params.confidence_threshold}
    kraken2 --db {params.db} --threads {threads} --output {output.rev_krak} \
    {params.gzipped_string} {input.rev_reads} --confidence {params.confidence_threshold}
    """

##### STEP FOUR - Run Bracken on filtered Kraken2 report.

rule bracken:
  input:
    krak_report = join(outdir, "classification/{samp}.krak.report.filtered"),
    good_to_go = join(outdir, "processed_filtered_kraken/{samp}.txt")
  output:
    brack_report_1 = join(outdir, "classification/{samp}.krak.report.filtered.bracken"),
    brack_report_2 = join(outdir, "classification/{samp}.krak.report_bracken_species.filtered")
  params:
    db = config['database'],
    readlen = config['read_length'],
    threshold = config['filter_thresh'],
    possible_1 = join(outdir, "classification/{samp}.krak.report.filtered.bracken.temp"),
    possible_2 = join(outdir, "classification/{samp}.krak.report_bracken_species.filtered.temp")
  threads: 1
  resources:
    mem_mb = 8*1024,
    runtime = 1*60,
  shell: """
    # protection against Bracken error
    [ -f {params.possible_1} ] && \
    ( cp {params.possible_1} {output.brack_report_1};
    cp {params.possible_2} {output.brack_report_2} ) \
    || bracken -d {params.db} -i {input.krak_report} \
    -o {output.brack_report_1} -r {params.readlen} \
    -l 'S' -t {params.threshold}
    """

rule samples_failed_bracken: # make a list of samples that failed Bracken
  input:
    expand(join(outdir, "processed_filtered_kraken/{samp}.txt"), samp=sample_names)
  output:
    failed=join(outdir, "classification/samples_that_failed_bracken.txt")
  params:
    classdir = join(outdir, "classification")
  shell: """
    count=`find {params.classdir} -maxdepth 1 -name "*.krak.report.filtered.bracken.scaled.temp" | wc -l`
    if [[ $count != 0 ]]
    then
      ls {params.classdir}/*.krak.report.filtered.bracken.scaled.temp | rev | \
      cut -d'/' -f 1 | rev | cut -d'.' -f 1 > {output.failed}
    else
      touch {output.failed}
    fi
    """

##### STEP FIVE - Merge final Bracken reports into usable tables - 1) counts table, 2) normalized counts - normalize by TOTAL READS in sample, 3) normalized counts - normalize by BRACKEN-CLASSIFIED READS in sample.

rule prepare_to_merge_counts: # do this rule for each Bracken report individually
  input:
    brack_report = join(outdir, "classification/{samp}.krak.report_bracken_species.filtered"),
    failed_file = join(outdir, "classification/samples_that_failed_bracken.txt")
  output:
    brack_to_merge = temp(join(outdir, "classification/{samp}.krak.report_bracken_species.filtered.to_merge"))
  params:
    repo_dir = config['pipeline_directory'],
    db = config['database']
  shell: """
    python {params.repo_dir}/pipeline_scripts/prep_to_merge_counts.py \
    {input.brack_report} {params.db} {input.failed_file} {wildcards.samp}
    """

# Now we will calculate the total reads in each sample, for normalization
# Also calculate the total reads ultimately classified by Bracken

rule tot_reads_per_sample: # will need for normalization
  input:
    brack_report = join(outdir, "classification/{samp}.krak.report_bracken_species.filtered")
  output:
    tot_reads_file = temp(join(outdir, "classification/{samp}_total_reads.txt"))
  params:
    repo_dir = config['pipeline_directory'],
    krak_brack_dir = join(outdir, "classification")
  shell: """
    bash {params.repo_dir}/pipeline_scripts/calc_total_reads.sh {wildcards.samp} \
    {params.krak_brack_dir} {output.tot_reads_file}
    """

rule tot_reads_all: # aggregate outputs of the previous rule
  input:
    expand(join(outdir, "classification/{samp}_total_reads.txt"), samp=sample_names)
  output:
    outf=join(outdir, "final_merged_outputs/total_reads.tsv")
  params:
    classdir=join(outdir, "classification")
  shell: """
    echo -e 'Samp_Name\tTot_Samp_Reads\tUnassigned_Step_One\tAssigned_Step_Three\tUnassigned_Step_Three' > {output.outf}
    cat {params.classdir}/*total_reads.txt >> {output.outf}
    """

rule prepare_to_merge_normed: # normalized versions of to_merge files produced in an earlier rule
  input:
    counts_files=expand(join(outdir, "classification/{samp}.krak.report_bracken_species.filtered.to_merge"), samp=sample_names),
    tot_reads_file=join(outdir, "final_merged_outputs/total_reads.tsv"),
    failed_file = join(outdir, "classification/samples_that_failed_bracken.txt")
  output:
    temp(expand(join(outdir, "classification/{samp}.krak.report_bracken_species.filtered.to_merge.norm_brack"), samp=sample_names))
  params:
    repo_dir = config['pipeline_directory'],
    classdir=join(outdir, "classification"),
    failed_file = join(outdir, "classification/samples_that_failed_bracken.txt")
  shell: """
    python {params.repo_dir}/pipeline_scripts/prep_to_merge_normed.py \
    {input.tot_reads_file} {params.classdir} {input.failed_file}
    """

rule merge_counts_normed: # make the 3 tables we want! :)
  input:
    expand(join(outdir, "classification/{samp}.krak.report_bracken_species.filtered.to_merge"), samp=sample_names),
    expand(join(outdir, "classification/{samp}.krak.report_bracken_species.filtered.to_merge.norm_brack"), samp=sample_names),
    failed_file = join(outdir, "classification/samples_that_failed_bracken.txt")
  output:
    list1=temp(join(outdir, "classification/counts_tables.txt")),
    list2=temp(join(outdir, "classification/norm_brack_tables.txt")),
    counts=join(outdir, "final_merged_outputs/counts.txt"),
    norm_brack=join(outdir, "final_merged_outputs/relative_read_abundance.txt")
  params:
    repo_dir = config['pipeline_directory'],
    classdir=join(outdir, "classification"),
    finaldir=join(outdir, "final_merged_outputs")
  shell: """
    ls {params.classdir}/*to_merge | rev | cut -d'/' -f 1 | rev > {params.classdir}/counts_tables.txt
    ls {params.classdir}/*norm_brack | rev | cut -d'/' -f 1 | rev > {params.classdir}/norm_brack_tables.txt
    Rscript {params.repo_dir}/pipeline_scripts/merging_bracken_tables.R {params.classdir} {params.finaldir} \
    {input.failed_file} {output.list1} {output.list2}
    """

##### STEP SIX - Correct species abundances reported by Bracken for genome length.

rule scale_bracken:
  input:
    bracken_report = join(outdir, "classification/{samp}.krak.report.filtered.bracken"),
    filtering_decisions = join(outdir, "classification/{samp}.krak.report.filtering_decisions.txt")
  output:
    scaled_bracken = join(outdir, "classification/{samp}.krak.report.filtered.bracken.scaled")
  params:
    repo_dir = config['pipeline_directory'],
    db = config['database'],
    readlen = config['read_length'],
    paired = paired_end,
    possible = join(outdir, "classification/{samp}.krak.report.filtered.bracken.scaled.temp")
  threads: 1
  shell: """
    [ -f {params.possible} ] && \
    cp {params.possible} {output.scaled_bracken} || \
    python {params.repo_dir}/pipeline_scripts/scale_bracken.py {params.db} {input.bracken_report} \
    {input.filtering_decisions} {params.readlen} {params.paired}
    """

rule merge_scaled_bracken:
  input:
    expand(join(outdir, "classification/{samp}.krak.report.filtered.bracken.scaled"), samp=sample_names),
    failed_file = join(outdir, "classification/samples_that_failed_bracken.txt")
  output:
    list=temp(join(outdir, "classification/scaled_reports.txt")),
    merged_temp=temp(join(outdir, "classification/relative_taxonomic_abundance_temp.txt")),
    merged_final=join(outdir, "final_merged_outputs/relative_taxonomic_abundance.txt")
  params:
    repo_dir = config['pipeline_directory'],
    classdir=join(outdir, "classification"),
    finaldir=join(outdir, "final_merged_outputs"),
    db=config['database']
  shell: """
    ls {params.classdir}/*scaled | rev | cut -d'/' -f 1 | rev > {params.classdir}/scaled_reports.txt
    Rscript {params.repo_dir}/pipeline_scripts/merging_scaled_bracken.R {params.classdir} {output.list} {input.failed_file}
    python {params.repo_dir}/pipeline_scripts/merging_scaled_bracken.py {params.db} {output.merged_temp} {params.finaldir}
    """

##### STEP SEVEN - Move and/or delete certain "intermediate" files.

rule deal_with_intermediate:
  input: # all the outputs except pipeline_completed.txt
    join(outdir, "final_merged_outputs/total_reads.tsv"),
    join(outdir, "final_merged_outputs/counts.txt"),
    join(outdir, "final_merged_outputs/relative_read_abundance.txt"),
    join(outdir, "final_merged_outputs/relative_taxonomic_abundance.txt")
  output:
    completed=join(outdir, "pipeline_completed.txt")
  params:
    delete=config['delete_intermediate'],
    classdir=join(outdir, "classification"),
    intdir=join(outdir, "classification/intermediate")
  shell: """
    mkdir {params.intdir}
    mv {params.classdir}/*krak {params.intdir}
    mv {params.classdir}/*krak.report {params.intdir}
    mv {params.classdir}/*krak.report.filtered {params.intdir}
    mv {params.classdir}/*krak.report.filtered.bracken {params.intdir}
    mv {params.classdir}/*krak.report.filtering_decisions.txt {params.intdir}
    if [[ {params.delete} == 'True' ]]; then
      rm -r {params.intdir}
    fi
    touch {output.completed}
    """
