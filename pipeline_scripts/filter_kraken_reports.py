"""
STEP ONE
Import required libraries.
Also get arguments from the command line
1. Kraken report file (full path)
2. Kraken database (full path)
3. The threshold for "maximum Kraken-based coverage of any genome in this strain" above which a BACTERIAL strain should be counted as a TRUE POSITIVE
4. same as #3 but for VIRAL strain
5. The threshold for "maximum unique minimizers in the sample mapped to any genome in this strain" above which a BACTERIAL strain should be counted as a TRUE POSITIVE
6. same as #5 but for VIRAL strain
7. Coverage threshold for archaeal strain
8. Coverage threshold for eukaryotic strain
9. Minimizer threshold for archaeal strain
10. Minimizer threshold for eukaryotic strain
Note: for each domain, BOTH the cov threshold and the minimizer threshold have to be "passed" for a strain to be counted as true positive.
"""

import sys
import pandas as pd
import numpy as np

kraken_report, kraken_db, max_cov_bacteria, max_cov_virus, \
max_minimizers_bacteria, max_minimizers_virus, max_cov_arc, \
max_cov_euk, max_minimizers_arc, max_minimizers_euk = \
sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), \
int(sys.argv[5]), int(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), \
int(sys.argv[9]), int(sys.argv[10])

"""
STEP TWO
Define some functions that will be required for STEP THREE
"""

def make_dicts(nodes_file):
  with open(nodes_file, 'r') as infile:
    # make a child to parent dictionary
    # and a taxid to rank dictionary
    child_to_parent = {}
    taxid_to_rank = {}
    for line in infile:
      line=line.rstrip('\n').split('\t')
      child, parent, rank = line[0], line[2], line[4]
      child_to_parent[child] = parent
      taxid_to_rank[child] = rank
  return child_to_parent, taxid_to_rank

def taxid_to_desired_rank(taxid, desired_rank, child_parent, taxid_rank):
  # look up the specific taxid,
  # build the lineage using the dictionaries
  # stop at the desired rank and return the taxid
  lineage = [[taxid, taxid_rank[taxid]]]
  if taxid_rank[taxid] == desired_rank:
    return taxid
  child, parent = taxid, None
  if child == '0':
    return 'unclassified'
  while not parent == '1':
    # print(child, parent)
    # look up child, add to lineage
    parent = child_parent[child]
    rank = taxid_rank[parent]
    lineage.append([parent, rank])
    if rank == desired_rank:
      return parent
    child = parent # needed for recursion
  return 'error - taxid above desired rank, or not annotated at desired rank'

"""
STEP THREE
Adapted from /labs/asbhatt/yishay/scripts/kraken_strain_filter_low_rank.py
Make a modified version of the Kraken report that:
1) Only includes lines representing strain and strains
2) Has two extra columns: 1) coverage, 2) strain level taxid
"""
# read in Kraken report, merge with inspect file to get info to calculate coverage

kraken_file = pd.read_csv(kraken_report, sep="\t", names=['fraction','fragments', 'assigned','minimizers','uniqminimizers', 'classification_rank','ncbi_taxa','scientific name'])
inspect_fname = kraken_db + '/inspect.out'
inspect_file = pd.read_csv(inspect_fname, sep="\t", names=['frac','minimizers_clade', 'minimizers_taxa', 'rank','ncbi_taxonomy','sci_name'])

kraken_file = kraken_file.merge(inspect_file,left_on='ncbi_taxa', right_on='ncbi_taxonomy')

# calculate coverage
kraken_file['cov'] = kraken_file['uniqminimizers']/kraken_file['minimizers_taxa']

# assign coverage as NA if there are fewer than 5 minimizers assigned to this taxon
kraken_file.loc[kraken_file.minimizers_taxa < 5,'cov']=np.nan

# filter kraken_file to strain only
strain_kraken = kraken_file.copy()[kraken_file['classification_rank'].str.startswith('S1', na=False)]

# add a column for strain level taxid - use functions defined above
child_parent, taxid_rank = make_dicts(kraken_db + '/taxonomy/nodes.dmp')
strain_kraken['strain_level_taxa'] = strain_kraken.apply(lambda x: taxid_to_desired_rank(str(x['ncbi_taxonomy']), 'strain', child_parent, taxid_rank), axis=1)
strain_kraken.to_csv(kraken_report + ".strain", sep="\t", index=False)

"""
STEP FOUR
Read in the modified Kraken report and make two dictionaries.
1) strain_to_superkingdom: strain taxid to superkingdom
2) strain_genome_info:
strain taxid to list of: [taxid, cov, rank, unique_minimizers] for every taxon
found by Kraken and falling under that strain taxid
"""
with open(kraken_report + '.strain', 'r') as infile:
  # figure out the column number for all the features of interest

  header=infile.readline().rstrip('\n').split('\t')

  unique_minimizers_col, cov_col, strain_col, rank_col, taxid_col = \
  None, None, None, None, None

  for i in range(len(header)):
    if header[i] == 'cov':
      cov_col = i
    elif header[i] == 'strain_level_taxa':
      strain_col = i
    elif header[i] == 'rank':
      rank_col = i
    elif header[i] == 'uniqminimizers':
      unique_minimizers_col = i
    elif header[i] == 'ncbi_taxa':
      taxid_col = i

  assert unique_minimizers_col != None
  assert cov_col != None
  assert strain_col != None
  assert rank_col != None
  assert taxid_col != None

  # now go through report!
  # first initialize the dictionaries
  strain_to_superkingdom = {}
  strain_genome_info = {}

  i = 0
  for line in infile: # skipped header already
    i += 1
    line=line.rstrip('\n').split('\t')

    unique_minimizers, cov, strain_taxid, rank, taxid = \
    int(line[unique_minimizers_col]), line[cov_col], line[strain_col], line[rank_col], \
    line[taxid_col]

    if cov == '':
      cov = 0
    else:
      cov = float(cov)

    # add to dictionaries
    if rank == 'S1':
      # can add to the strain_to_superkingdom dictionary
      # look up superkingdom
      superkingdom = taxid_to_desired_rank(strain_taxid, 'superkingdom', child_parent, taxid_rank)
      strain_to_superkingdom[strain_taxid] = superkingdom
      # if i > 100 and i < 200:
        # print('superkingdom', i)
        # print(strain_to_superkingdom[strain_taxid])

    # regardless, add to the strain_genome_info dictionary
    if not strain_taxid in strain_genome_info:
      # initialize
      strain_genome_info[strain_taxid] = [[taxid, cov, rank, unique_minimizers]]
    else:
      strain_genome_info[strain_taxid].append([taxid, cov, rank, unique_minimizers])

#    if i > 100 and i < 200:
      # print('other', i)
      # print(strain_genome_info[strain_taxid])

"""
STEP FIVE
Go through strain_genome_info and figure out - using strain_to_superkingdom as well -
which strain should be "kept" in the Kraken report.
Output a file containing info about these decisions.
Here we will use the parameters that were passed in:
max_cov_bacteria, max_cov_virus, max_minimizers_bacteria, max_minimizers_virus
Recall: if max_cov is passed, definitely keep. If not, keep if max_minimizers threshold is passed.
Recall - we only care about the coverages/minimizer values for the GENOMES in
each strain - i.e., the entries in strain_genome_info that don't have a 'lower rank'
right afterwards.
Note - will also make a dictionary that stores this info - taxid_to_lowest_rank.
"""
strain_to_keep = set()
taxid_to_lowest_rank = {}

out_fname = kraken_report + '.filtering_decisions.txt'

with open(out_fname, 'w') as outfile:

  # write out a header first
  header = ['strain_taxid', 'superkingdom', 'max_cov', 'max_uniq_minimizers', 'kept']
  outfile.write('\t'.join(header) + '\n')

  for strain in strain_genome_info:
    lines_to_consider = strain_genome_info[strain]
    relevant_coverages, relevant_unique_minimizer_vals = [], []
    for i in range(len(lines_to_consider)):

      taxid_to_lowest_rank[lines_to_consider[i][0]] = 'False'

      if i == (len(lines_to_consider) - 1): # special case - last line in the Kraken report for this strain taxid
        relevant_coverages.append(lines_to_consider[i][1])
        relevant_unique_minimizer_vals.append(lines_to_consider[i][3])
        # update taxid_to_lowest_rank
        taxid_to_lowest_rank[lines_to_consider[i][0]] = 'True'
      else: # not the last line in the Kraken report for this strain taxid
        this_rank = lines_to_consider[i][2]
        next_rank = lines_to_consider[i+1][2]
        if this_rank == 'S1':
          continue # since it's not the last line in the Kraken report, next_rank must be "lower"
          # sanity check
          assert len(next_rank) > len(this_rank)
        else: # we're dealing with S1, S2, etc.
          if int(this_rank[1:]) < int(next_rank[1:]): # e.g., S1 < S2
            continue
          else: # e.g., S3 > S2
            relevant_coverages.append(lines_to_consider[i][1])
            relevant_unique_minimizer_vals.append(lines_to_consider[i][3])
            # update taxid_to_lowest_rank
            taxid_to_lowest_rank[lines_to_consider[i][0]] = 'True'

    # now we can get the maximum coverage from the relevant coverages
    max_cov = max(relevant_coverages)
    # also get the maximum "# unique minimizers"
    max_minimizers = max(relevant_unique_minimizer_vals)

    # should we keep this strain?
    # figure out based on superkingdom
    superkingdom = strain_to_superkingdom[strain]

    # superkingdom
    # Bacteria        609216830
    # Archaea         439684927
    # Bamfordvirae    57932934
    # Heunggongvirae  1238430944
    # Loebvirae       94436553
    # Orthornavirae   104708768
    # Pararnavirae    1921562045
    # r__Monodnaviria_Unclassified    414168241
    # r__Unclassified 1346054397
    # Sangervirae     107962225
    # Shotokuvirae    1793913686

    #if superkingdom == '2': # bacteria
    if superkingdom == '609216830': # bacteria
      if (max_cov >= max_cov_bacteria) and (max_minimizers >= max_minimizers_bacteria):
        strain_to_keep.add(strain)

    #elif superkingdom == '2157': # archaea
    elif superkingdom == '439684927': # archaea
      if (max_cov >= max_cov_arc) and (max_minimizers >= max_minimizers_arc):
        strain_to_keep.add(strain)

    #elif superkingdom == '2759': # eukaryotes
    #  if (max_cov >= max_cov_euk) and (max_minimizers >= max_minimizers_euk):
    #    strain_to_keep.add(strain)

    #elif superkingdom == '10239': # viruses
    elif superkingdom in ['57932934', '1238430944', '94436553', '104708768', '1921562045', '414168241', '1346054397', '107962225', '1793913686']: # viruses
      if (max_cov >= max_cov_virus) and (max_minimizers >= max_minimizers_virus):
        strain_to_keep.add(strain)

    else:
      strain_to_keep.add(strain)

    if strain in strain_to_keep:
      keep = 'True'
    else:
      keep = 'False'
    outfile.write('\t'.join([strain, superkingdom, str(max_cov), str(max_minimizers), keep]) + '\n')

"""
STEP SIX
Go through the original Kraken report and make a new filtered version,
where lines from strain to filter out are removed.
"""
with open(kraken_report, 'r') as infile:
  with open(kraken_report + '.filtered', 'w') as outfile:
    for line in infile:
      line=line.rstrip('\n').split('\t')
      if line[5].startswith('S1'):
        # look up the taxonomy id - what is the strain level taxid?
        strain_taxid = taxid_to_desired_rank(line[6], 'strain', child_parent, taxid_rank)
        # do we care about this strain?
        if strain_taxid in strain_to_keep:
          outfile.write('\t'.join(line) + '\n')
      else:
        outfile.write('\t'.join(line) + '\n')
