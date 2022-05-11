# command line arguments
# 1) the original merged output to pare down
# must be one of: counts.txt, counts_norm_out_of_tot.txt, counts_norm_out_of_brack_classified.txt
# 2) the desired rank (e.g., genus, species)

import sys

orig, rank = sys.argv[1], sys.argv[2]

# construct outfile name
orig_prefix = orig[0:orig.rfind('.')]
outfname = orig_prefix + '_' + rank + '.txt'

with open(orig, 'r') as infile:
  with open(outfname, 'w') as outfile:
    for line in infile:
      line_split=line.split('\t')
      tax_name = line_split[0].split('|')
      tax_rank = tax_name[len(tax_name)-1].split('_')[0]
      if tax_rank == rank:
        outfile.write(line)
