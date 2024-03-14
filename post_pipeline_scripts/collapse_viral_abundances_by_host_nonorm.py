import pandas as pd
import sys
#usage: calculate_host_abundance.py <merged counts table>  <host prediction file>  <path to outfile>
#usage: calculate_host_abundance.py counts_table.txt host_prediction.tsv host_counts.txt
#input###

def main():
    #input###
    input_counts = sys.argv[1] #merged counts
    host_file = sys.argv[2] # microbial host per species level taxonomy
    out = sys.argv[3]

    #reading inputs
    profile = pd.read_csv(input_counts, sep="\t")
    samples_list = list(profile.columns)[2:]

    host_prediction = pd.read_csv(host_file, sep="\t")

    #filteing-in only viral species
    #profile_viral_species = profile[(profile.Taxon_Lineage_with_Names.str.contains("superkingdom_Viruses")) & (profile.Taxon_Lineage_with_Names.str.contains("species"))]
    profile_viral_species = profile[(profile.Taxon_Lineage_with_Names.str.contains("root_Viruses")) & (profile.Taxon_Lineage_with_Names.str.contains("species"))]

    if profile_viral_species.empty:
        header = "Taxon," + ",".join(profile.columns[2:])
        with open(out, "w") as outf:
            sys.stderr.write("No viral species present; exiting\n")
            outf.write(header + "\n")
        sys.exit()

    #merging dataframes by species taxonomy id
    species_taxa = profile_viral_species.Taxon_Lineage_with_IDs.str.split("_").str[-1]
    profile_viral_species.insert(loc=2, column='species_taxonomy', value=species_taxa.astype(int))

    profile_with_host = pd.merge(profile_viral_species, host_prediction, how='left', right_on='species_taxa', left_on='species_taxonomy')
    #sum counts/abundance per level
    profile_with_host['Host genus'].fillna("d__Unknown;p__Unknown;c__Unknown;o__Unknown;f__Unknown;g__unknown", inplace=True)
    profile_with_host.drop(columns=['species_taxonomy', 'species_taxa', 'Taxon_Lineage_with_Names', 'Taxon_Lineage_with_IDs'], inplace=True)

    profile_d = profile_with_host.copy()
    profile_p = profile_with_host.copy()
    profile_c = profile_with_host.copy()
    profile_o = profile_with_host.copy()
    profile_f = profile_with_host.copy()
    profile_g = profile_with_host.copy()

    profile_d['Taxon'] = profile_d['Host genus'].str.split(";").str[0:1].str.join(";")
    profile_p['Taxon'] = profile_p['Host genus'].str.split(";").str[0:2].str.join(";")
    profile_c['Taxon'] = profile_c['Host genus'].str.split(";").str[0:3].str.join(";")
    profile_o['Taxon'] = profile_o['Host genus'].str.split(";").str[0:4].str.join(";")
    profile_f['Taxon'] = profile_f['Host genus'].str.split(";").str[0:5].str.join(";")
    profile_g['Taxon'] = profile_g['Host genus'].str.split(";").str[0:6].str.join(";")

    profile_d.drop(columns=['Host genus'], inplace=True)
    profile_p.drop(columns=['Host genus'], inplace=True)
    profile_c.drop(columns=['Host genus'], inplace=True)
    profile_o.drop(columns=['Host genus'], inplace=True)
    profile_f.drop(columns=['Host genus'], inplace=True)
    profile_g.drop(columns=['Host genus'], inplace=True)

    #creating out file
    headers = ["Taxon"] + samples_list
    d = profile_d.groupby('Taxon').sum().reset_index().loc[:, headers]
    p = profile_p.groupby('Taxon').sum().reset_index().loc[:, headers]
    c = profile_c.groupby('Taxon').sum().reset_index().loc[:, headers]
    o = profile_o.groupby('Taxon').sum().reset_index().loc[:, headers]
    f = profile_f.groupby('Taxon').sum().reset_index().loc[:, headers]
    g = profile_g.groupby('Taxon').sum().reset_index().loc[:, headers]

    #creating out file
    taxon = pd.concat([d, p, c, o, f, g])
    #if taxon.iloc[1].dtype=='float': #normalize to 1 if relative abundance was given
    #    taxon = pd.concat([d/d.sum() , p/p.sum(), c/c.sum() , o/o.sum(), f/f.sum(), g/g.sum()])

    #taxon.index.name = "Taxon"
    taxon.to_csv(out, index=False)


if __name__ == "__main__":
    main()
