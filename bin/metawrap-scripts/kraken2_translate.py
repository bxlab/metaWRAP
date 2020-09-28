#!/usr/bin/env python
import sys
import os


def get_full_name(taxid, names_map, ranks_map):
    """
    Generate full taxonomic lineage from a taxid
    Args:
        taxid (str): taxid
        names_map (dict[str:str]): the taxonomic names for each taxid node
        ranks_map (dict[str:str]): the parent node for each taxid
    Returns:
        taxonomy (str): full taxonomic lineage string

    """
    taxid_lineage = list()
    while True:
        taxid_lineage.append(taxid)
        if taxid not in ranks_map:
            break
        new_taxid = ranks_map[taxid]
        if taxid == new_taxid:
            break
        else:
            taxid = new_taxid

    names_lineage = list()
    for taxid in taxid_lineage:
        name = names_map[taxid]
        names_lineage.append(name)

    taxonomy = ";".join(reversed(names_lineage))
    return taxonomy


def load_kraken_db_metadata(kraken2_db):
    """
    Load NCBI taxonomic name mappings to be able to convert taxids into taxonomic strings
    Args:
        kraken2_db (str): path to kraken2 standard database location
    Returns:
        names_map (dict[str:str]): the taxonomic names for each taxid node
        ranks_map (dict[str:str]): the parent node for each taxid

    """
    print("Loading NCBI node names")
    names_path = os.path.join(kraken2_db, "taxonomy", "names.dmp")
    names_map = dict()
    with open(names_path) as input:
        for line in input:
            cut = line.rstrip().split("\t")
            taxid = cut[0]
            name = cut[2]
            type = cut[6]
            if type == "scientific name":
                names_map[taxid] = name

    print("Loading NCBI taxonomic ranks")
    ranks_path = os.path.join(kraken2_db, "taxonomy", "nodes.dmp")
    ranks_map = dict()
    with open(ranks_path) as input:
        for line in input.readlines():
            cut = line.rstrip().split("\t")
            taxid = cut[0]
            parent_taxid = cut[2]
            rank = cut[4]
            ranks_map[taxid] = parent_taxid
    return names_map, ranks_map


def translate_kraken2_annotations(annotation_file=None, kraken2_db=None, output=None):
    """
    Translate the kraken2 annotations from taxids to taxonomy strings
    Args:
        annotation_file (str): kraken2 output file
        kraken2_db (str): path to kraken2 standard database
        output (str): path to final file with translated annotations
    Returns: None

    """
    print("Translating kraken2 annotations")
    if os.path.isfile(output):
        print("Looks like the translated kraken2 output already exists. Skipping...")
        return None

    names_map, ranks_map = load_kraken_db_metadata(kraken2_db)
    print("Writing translated taxonomy names to %s" % output)
    with open(output, "w") as output:
        with open(annotation_file) as input:
            for line in input:
                cut = line.rstrip().split("\t")
                contig = cut[1]
                if cut[0] == "U":
                    taxonomy = ""
                else:
                    taxid = cut[2].split()[-1][:-1]
                    taxonomy = get_full_name(taxid, names_map, ranks_map)
                output.write("%s\t%s\n" % (contig, taxonomy))


def main():
    """
    Load in arguments and start analysis
    Returns: None

    """
	
    database_location = sys.argv[1]
    kraken_file = sys.argv[2]
    output_file = sys.argv[3]
    print("Translating kraken2 annotations from %s, using metadata from the kraken2 database in %s; saving to %s" \
          % (kraken_file, database_location, output_file))
    translate_kraken2_annotations(annotation_file=kraken_file, kraken2_db=database_location, output=output_file)
    

if __name__ == '__main__':
    """ Launch script
    """
    main()

