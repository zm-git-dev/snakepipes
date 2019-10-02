#!/usr/bin/env python

import sys
import pyBigWig
import requests
import re
import argparse
from argparse import RawTextHelpFormatter
import itertools
import gzip


def parse_arguments(defaults):
    parser = argparse.ArgumentParser(description="Convert chromosome names for bigwig files between ensembl, gencode and UCSC naming schemes\n"
                                     "Per default it writes to the same location as original file, however with a modified filename:\n"
                                     "eg. test.bw --> test.[toFormat]_chroms.bw\n"
                                     "Change this with the -o option!\n\n"
                                     "Mapping tables are taken from https://github.com/dpryan79/ChromosomeMappings\n\n"
                                     "Provided mapping options need to exactly match an existing file\n"
                                     "[GENOME]_[FROM_FORMAT]2[TO_FORMAT].txt in this repo!",
                                     usage='$ convertChroms BIGWIG', formatter_class=RawTextHelpFormatter)

    parser.add_argument('bw_in_filename',
                        nargs='+',
                        metavar='BIGWIG',
                        help='single bigwig file or list of files that will be converted')

    parser.add_argument('--genome', '-g',
                        action='store',
                        dest='genome',
                        help="Genome version of original bigwig \n"
                        "(GRCm37|GRCm38|GRCh37|GRCh38|BDGP6|dm3|GRCz10|GRCz11|\n"
                        "JGI_4.2|MEDAKA1|R64-1-1|WBcel235|Zv9|galGal4|rn5|rn6)\n"
                        "(default: %(default)s)",
                        default=defaults["genome"])

    parser.add_argument('--fromFormat', '-f',
                        action='store',
                        dest='from_format',
                        help='Chr naming format of original bigwig (ensembl|gencode|UCSC) (default: %(default)s)',
                        default=defaults["fromFormat"])

    parser.add_argument('--toFormat', '-t',
                        action='store',
                        dest='to_format',
                        help='Chr naming format of converted bigwig (ensembl|gencode|UCSC) (default: %(default)s)',
                        default=defaults["toFormat"])

    parser.add_argument('--outFileName', '-o',
                        action='store',
                        nargs='*',
                        dest='bw_out_filename',
                        help='Output filename (default: %(default)s)',
                        default=defaults["bw_out_filename"])

    parser.add_argument('--baseURL', '-u',
                        action='store',
                        dest='base_url',
                        help="base url where the mapping tables can be found (default: %(default)s)\n"
                        "Local files can be given with \'file://[BASE_DIR]/\'",
                        default=defaults["base_url"])

    parser.add_argument('--verbose', '-v',
                        action='store_true',
                        dest='verbose',
                        help='Be more verbose where possible (default: %(default)s)',
                        default=defaults["verbose"])
    
    parser.add_argument('--gzip',
                        action='store_true',
                        dest='gzipped',
                        help='File is a gzipped text file (default: %(default)s)',
                        default=defaults["gzipped"])
    

    return parser


def get_chromosome_mapping(genome="GRCm38", from_format="ensembl", to_format="UCSC", verbose=True, base_url='https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/'):
    """
    creates a dict with chromosome name mappings according to provided conversion formats

    default base URL access a github repo with conversion files,
    but you can also give eg. a path to local directory
    """
    mapping_file = genome + '_' + from_format + '2' + to_format + '.txt'

    if re.match('^file:[/]+.*', base_url):
        base_url = re.sub("file:[/]*(/.*)", "\\1", base_url)

    if verbose:
        print("load mapping table (" + mapping_file + ') from ' + base_url)

    tab = None

    if re.match('^https?://.*', base_url):
        try:
            r = requests.get(base_url + '/' + mapping_file)
            r.raise_for_status()
        except requests.exceptions.RequestException as e:
            print("\n", e, "\n\nPlease provide correct name (GENOME, FROM_FORMAT, TO_FORMAT) for a mapping table!\n")
            sys.exit(1)
        tab = r.text
    elif re.match('^[/]+.*', base_url):
        try:
            tab = open(base_url + '/' + mapping_file).read()
        except IOError as e:
            print("\n", e, "\n\nPlease provide a correct name (GENOME, FROM_FORMAT, TO_FORMAT) for a mapping table!\n")
            sys.exit(1)
    else:
        print("\nPlease provide a correct BASE_URL for a mapping table!\n")
        sys.exit(1)

    mapping_table = {}
    for ent in tab.split("\n"):
        if len(ent) == 0:
            continue
        pair = ent.split("\t")
        if (len(pair[1]) <= 0):
            # if (verbose):
            #    print("skip chrom \'" + pair[0] + "\' - cannot be mapped to "+to_format)
            continue
        mapping_table[pair[0]] = pair[1]

    return mapping_table

def convert_tsv(mapping_table, text_in_filename, text_out_filename, verbose=False):
    """
    convert chromosome names of a bigwig file according to given mapping_table

    it checks which chromosome names that can correctly mapped, all other chromosomes are skipped
    """
    print("\nstart tsv conversion " + text_in_filename + " --> " + text_out_filename + "...")
    
    matched_chroms = {}
    chrom_field = 0

    if text_in_filename.endswith(".gz"):
        print("zipped file: "+text_in_filename)
        orig = gzip.open(text_in_filename,'rt')
        out = gzip.open(text_out_filename,'wt')
    else:
        orig = open(text_in_filename,'r')
        out = open(text_out_filename,'w')
        
    for line in orig:
        fields = line.split("\t")
        if fields[chrom_field] in mapping_table:
            fields[chrom_field] = mapping_table[fields[chrom_field]]
            out.write("\t".join(fields))
            matched_chroms[fields[chrom_field]] = 1
        elif re.match("^#",fields[chrom_field]) or re.match("^track",fields[chrom_field]):
            out.write("\t".join(fields))
            print("unchanged line: " + fields[chrom_field])
        else:
            if fields[0] not in matched_chroms:
                matched_chroms[fields[chrom_field]] = 0    
                    
    orig.close()
    out.close()
    
    for i in matched_chroms:
        if matched_chroms[i] == 1:
            print(i + " --> OK")
        else:
            print(i + " --> Not mapped")

    if (verbose):
        print("\ntsv conversion " + text_in_filename + " --> " + text_out_filename + " done!")


def convert_bigwig(mapping_table, bw_in_filename, bw_out_filename, verbose=False):
    """
    convert chromosome names of a bigwig file according to given mapping_table

    it checks which chromosome names that can correctly mapped, all other chromosomes are skipped
    """
    print("\nstart bigwig conversion " + bw_in_filename + " --> " + bw_out_filename + "...")
    bw = pyBigWig.open(bw_in_filename)
    curr_chroms = bw.chroms()

    final_mapping_table = {}
    new_chroms = {}

    for c in curr_chroms:
        if c not in mapping_table:
            if (verbose):
                print("skip original chrom \'" + c + "\' - cannot be found in mapping table! Right GENOME & FROM_FORMAT?")
            continue
        final_mapping_table[c] = mapping_table[c]
        new_chroms[mapping_table[c]] = curr_chroms[c]

    if (len(new_chroms) <= 0):
        print("No chromosomes found for mapping! Wrong 'FROM_FORMAT'?")
        print("skip this file...")
        return

    bw_out = pyBigWig.open(bw_out_filename, "w")
    bw_out.addHeader(list(new_chroms.items()))

    for c in final_mapping_table:
        c_int = bw.intervals(c)
        c_map = final_mapping_table[c]
        if verbose:
            print("convert chromosome: ", c, " --> ", c_map)
        bw_out.addEntries(list(itertools.repeat(c_map, len(c_int))), [x[0] for x in c_int], ends=[x[1] for x in c_int], values=[x[2] for x in c_int])

    bw_out.close()
    bw.close()

    if (verbose):
        print("\nbigwig conversion " + bw_in_filename + " --> " + bw_out_filename + " done!")


def main(args=None):

    defaults = {
        'genome': 'GRCm38',
        'fromFormat': 'ensembl',
        'toFormat': 'UCSC',
        'verbose': False,
        'bw_out_filename': None,
        'base_url': 'https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/',
        'gzipped': False
    }

    args = parse_arguments(defaults).parse_args(args)

    mapping_table = get_chromosome_mapping(genome=args.genome, from_format=args.from_format, to_format=args.to_format, verbose=args.verbose, base_url=args.base_url)

    ## https://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python
    textchars = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))

    for curr_file in args.bw_in_filename:
        if args.bw_out_filename is not None and len(args.bw_out_filename) != args.bw_in_filename:
            sys.exit("Please use same number of arguments for --outFileName/-o as well as you have input files!")
        elif args.bw_out_filename is not None and len(args.bw_out_filename) == len(args.bw_in_filename):
            bw_out_filename = args.bw_out_filename.index(curr_file)
        elif curr_file.endswith(".gz"):
            bw_out_filename = re.sub("(.[^.]+.gz)$", ".%s\\1" % (args.to_format + "_chroms"), curr_file)
        else: 
            bw_out_filename = re.sub("(.[^.]+)$", ".%s\\1" % (args.to_format + "_chroms"), curr_file)
        
        if curr_file.endswith(".gz"):
            convert_tsv(mapping_table, curr_file, bw_out_filename, args.verbose)
        elif is_binary_string(open(curr_file, 'rb').read(1024)):
            convert_bigwig(mapping_table, curr_file, bw_out_filename, args.verbose)
        else:
            convert_tsv(mapping_table, curr_file, bw_out_filename, args.verbose)


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)
