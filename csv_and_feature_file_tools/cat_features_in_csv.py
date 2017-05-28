#!/usr/bin/env python
"""
    Concatenate additional features into a matrix (CSV) where the rows represent
    "samples" (patients, nodules, etc.) and the columns represent "feature" values.
    The "new-features.csv" file must contain a row corresponding to every row in the
    "input.csv" file.
"""
from __future__ import print_function
"""
    Usage:
        {} [-h] input.csv new-features.csv [output.csv]

        input.csv          the input CSV file
        new-features.csv   file containing news features to concatenate with "input.csv"
        output.csv         the output CSV file (optional); if omitted, output to stdout

        Options:
            -H       "new-features.csv" incudes a "header" column that must be ignored
"""
import os, sys, csv, argparse

def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=__doc__)
    parser.add_argument(
        "-H", dest='header_col', action="store_true",
        help="\"new-features.csv\" incudes a \"header\" column that must be ignored")
    parser.add_argument(
        "input_file", metavar='input.csv',
        help="The input CSV file name.")
    parser.add_argument(
        "new_features_file", metavar='new-features.csv',
        help="file containing news features to concatenate with \"input.csv\"")
    parser.add_argument(
        "output_file", metavar='output.csv', nargs='?',
        help="Output CSV file name. (output to stdout if omitted)")
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = process_command_line(sys.argv)
    input_file         = args.input_file
    output_file        = args.output_file
    new_features_file  = args.new_features_file

    dialect = None
    with open(input_file, 'rU') as fin:
        try:
            dialect = csv.Sniffer().sniff(fin.read(1024*100))
        except:
            dialect = csv.excel
        if dialect.delimiter.isalpha():  # Sniffer finds 'I' sometimes; it should default to ',' if this happens.
            dialect.delimiter = ','

    fout = sys.stdout
    if output_file is not None:
        fout = open(output_file, 'wb')
    writer = csv.writer(fout, dialect=dialect)

    write_error = False
    with open(input_file, 'rU') as fin, open(new_features_file, 'rU') as nfin:
        n_dialect = None
        try:
            n_dialect = csv.Sniffer().sniff(nfin.read(1024*100))
        except:
            n_dialect = csv.excel
        if n_dialect.delimiter.isalpha():  # Sniffer finds 'I' sometimes; it should default to ',' if this happens.
            n_dialect.delimiter = ','
        nfin.seek(0)
        i_reader = csv.reader(fin, dialect=dialect)
        n_reader = csv.reader(nfin, dialect=n_dialect)
        # For each line in the input file, concatenate the corresponding line from the "new-features" file and write.
        for row in i_reader:
            n_row = None
            try:
                n_row = n_reader.next()
            except StopIteration:
                print("Error: Number of rows in '{}' and '{}' do not match.".format(input_file, new_features_file), file=sys.stderr)
                write_error = True
                break
            if args.header_col:
                n_row = n_row[1:]
            out_row = row
            out_row.extend(n_row)
            writer.writerow(out_row)

    if output_file is not None:
        fout.close()
        if write_error:
            os.unlink(output_file)
            sys.exit(1)

if __name__ == '__main__':
    main()
