#!/usr/bin/env python
"""
    Adds header row, column, or both to a CSV file given files, using
    files for the source of the headers.
    Fails if the number of header values given doesn't match the number
    of rows/columns they are being applied to.
"""
from __future__ import print_function
import os, sys, csv, argparse, re

def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description=__doc__)
    
    parser.add_argument(
        "-r", "--row", dest='row_file', metavar='row-headers-file',
        help="File containing one entry per line for each row in input-csv.")
    parser.add_argument(
        "-c", "--col", dest='col_file', metavar='col-headers-file',
        help="File containing one entry per line for each column in input-csv.")
    parser.add_argument(
        "input_file", metavar='input-csv',
        help="The input CSV file name.")
    parser.add_argument(
        "output_file", metavar='output.csv', nargs='?',
        help="Output CSV file name. If omitted, the input file will be changed in place.")
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = process_command_line(sys.argv)
    input_file  = args.input_file
    output_file = args.output_file
    col_file    = args.col_file
    row_file    = args.row_file

    row_headings = None
    col_headings = None
    if row_file is not None:
        row_headings = [s.strip() for s in open(row_file, 'rU').readlines()]
    if col_file is not None:
        col_headings = [s.strip() for s in open(col_file, 'rU').readlines()]

    data    = []
    dialect = None
    n_cols  = None
    with open(input_file, 'rU') as fin:
        try:
            dialect = csv.Sniffer().sniff(fin.read(1024*100))
        except:
            dialect = csv.excel
        fin.seek(0)
        if dialect.delimiter.isalnum() or dialect.delimiter == '.': # Sniffer finds 'I' or '0' a decimal point sometimes; it should default to ',' if this happens.
            dialect.delimiter = ','
        reader = csv.reader(fin, dialect=dialect)
        if col_headings is not None:
            if len(col_headings) != len(reader[0]):
                print("Error: Number of column headings provided did not match number of columns in the CSV file.", file=sys.stderr)
                sys.exit(1)
            else:
                data.append(col_headings)
        for idx, line in enumerate(reader):
            if n_cols is None:
                n_cols = len(line)
            if row_headings is not None:
                if len(row_headings) <= idx:
                    print("Error: Number of row headings provided did not match number of rows in the CSV file.", file=sys.stderr)
                    sys.exit(1)
                line.insert(0, row_headings[idx])
            data.append(line)
    fout = sys.stdout
    if output_file is None:
        output_file = input_file
    if output_file != '':
        fout = open(output_file, 'wb')
    if n_cols == 1:
        dialect.delimiter = '\t'
        dialect.quoting   = csv.QUOTE_MINIMAL
    writer = csv.writer(fout, dialect=dialect)
    for idx, row in enumerate(data):
        writer.writerow(row)

if __name__ == '__main__':
    main()
