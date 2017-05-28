#!/usr/bin/env python
"""
    Re-order the rows in a CSV so that the first column matches the order of rows in the 
    "order-file"; rows not included in the "order-file" are either removed or moved to the
    end (based on -x flag).
    Usage:
        {} [-x] input.csv order-file [output.csv]

        input.csv    the input CSV file
        order-file   list of rows (values from first column) to match in output
        output.csv   the output CSV file (optional); if omitted, output to stdout

        Options:
            -x       remove rows that don't match "order-file"; without this flag, the
                     behavior is to move the non-matching lines to the end.
            -d       use dots as punctuation for matching (R-compatibility)
            -H       input file has a header row; ignore and place it in the output
            -w       warn on missing rows, but continue

"""
from __future__ import print_function
import os, sys, csv, argparse, re

def process_command_line(argv):
    '''Parse the command line and do a first-pass on processing them into a
    format appropriate for the rest of the script.'''

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Re-order the rows in a CSV so that the first column matches the "
                                              "order of rows\nin the  \"order-file\"; rows not included in the "
                                              "\"order-file\" are either removed\nor moved to the end "
                                              "(based on -x flag).")
    parser.add_argument(
        "-x", dest='remove_extras', action="store_true",
        help="Remove non-matching rows.")
    parser.add_argument(
        "-d", dest='use_dots', action="store_true",
        help="Replace punctuation in \"order-file\" values with dots (for R-compatibility)")
    parser.add_argument(
        "-H", dest='header_row', action='store_true',
        help="Input file has a header row; ignore first row and copy it to the output.")
    parser.add_argument(
        "-w", dest='warn_missing', action="store_true",
        help="Warn about missing rows, but continue.")
    parser.add_argument(
        "input_file", metavar='input.csv',
        help="The input CSV file name.")
    parser.add_argument(
        "order_file", metavar='order-file',
        help="Specification of order for output rows.")
    parser.add_argument(
        "output_file", metavar='output.csv', nargs='?',
        help="Output CSV file name.")
    args = parser.parse_args(argv[1:])
    return args

def main():
    args = process_command_line(sys.argv)
    input_file  = args.input_file
    output_file = args.output_file
    order_file  = args.order_file

    order_spec  = []
    with open(order_file, 'rU') as fin:
        for line in fin:
            line = line.strip()
            order_spec.append(line)

    if args.use_dots:
        order_spec = [re.sub(r'(\W|_)', '.', s) for s in order_spec]

    data    = []
    dialect = None
    with open(input_file, 'rU') as fin:
        try:
            dialect = csv.Sniffer().sniff(fin.read(1024*100))
        except:
            dialect = csv.excel
        fin.seek(0)
        if dialect.delimiter.isalpha():  # Sniffer finds 'I' sometimes; it should default to ',' if this happens.
            dialect.delimiter = ','
        reader = csv.reader(fin, dialect=dialect)
        for line in reader:
            data.append(line)
    fout = sys.stdout
    if output_file is not None:
        fout = open(output_file, 'wb')
    writer = csv.writer(fout, dialect=dialect)

    if args.header_row:                 # if there is a header, just copy it to the output file and remove from `data`
        writer.writerow(data[0])
        data = data[1:]
    
    had_error = False
    for next_select in order_spec:
        did_write   = False
        if args.use_dots:
                next_select = re.sub(r'(\W|_)', '.', next_select)
        for idx, row in enumerate(data):
            header = row[0].strip()
            if args.use_dots:
                header = re.sub(r'(\W|_)', '.', header)
                row[0] = header
            if header == next_select:
                writer.writerow(row)
                did_write = True
                del data[idx]
                break
        if not did_write:
            print("{}: required output row \"{}\" missing in input file.".format('Error' if not args.warn_missing else 'Warning', next_select), file=sys.stderr)
            had_error = True
            if not args.warn_missing:
                break
    # If there are more rows, output them now, unless the -x flag was used.
    continue_output = (args.warn_missing or not had_error)
    if not args.remove_extras and continue_output:
        for row in data:
            if args.use_dots:
                row[0] = re.sub(r'(\W|_)', '.', row[0])
            writer.writerow(row)
    if output_file is not None:
        fout.close()
        if not args.warn_missing and had_error:
            os.unlink(output_file)

if __name__ == '__main__':
    main()
