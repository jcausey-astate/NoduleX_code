#!/usr/bin/env python
"""
    Transpose an entire CSV file so that columns in the input become rows in the output.
    Usage:
        {} input.csv [output.csv]

        input.csv    the input CSV file
        output.csv   the output CSV file (optional); if omitted, output to stdout
"""
from __future__ import print_function
import os, sys, csv

def main():
    if len(sys.argv) < 2:
        print("Error: input CSV required.\n\n{}".format(__doc__.format(os.path.basename(sys.argv[0]))))
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None

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
    for row in zip(*data):
        writer.writerow(row)
    if output_file is not None:
        fout.close()

if __name__ == '__main__':
    main()