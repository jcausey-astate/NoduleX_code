#!/usr/bin/env python
"""
    Collect features from individual (one per nodule) CSV files into
    a matrix with one column per nodule and one row per feature
    (where features with multiple values are flattened into one row
    per value)

    For usage info, run with the --help option.
"""

import sys, os, csv, numpy as np

def main():
    output_file, csv_files, rows, value_delim, heading_delim, padding_value = process_args()
    rows          = check_selected_rows(rows)
    input_matrix  = read_csv_files(csv_files, selected_rows=rows, value_delim=value_delim, heading_delim=heading_delim, padding_value=padding_value)
    feature_wise  = np.array(input_matrix).transpose()
    write_output_matrix(output_file, feature_wise)

def get_usage():
    return __doc__.format(os.path.basename(sys.argv[0])).strip()

def usage():
    print(get_usage)

def process_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('output_file', metavar='output-file', type=str,
                        help='name of output (CSV) file')
    parser.add_argument('input_files', metavar='input-file', type=str,
                        nargs='+', help='name(s) of input file(s) (CSV-like)')
    parser.add_argument('--rows', dest='rows', type=int, nargs='*',
                        help='rows to include in output matrix (numbering from 1)')
    parser.add_argument('--remove-rows', dest='remove_rows', type=int, nargs='*',
                        help='rows to remove from output matrix (numbering from 1)')
    parser.add_argument('--delim', dest='value_delim', type=str, default=',',
                        help='delimiter separating values (for features with more than one)')
    parser.add_argument('--heading-delim', dest='heading_delim', type=str, default=':',
                        help='delimiter separating row headings from values')
    parser.add_argument('--padding', dest='padding_value', type=str, default='0', 
                        help='value to pad uneven features (with more than one but varying # of values)')
    args = parser.parse_args()
    padding_value = float(args.padding_value)
    if int(padding_value) == padding_value:
        padding_value = int(padding_value)
    if args.remove_rows is not None and args.rows is not None:
        print("Error:  Cannot use --remove-rows and --rows simultaneously.")
        parser.usage()
        sys.exit(1)
    rows = args.rows
    if args.remove_rows is not None:
        rows = [-i for i in args.remove_rows]
    return args.output_file, args.input_files, rows, args.value_delim, args.heading_delim, padding_value

def flatten(orig_list):
    flattened = [item for sublist in orig_list for item in sublist]
    return flattened

def fix_length(orig_list, length, padding_value=0):
    fixed_list = orig_list
    if (length is not None) and (len(fixed_list) < length):
        fixed_list.extend([padding_value] * (length - len(fixed_list)))
    return fixed_list

def read_csv_files(file_list, selected_rows='all', value_delim=',', heading_delim=':', padding_value=0):
    matrix        = []
    value_lengths = []
    if not isinstance(selected_rows, list):
        selected_rows = None
    remove_rows   = False
    if selected_rows is not None and any([i < 0 for i in selected_rows]):
        remove_rows = True
    for input_file in file_list:
        with open(input_file, 'rU') as fin:
            csvin = csv.reader(fin, delimiter=heading_delim)
            matrix_row = []
            for idx, row in enumerate(csvin):
                if (selected_rows is None) \
                    or ((not remove_rows) and ( (idx+1)     in selected_rows)) \
                    or (remove_rows       and (-(idx+1) not in selected_rows)):
                    values = flatten([i.split(value_delim) for i in row[1:]])
                    matrix_row.append(values)
                    if idx < len(value_lengths):
                        value_lengths[idx] = max(value_lengths[idx], len(values))
                    else:
                        value_lengths.append(len(values))
            matrix.append(matrix_row)
    for idx, row in enumerate(matrix):
        for i in range(len(row)):
            row[i]  = fix_length(row[i], value_lengths[i], padding_value=padding_value)
        matrix[idx] = flatten(row)
    return matrix

def write_output_matrix(output_filename, matrix, delim=','):
    with open(output_filename, 'w') as fout:
        csvout = csv.writer(fout, delimiter=delim)
        for row in matrix:
            csvout.writerow(row)

def check_selected_rows(rows):
    if rows is not None:
        if any([r < 0 for r in rows]):
            print("Removing rows {0}".format([-i for i in rows]))
            if any([r >= 0 for r in rows]):
                raise(Exception('Cannot use "include" and "exclude" rows simultaneously.'))
    else:
        rows = []
    return rows

if __name__ == '__main__':
    main()
