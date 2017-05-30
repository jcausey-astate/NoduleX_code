function results_to_csv(input_matfile, output_prefix, separator=',')
    results = load(input_matfile);
    field_name = fieldnames(results){1};
    cellvalues = results.(field_name);
    for i = 1:columns(cellvalues)
        cell2csv([output_prefix, sprintf('_part_%05d', i), '.csv'], cellvalues{i}, separator);
    end
endfunction