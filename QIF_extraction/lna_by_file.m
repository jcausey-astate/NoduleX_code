function lna_by_file(binary_file, gray_file, patient, output_dir)

    addpath('Matlab_Source') ;
    addpath('Matlab_Source/boxcount') ;
    pkg load image


    file_name_bin   = binary_file ;
    file_name_gray  = gray_file ;
    metrics_2d_file = [output_dir '/' patient '_features2d'] ;
    metrics_3d_file = [output_dir '/' patient '_features3d'] ;

    segmentation_num = 01 ;
    kernel = 'B35f' ;
    gray_threshold = -999 ;
    % gray_threshold = -500 ;
    fast_only = 1 ;

    disp(patient)
                   
    [metrics_3d, metrics_2d] = lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, gray_threshold, fast_only ) ;

    save([metrics_2d_file '.hd5'], '-hdf5',  'metrics_2d');
    save([metrics_3d_file '.hd5'], '-hdf5',  'metrics_3d');

end
