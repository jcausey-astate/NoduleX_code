label_f = importdata( label_list );
unseg_f = importdata( unseg_list );

for i = 1:10
    
    file_name_bin  = label_f{i}
    file_name_gray = unseg_f{i}
    
  segmentation_num = 01 ;
  kernel = 'B35f' ;
  gray_threshold = -999 ;
  fast_only = 0 ;
       
	[metrics_3d, metrics_2d] = lna(file_name_bin , ...
                                   file_name_gray, ...
                                   segmentation_num, ...
                                   kernel, ...
                                   gray_threshold, ...
                                   fast_only )
    disp( metrics_3d )
    disp( metrics_2d )
    
end
