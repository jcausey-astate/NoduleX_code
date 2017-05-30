%**************************************************************************
% Process the acrylic cylinder. This is a test of version 1.1 of the Lung
% Nodule Analysis code to see if the output can be imported into the REDCAP
% database.
%**************************************************************************

% Process the patients that have been segmented by Slicer as of 2012-02-28.

file_name_bin  = '/media/USB20FD/Slicer_Output/001-FGR/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/001-FGR_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/001-FGR/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/001-FGR_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/003-DLR/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/003-DLR_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/003-DLR/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/003-DLR_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/003-DLR/NODULE_SCAN_0_6_B35F_0001/Segmentation_02/003-DLR_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/003-DLR/NODULE_SCAN_0_6_B35F_0001/Segmentation_02/003-DLR_B35F_cropped' ;
segmentation_num = 02 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/005-D-C/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/005-D-C_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/005-D-C/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/005-D-C_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)



file_name_bin  = '/media/USB20FD/Slicer_Output/006a-RGO/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/006a-RGO_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/006a-RGO/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/006a-RGO_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/008-DRQ/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/008-DRQ_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/008-DRQ/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/008-DRQ_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/010-TBS/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/010-TBS_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/010-TBS/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/010-TBS_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/011-H-S/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/011-H-S_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/011-H-S/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/011-H-S_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/013-JKD/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/013-JKD_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/013-JKD/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/013-JKD_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/014-JCB/NODULE_SCAN_0_6_B35F_0002/Segmentation_01/014-JCB_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/014-JCB/NODULE_SCAN_0_6_B35F_0002/Segmentation_01/014-JCB_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/015-MAC/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/015-MAC_B35F_cropped-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/015-MAC/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/015-MAC_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/017-WHF/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/017-WHF_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/017-WHF/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/017-WHF_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/021-MSH/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/021-MSH_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/021-MSH/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/021-MSH_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)


file_name_bin  = '/media/USB20FD/Slicer_Output/024-GRO/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/024_GRO_B35F_cropped1-1-label'
file_name_gray = '/media/USB20FD/Slicer_Output/024-GRO/NODULE_SCAN_0_6_B35F_0001/Segmentation_01/024_GRO_B35F_cropped' ;
segmentation_num = 01 ;
kernel = 'B35f' ;
gray_threshold = -999 ;
fast_only = 0 ;
lna01_4(file_name_bin, file_name_gray, segmentation_num, kernel, ...
    gray_threshold, fast_only)

