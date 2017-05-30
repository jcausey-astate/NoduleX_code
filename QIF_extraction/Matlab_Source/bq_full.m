function [] = bq_full(filename_segm, filename_gray, bone_num_vec)

% bq_full.m (Calculates the 3-D distance to surface metric that takes so
% long.)
%
% Calculate metrics for some lung nodules in support of an ICTS grant
% proposal due on Monday, November 16, 2009.
%
% Metrics have been proposed for 3-D image volumes and for 2-D image
% slices.

% Read in the binary-threshholded 3-D image with the nodule only, which is
% unsigned 8-bit.

% Also read in the grayscale 3-D image which contains the nodule. It is
% signed 16-bit.
%
% Software has been "re-purposed" to analyze bones rather than lung
% nodules. DGP 2010-04-29

disp(' ') ;
disp(' ') ;
disp(' ') ;
% filename_segm = input('Filename of binary segmentation without extension:    ', 's') ;
% filename_gray = input('Filename of grayscale image without extension:        ', 's') ;
% bone_num_vec  = input('Bone number(s) to be processed (vector of integers):  ') ;

for jjj = 1:length(bone_num_vec)
    
    bone_index = bone_num_vec(jjj) ;
    
    disp_string = ['Working on bone number: ' num2str(bone_index)] ;
    disp(disp_string) ;

    diary_name = [filename_segm '_' num2str(bone_index, '%02d') '.diary'] ;
    diary(diary_name) ;

    header_name = [filename_segm '.hdr'] ;
    header = analyze75info(header_name) ;
    dy = double(header.PixelDimensions(1)) ;
    dx = double(header.PixelDimensions(2)) ;
    dz = double(header.PixelDimensions(3)) ;
    image_segm_name = [filename_segm '.img'] ;
    a = int8(analyze75read(image_segm_name)) ;
    temp = size(a) ;
    width = temp(2) ;
    height = temp(1) ;
    depth = temp(3) ;
    s_3d_orig = double(a) ;
    clear a temp ;
    disp_string = ['# rows = ' num2str(height) '   # columns = ' num2str(width) '   # planes = ' num2str(depth)] ;
    disp(disp_string) ;

    image_gray_name = [filename_gray '.img'] ;
    g = analyze75read(image_gray_name) ;
    g_3d = double(g) ;
    clear g ;
    
    s_3d = zeros(size(s_3d_orig)) ;
    qqq = find(s_3d_orig == bone_index) ;
    s_3d(qqq) = 1 ;

    %**************************************************************************
    %**************************************************************************
    %**************************************************************************
    % CALCULATE 3-D METRICS
    %**************************************************************************
    %**************************************************************************
    %**************************************************************************

    disp('CALCULATING 3-D METRICS ...') ;

    %**************************************************************************
    % Calculate size metrics.
    %**************************************************************************

    disp('   Calculating size metrics ...') ;

    % Calculate the volume.

    disp('      Calculating volume ...') ;
    number_voxels_inside_3d = sum(sum(sum(s_3d))) ;
    volume_3d = number_voxels_inside_3d * dx * dy * dz ;
    disp('      Finished calculating volume.') ;

    % Calculate the surface area.

    disp('      Calculating surface area ...') ;
    type_a_surfaces_3d = sum(sum(sum(abs(s_3d(:,:,1:(end-1)) - s_3d(:,:,2:end))))) ;
    type_b_surfaces_3d = sum(sum(sum(abs(s_3d(:,1:(end-1),:) - s_3d(:,2:end,:))))) ;
    type_c_surfaces_3d = sum(sum(sum(abs(s_3d(1:(end-1),:,:) - s_3d(2:end,:,:))))) ;
    total_surfaces_3d  = type_a_surfaces_3d + type_b_surfaces_3d + type_c_surfaces_3d ;

    surface_area_3d = (type_a_surfaces_3d * dx * dy) ...
                    + (type_b_surfaces_3d * dx * dz) ...
                    + (type_c_surfaces_3d * dy * dz) ;
    disp('      Finished calculating surface area.') ;

    disp('   Finished calculating size metrics.') ;

    %**************************************************************************
    % Calculate shape metrics.
    %**************************************************************************

    disp('   Calculating shape metrics ...') ;

    % Calculate the sphericity.

    disp('      Calculating sphericity ...') ;
    sphericity_3d = (pi^(1/3)) * ((6 * volume_3d)^(2/3)) / surface_area_3d ;
    disp('      Finished calculating sphericity.') ;

    % Calculate the compactness.

    disp('      Calculating compactness ...') ;
    compactness_3d = (36 * pi * (volume_3d^2)) / (surface_area_3d^3) ;
    disp('      Finished calculating compactness.') ;

    % Calculate the principal rotational moments.
    % This calculation uses the binary segmentation image.
    % See http://en.wikipedia.org/Moment_of_inertia for a definition of the
    % tensor matrix that is calculated.
    % The tensor matrix is then diagonalized using a singular-value
    % decomposition (SVD) and normalized to have unit energy so that it is
    % independent of size.

    disp('      Calculating principal rotational moments ...') ;

    % Find the centroid. This will be a tuple.

    x = 1:width ;
    y = 1:height ;
    z = 1:depth ;
    [X, Y, Z] = meshgrid(x, y, z) ;
    x_c = sum(sum(sum(s_3d .* X))) / number_voxels_inside_3d ;
    y_c = sum(sum(sum(s_3d .* Y))) / number_voxels_inside_3d ;
    z_c = sum(sum(sum(s_3d .* Z))) / number_voxels_inside_3d ;
    clear X Y Z ;

    x_shift = x - x_c ;
    y_shift = y - y_c ;
    z_shift = z - z_c ;
    [X_shift, Y_shift, Z_shift] = meshgrid(x_shift, y_shift, z_shift) ;
    clear x y z x_c y_c z_c x_shift y_shift z_shift ;

    x_sq_term_3d = sum(sum(sum((s_3d .* X_shift).^2))) / number_voxels_inside_3d ;
    y_sq_term_3d = sum(sum(sum((s_3d .* Y_shift).^2))) / number_voxels_inside_3d ;
    z_sq_term_3d = sum(sum(sum((s_3d .* Z_shift).^2))) / number_voxels_inside_3d ;

    xy_term_3d   = sum(sum(sum((s_3d .* X_shift) .* (s_3d .* Y_shift)))) / number_voxels_inside_3d ;
    xz_term_3d   = sum(sum(sum((s_3d .* X_shift) .* (s_3d .* Z_shift)))) / number_voxels_inside_3d ;
    yz_term_3d   = sum(sum(sum((s_3d .* Y_shift) .* (s_3d .* Z_shift)))) / number_voxels_inside_3d ;

    clear X_shift Y_shift Z_shift ;

    % Define the tensor as shown on the wikipedia web page.

    I11_3d = y_sq_term_3d + z_sq_term_3d ;
    I22_3d = x_sq_term_3d + z_sq_term_3d ;
    I33_3d = x_sq_term_3d + y_sq_term_3d ;

    I12_3d = xy_term_3d ;
    I21_3d = I12_3d ;

    I13_3d = xz_term_3d ;
    I31_3d = I13_3d ;

    I23_3d = yz_term_3d ;
    I32_3d = I23_3d ;

    I_matrix_3d = [ I11_3d -I12_3d -I13_3d ; ...
                   -I21_3d  I22_3d -I23_3d ; ...
                   -I31_3d -I32_3d  I33_3d] ;

    % Diagonalize it. The r_i_3d are the answers we seek.

    [U, S, V] = svd(I_matrix_3d) ;
    norm_fac = sqrt(sum(diag(S) .^ 2)) ;
    r_1_3d = S(1,1) / norm_fac ;
    r_2_3d = S(2,2) / norm_fac ;
    r_3_3d = S(3,3) / norm_fac ;
    clear norm_fac U S V ;
    clear I11_3d I22_3d I33_3d I12_3d I21_3d I13_3d I31_3d I23_3d I32_3d ;

    disp('      Finished calculating principal rotational moments.') ;

    disp('   Finished calculating shape metrics.') ;

    %**************************************************************************
    % Calculate texture metrics.
    %**************************************************************************

    disp('   Calculating texture metrics ...') ;

    % Calculate statistics of the values of the voxels.

    disp('      Calculating statistics of gray levels ...') ;
    i = find(s_3d == 1) ;
    mean_gray_3d     = mean(g_3d(i)) ;
    std_gray_3d      = std(g_3d(i))  ;
    var_gray_3d      = var(g_3d(i))  ;
    skewness_gray_3d = skewness(g_3d(i)) ;
    kurtosis_gray_3d = kurtosis(g_3d(i)) ;
    disp('      Finished calculating statistics of gray levels.') ;

    % Now do the same thing after differencing with the nearest neighbors.

    disp('      Calculating statistics of difference images of gray levels ...') ;
    kern = zeros(3, 3, 3) ;
    q = -1/6 ;
    kern(:,:,1) = [0 0 0 ; ...
                   0 q 0 ; ...
                   0 0 0] ;
    kern(:,:,2) = [0 q 0 ; ...
                   q 1 q ; ...
                   0 q 0] ;
    kern(:,:,3) = [0 0 0 ; ...
                   0 q 0 ; ...
                   0 0 0] ;
    g_3d_diff = convn(g_3d, kern, 'same') ;
    mean_gray_diff_3d     = mean(g_3d_diff(i)) ;
    std_gray_diff_3d      = std(g_3d_diff(i)) ;
    var_gray_diff_3d      = var(g_3d_diff(i))  ;
    skewness_gray_diff_3d = skewness(g_3d_diff(i)) ;
    kurtosis_gray_diff_3d = kurtosis(g_3d_diff(i)) ;
    clear i q kern ;
    disp('      Finished calculating statistics of difference images of gray levels.') ;

    disp('   Finished calculating texture metrics.') ;

    %**************************************************************************
    % Calculate margin metrics.
    %**************************************************************************

    % DO NOT COMMENT OUT THE FOLLOWING 4 LINES!
    summed_distance_3d      = 0 ;
    mean_distance_3d        = 0 ;
    norm_summed_distance_3d = 0 ;
    norm_mean_distance_3d   = 0 ;
    % DO NOT COMMENT OUT THE PRECEDING 4 LINES!

    disp('   Calculating margin metrics ...') ;

    % Calculate the mean distance to surface.

    disp('      Calculating distance to surface ...') ;
    dist_img_name = [filename_segm '_3d_' num2str(bone_index, '%02d') '.dst'] ;
    fid = fopen(dist_img_name, 'w') ;

    disp('         Finding border voxels ...') ;
    kern = ones(3,3,3) ;
    border = sign((1 - s_3d) .* convn(s_3d, kern, 'same')) ;

    [iii_list, jjj_list, kkk_list] = ind2sub(size(s_3d), find(border > 0)) ;
    clear border ;
    disp('         Finished finding border voxels.') ;

    disp('         Allocating 3D image for distance to surface ...') ;
    dist_exterior_3d = single(zeros(height, width)) ;
    disp('         Finished allocating 3D image for distance to surface.') ;

    summed_distance_3d = 0.0 ;

    for k = 1:depth
        disp_string = ['         Plane ' num2str(k) ' out of ' num2str(depth) '. ' datestr(now())] ;
        disp(disp_string) ;
        term_3 = ((k - kkk_list) * dz).^2 ;
        for j = 1:width
            term_2 = ((j - jjj_list) * dx).^2 ;
            for i = 1:height
                if (s_3d(i,j,k) == 1)
                    term_1 = ((i - iii_list) * dy).^2 ;
                    dist_sq_vec = term_1 + term_2 + term_3 ;
                    dist_exterior_3d(i, j) = sqrt(min(dist_sq_vec)) ;
                    summed_distance_3d = summed_distance_3d + dist_exterior_3d(i, j) ;
                end
            end
        end
        count = fwrite(fid, dist_exterior_3d, 'single') ;
    end
    mean_distance_3d = summed_distance_3d / number_voxels_inside_3d ;
    norm_summed_distance_3d = summed_distance_3d / (volume_3d^(1/3)) ;
    norm_mean_distance_3d = mean_distance_3d / (volume_3d^(1/3)) ;

    fclose(fid) ;
    disp('      Finished calculating distance to surface.') ;

    disp('   Finished calculating margin metrics.') ;

    disp('FINISHED CALCULATING 3-D METRICS.') ;




    %**************************************************************************
    %**************************************************************************
    %**************************************************************************
    % CALCULATE 2-D METRICS
    %**************************************************************************
    %**************************************************************************
    %**************************************************************************

    disp('CALCULATING 2-D METRICS ...') ;

    max_inside_2d = 0 ;
    plane_index_of_maximal_area = 0 ;
    for i = 1:depth
       test = sum(sum(s_3d(:,:,i))) ;
       if (test > max_inside_2d)
          max_inside_2d = test ;
          plane_index_of_maximal_area = i ;
       end
    end
    s_2d = s_3d(:,:,plane_index_of_maximal_area) ;
    g_2d = g_3d(:,:,plane_index_of_maximal_area) ;
    clear max_inside_2d index_best i test ;

    %**************************************************************************
    % Calculate size metrics ...') ;
    %**************************************************************************

    disp('   Calculating size metrics ...') ;

    % Calculate the area.

    disp('      Calculating area ...') ;
    number_pixels_inside_2d = sum(sum(s_2d)) ;
    area_2d = number_pixels_inside_2d * dx * dy ;
    disp('      Finished calculating area.') ;

    % Calculate the perimeter.

    disp('      Calculating perimeter ...') ;
    type_b_edges_2d = sum(sum(abs(s_2d(:,1:(end-1)) - s_2d(:,2:end)))) ;
    type_c_edges_2d = sum(sum(abs(s_2d(1:(end-1),:) - s_2d(2:end,:)))) ;
    total_edges_2d  = type_b_edges_2d + type_c_edges_2d ;

    perimeter_2d = (type_b_edges_2d * dx) ...
                 + (type_c_edges_2d * dy) ;
    disp('      Finished calculating perimeter.') ;

    disp('   Finished calculating size metrics.') ;

    %**************************************************************************
    % Calculate shape metrics.
    %**************************************************************************

    disp('   Calculating shape metrics ...') ;

    % Calculate the circularity (McNitt-Gray's definition).

    disp('      Calculating circularity ...') ;
    circularity_2d = (4 * pi * area_2d) / (perimeter_2d^2) ;
    disp('      Finished calculating circularity.') ;

    % Calculate the principal rotational moments.
    % This calculation uses the binary segmentation image.
    % See http://en.wikipedia.org/Moment_of_inertia for a definition of the
    % tensor matrix that is calculated (modified for 2D).
    % The tensor matrix is then diagonalized using a singular-value
    % decomposition (SVD) and normalized to have unit energy so that it is
    % independent of size.

    disp('      Calculating principal rotational moments ...') ;

    % Find the centroid. This will be a tuple.

    x = 1:width ;
    y = 1:height ;
    [X, Y] = meshgrid(x, y) ;
    x_c = sum(sum(s_2d .* X)) / number_pixels_inside_2d ;
    y_c = sum(sum(s_2d .* Y)) / number_pixels_inside_2d ;
    clear X Y ;

    x_shift = x - x_c ;
    y_shift = y - y_c ;
    [X_shift, Y_shift] = meshgrid(x_shift, y_shift) ;
    clear x y x_c y_c x_shift y_shift ;

    x_sq_term_2d = sum(sum((s_2d .* X_shift).^2)) / number_pixels_inside_2d ;
    y_sq_term_2d = sum(sum((s_2d .* Y_shift).^2)) / number_pixels_inside_2d ;

    xy_term_2d   = sum(sum((s_2d .* X_shift) .* (s_2d .* Y_shift))) / number_pixels_inside_2d ;

    clear X_shift Y_shift Z_shift ;

    % Define the tensor as shown on the wikipedia web page, but modified for
    % 2D.

    I11_2d = y_sq_term_2d ;
    I22_2d = x_sq_term_2d ;

    I12_2d = xy_term_2d ;
    I21_2d = I12_2d ;

    I_matrix_2d = [ I11_2d -I12_2d ; ...
                   -I21_2d  I22_2d] ;

    % Diagonalize it. The r_i_2d are the answers we seek.

    [U, S, V] = svd(I_matrix_2d) ;
    norm_fac = sqrt(sum(diag(S) .^ 2)) ;
    r_1_2d = S(1,1) / norm_fac ;
    r_2_2d = S(2,2) / norm_fac ;
    clear norm_fac U S V ;
    clear I11_2d I22_2d I12_2d I21_2d ;

    disp('      Finished calculating principal rotational moments.') ;

    disp('   Finished calculating shape metrics.') ;

    %**************************************************************************
    % Calculate texture metrics.
    %**************************************************************************

    disp('   Calculating texture metrics ...') ;

    % Calculate statistics of the values of the voxels.

    disp('      Calculating statistics of gray levels ...') ;
    i = find(s_2d == 1) ;
    mean_gray_2d     = mean(g_2d(i)) ;
    std_gray_2d      = std(g_2d(i))  ;
    var_gray_2d      = var(g_2d(i))  ;
    skewness_gray_2d = skewness(g_2d(i)) ;
    kurtosis_gray_2d = kurtosis(g_2d(i)) ;
    disp('      Finished calculating statistics of gray levels.') ;

    % Now do the same thing after differencing with the nearest neighbors.

    disp('      Calculating statistics of difference images of gray levels ...') ;
    kern = zeros(3, 3) ;
    q = -1/4 ;
    kern = [0 q 0 ; ...
            q 1 q ; ...
            0 q 0] ;
    g_2d_diff = convn(g_2d, kern, 'same') ;
    mean_gray_diff_2d     = mean(g_2d_diff(i)) ;
    std_gray_diff_2d      = std(g_2d_diff(i)) ;
    var_gray_diff_2d      = var(g_2d_diff(i)) ;
    skewness_gray_diff_2d = skewness(g_2d_diff(i)) ;
    kurtosis_gray_diff_2d = kurtosis(g_2d_diff(i)) ;
    clear i q kern ;
    disp('      Finished calculating statistics of difference images of gray levels.') ;

    disp('   Finished calculating texture metrics.') ;

    %**************************************************************************
    % Calculate margin metrics.
    %**************************************************************************

    disp('   Calculating margin metrics ...') ;

    % Calculate the mean distance to surface.

    disp('      Calculating distance to surface ...') ;
    dist_img_name = [filename_segm '_2d_' num2str(bone_index, '%02d') '.dst'] ;
    fid = fopen(dist_img_name, 'w') ;

    disp('         Finding border voxels ...') ;
    kern = ones(3,3) ;
    border = sign((1 - s_2d) .* convn(s_2d, kern, 'same')) ;

    [iii_list, jjj_list] = ind2sub(size(s_2d), find(border > 0)) ;
    clear border ;
    disp('         Finished finding border voxels.') ;

    disp('        Allocating 2D image for distance to surface ...') ;
    dist_exterior_2d = single(zeros(height, width)) ;
    disp('        Finished allocating 2D image for distance to surface.') ;

    summed_distance_2d = 0.0 ;

    for j = 1:width
       term_2 = ((j - jjj_list) * dx).^2 ;
       for i = 1:height
          if (s_2d(i,j) == 1)
             term_1 = ((i - iii_list) * dy).^2 ;
             dist_sq_vec = term_1 + term_2 ;
             dist_exterior_2d(i, j) = sqrt(min(dist_sq_vec)) ;
             summed_distance_2d = summed_distance_2d + dist_exterior_2d(i, j) ;
          end
       end
    end
    count = fwrite(fid, dist_exterior_2d, 'single') ;
    mean_distance_2d = summed_distance_2d / number_pixels_inside_2d ;
    norm_summed_distance_2d = summed_distance_2d / (area_2d^(1/2)) ;
    norm_mean_distance_2d = mean_distance_2d / (area_2d^(1/2)) ;

    fclose(fid) ;
    disp('      Finished calculating distance to surface.') ;

    disp('   Finished calculating margin metrics.') ;

    disp('FINISHED CALCULATING 2-D METRICS.') ;

    % Save the results.

    mat_name = [filename_segm '.mat'] ;
    eval_string = ['save ''' filename_segm ''' height width depth number_voxels_inside_3d volume_3d surface_area_3d sphericity_3d compactness_3d r_1_3d r_2_3d r_3_3d mean_gray_3d std_gray_3d var_gray_3d skewness_gray_3d kurtosis_gray_3d mean_gray_diff_3d std_gray_diff_3d var_gray_diff_3d skewness_gray_diff_3d kurtosis_gray_diff_3d summed_distance_3d mean_distance_3d norm_summed_distance_3d norm_mean_distance_3d plane_index_of_maximal_area number_pixels_inside_2d area_2d perimeter_2d circularity_2d r_1_2d r_2_2d mean_gray_2d std_gray_2d var_gray_2d skewness_gray_2d kurtosis_gray_2d mean_gray_diff_2d std_gray_diff_2d var_gray_diff_2d skewness_gray_diff_2d kurtosis_gray_diff_2d'] ;
    eval(eval_string) ;

    clear disp_string eval_string ;

    % Print the results on the screen.

    format long ;
    disp(' ') ;
    disp(['Filename = ' filename_segm]) ;
    disp(['     Height of input image (# rows) = ' num2str(height)]) ;
    disp(['   Width of input image (# columns) = ' num2str(width)]) ;
    disp(['    Depth of input image (# planes) = ' num2str(depth)]) ;
    disp([' ']) ;
    disp(['                         Bone Index = ' num2str(bone_index)]) ;
    disp(['****************************** STATISTICS FOR 3D ******************************']) ;
    disp(['********** Size Metrics **********']) ;
    disp(['                                      Number of voxels inside = ' num2str(number_voxels_inside_3d)]) ;
    disp(['                                                       Volume = ' num2str(volume_3d)]) ;
    disp(['                                                 Surface Area = ' num2str(surface_area_3d)]) ;
    disp(['********** Shape Metrics **********']) ; 
    disp(['                                                   Sphericity = ' num2str(sphericity_3d)]) ;
    disp(['                                                  Compactness = ' num2str(compactness_3d)]) ;
    disp(['                                       Rotational moment, r_1 = ' num2str(r_1_3d)]) ;
    disp(['                                       Rotational moment, r_2 = ' num2str(r_2_3d)]) ;
    disp(['                                       Rotational moment, r_3 = ' num2str(r_3_3d)]) ;
    disp(['********** Texture Metrics **********']) ;
    disp(['                                     Mean of HU inside nodule = ' num2str(mean_gray_3d)]) ;
    disp(['                       Standard Deviation of HU inside nodule = ' num2str(std_gray_3d)]) ;
    disp(['                                 Variance of HU inside nodule = ' num2str(var_gray_3d)]) ;
    disp(['                                 Skewness of HU inside nodule = ' num2str(skewness_gray_3d)]) ;
    disp(['                                 Kurtosis of HU inside nodule = ' num2str(kurtosis_gray_3d)]) ;
    disp(['                 Mean of difference image of HU inside nodule = ' num2str(mean_gray_diff_3d)]) ;
    disp(['   Standard Deviation of difference image of HU inside nodule = ' num2str(std_gray_diff_3d)]) ;
    disp(['             Variance of difference image of HU inside nodule = ' num2str(var_gray_diff_3d)]) ;
    disp(['             Skewness of difference image of HU inside nodule = ' num2str(skewness_gray_diff_3d)]) ;
    disp(['             Kurtosis of difference image of HU inside nodule = ' num2str(kurtosis_gray_diff_3d)]) ;
    disp(['********** Margin Metrics **********']) ;
    disp(['                                              Summed distance = ' num2str(summed_distance_3d)]) ;
    disp(['                                                Mean distance = ' num2str(mean_distance_3d)]) ;
    disp(['                                   Normalized summed distance = ' num2str(norm_summed_distance_3d)]) ;
    disp(['                                     Normalized mean distance = ' num2str(norm_mean_distance_3d)]) ;
    disp(['****************************** STATISTICS FOR 2D ******************************']) ;
    disp(['                           Transverse plane with maximal area = ' num2str(plane_index_of_maximal_area)]) ;
    disp(['********** Size Metrics **********']) ;
    disp(['                                      Number of pixels inside = ' num2str(number_pixels_inside_2d)]) ;
    disp(['                                                         Area = ' num2str(area_2d)]) ;
    disp(['                                                    Perimeter = ' num2str(perimeter_2d)]) ;
    disp(['********** Shape Metrics **********']) ;
    disp(['                                                  Circularity = ' num2str(circularity_2d)]) ;
    disp(['                                       Rotational moment, r_1 = ' num2str(r_1_2d)]) ;
    disp(['                                       Rotational moment, r_2 = ' num2str(r_2_2d)]) ;
    disp(['********** Texture Metrics **********']) ;
    disp(['                                     Mean of HU inside nodule = ' num2str(mean_gray_2d)]) ;
    disp(['                       Standard Deviation of HU inside nodule = ' num2str(std_gray_2d)]) ;
    disp(['                                 Variance of HU inside nodule = ' num2str(var_gray_2d)]) ;
    disp(['                                 Skewness of HU inside nodule = ' num2str(skewness_gray_2d)]) ;
    disp(['                                 Kurtosis of HU inside nodule = ' num2str(kurtosis_gray_2d)]) ;
    disp(['                 Mean of difference image of HU inside nodule = ' num2str(mean_gray_diff_2d)]) ;
    disp(['   Standard Deviation of difference image of HU inside nodule = ' num2str(std_gray_diff_2d)]) ;
    disp(['             Variance of difference image of HU inside nodule = ' num2str(var_gray_diff_2d)]) ;
    disp(['             Skewness of difference image of HU inside nodule = ' num2str(skewness_gray_diff_2d)]) ;
    disp(['             Kurtosis of difference image of HU inside nodule = ' num2str(kurtosis_gray_diff_2d)]) ;
    disp(['********** Margin Metrics **********']) ;
    disp(['                                              Summed distance = ' num2str(summed_distance_2d)]) ;
    disp(['                                                Mean distance = ' num2str(mean_distance_2d)]) ;
    disp(['                                   Normalized summed distance = ' num2str(norm_summed_distance_2d)]) ;
    disp(['                                     Normalized mean distance = ' num2str(norm_mean_distance_2d)]) ;


    diary off ;
    
    % Create comma-separated values list of output for easy importation
    % into Excel.
    
    diary_name_csv = [filename_segm '_' num2str(bone_index, '%02d') '_csv.diary'] ;
    diary(diary_name_csv) ;
    
%     if(jjj == 1)
        disp_string = [ ...
            filename_segm ',bone_index,height,width,depth,number_voxels_inside_3d,' ...
            'volume_3d,surface_area_3d,sphericity_3d,compactness_3d,r_1_3d,' ...
            'r_2_3d,r_3_3d,mean_gray_3d,std_gray_3d,var_gray_3d,skewness_gray_3d,' ...
            'kurtosis_gray_3d,mean_gray_diff_3d,std_gray_diff_3d,var_gray_diff_3d,skewness_gray_diff_3d,' ...
            'kurtosis_gray_diff_3d,summed_distance_3d,mean_distance_3d,norm_summed_distance_3d,' ...
            'norm_mean_distance_3d,plane_index_of_maximal_area,number_pixels_inside_2d,area_2d,' ...
            'perimeter_2d,circularity_2d,r_1_2d,r_2_2d,mean_gray_2d,' ...
            'std_gray_2d,var_gray_2d,skewness_gray_2d,kurtosis_gray_2d,' ...
            'mean_gray_diff_2d,std_gray_diff_2d,var_gray_diff_2d,skewness_gray_diff_2d,' ...
            'kurtosis_gray_diff_2d,summed_distance_2d,mean_distance_2d,' ...
            'norm_summed_distance_2d,norm_mean_distance_2d'] ;
        disp(disp_string)
%     end

    disp_string = [ ...
        filename_segm ',' num2str(bone_index,'%02d') ',' num2str(height) ',' num2str(width) ',' num2str(depth) ',' num2str(number_voxels_inside_3d) ',' ...
        num2str(volume_3d) ',' num2str(surface_area_3d) ',' num2str(sphericity_3d) ',' num2str(compactness_3d) ',' num2str(r_1_3d) ',' ...
        num2str(r_2_3d) ',' num2str(r_3_3d) ',' num2str(mean_gray_3d) ',' num2str(std_gray_3d) ',' num2str(var_gray_3d) ',' num2str(skewness_gray_3d) ',' ...
        num2str(kurtosis_gray_3d) ',' num2str(mean_gray_diff_3d) ',' num2str(std_gray_diff_3d) ',' num2str(var_gray_diff_3d), ',' num2str(skewness_gray_diff_3d) ',' ...
        num2str(kurtosis_gray_diff_3d) ',' num2str(summed_distance_3d) ',' num2str(mean_distance_3d) ',' num2str(norm_summed_distance_3d) ',' ...
        num2str(norm_mean_distance_3d) ',' num2str(plane_index_of_maximal_area) ',' num2str(number_pixels_inside_2d) ',' num2str(area_2d) ',' ...
        num2str(perimeter_2d) ',' num2str(circularity_2d) ',' num2str(r_1_2d) ',' num2str(r_2_2d) ',' num2str(mean_gray_2d) ',' ...
        num2str(std_gray_2d) ',' num2str(var_gray_2d) ',' num2str(skewness_gray_2d) ',' num2str(kurtosis_gray_2d) ',' ...
        num2str(mean_gray_diff_2d) ',' num2str(std_gray_diff_2d) ',' num2str(var_gray_diff_2d) ',' num2str(skewness_gray_diff_2d) ',' ...
        num2str(kurtosis_gray_diff_2d) ',' num2str(summed_distance_2d) ',' num2str(mean_distance_2d) ',' ...
        num2str(norm_summed_distance_2d) ',' num2str(norm_mean_distance_2d)] ;
    disp(disp_string)
    
    diary off ;
    
end