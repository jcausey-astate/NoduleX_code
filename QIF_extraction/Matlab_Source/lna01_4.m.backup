function[metrics_3d, metrics_2d] = lna(binary_ana, grayscale_ana, segment_num, kernel, gray_seg_thresh, fast_only)

   
   if (exist('fast_only', 'var') == 0)
       printf('fast_only undefined');
       fast_only = 0;
   endif

   tic ;
   disp('Entering lna01_4') ;
   lna_version = 1.4 ;

   % Version 1.2 incorporates the gradient margin metric and the lacunarity
   % calculation. The only known missing  metric is "correlation."
   %
   % Version 1.3 incorporates 5 new texture metrics from:
   %    Moses Amadasun and Robert T. King, "Textural Features Corresponding to
   %    Textural Properties," IEEE Transactions on Systems, Man, and
   %    Cybernetics, Vol. 19, No. 5, September/October 1989.
   %
   % Version 1.4 adds a new statistics, the median HU inside the nodule to the
   %    2D and 3D calculations.

   % Should I also input a path to the data and a path to the results as well
   % as the filenames themselves??????????????

   % bq_full.m (Calculates the 3-D distance to surface metric that takes so
   % long.)
   %
   % Calculate metrics for some lung nodules in support of an ICTS grant
   % proposal due on Monday, November 16, 2009.
   %
   % Metrics have been proposed for 3-D image volumes and for 2-D image
   % slices.

   % Read in the binary-thresholded 3-D image with the nodule only, which is
   % unsigned 8-bit.

   % Also read in the grayscale 3-D image which contains the nodule. It is
   % signed 16-bit.
   %
   % Software has been "re-purposed" to analyze bones rather than lung
   % nodules. DGP 2010-04-29
   %
   % This version includes the variable names and descriptions devised by
   % David Gierada and David Politte for use with the REDCAP database.

   disp(' ') ;
   disp(' ') ;
   disp(' ') ;

   filename_root = strrep(binary_ana, '_bin', '') ;
   diary_name = [filename_root '.diary'] ;
   diary(diary_name) ;

   header_name = [binary_ana '.hdr'] ;
   header = analyze75info(header_name) ;
   dy_mm = double(header.PixelDimensions(1)) ;
   dx_mm = double(header.PixelDimensions(2)) ;
   dz_mm = double(header.PixelDimensions(3)) ;
   image_bin_name = [binary_ana '.img'] ;
   a = int8(analyze75read(image_bin_name)) ;
   temp = size(a) ;
   d3_long_pix = temp(2) ;
   d3_trans_pix = temp(1) ;
   d3_z_pix = temp(3) ;
   s_3d_orig = double(a) ;
   disp_string = ['# rows = '    num2str(d3_trans_pix) ...
               '   # columns = ' num2str(d3_long_pix) ...
               '   # planes = '  num2str(d3_z_pix)] ;
   disp(disp_string) ;
   clear a temp disp_string ;

   image_gray_name = [grayscale_ana '.img'] ;
   g = analyze75read(image_gray_name) ;
   g_3d = double(g) ;
   clear g ;

   if (~isempty(find(mod(s_3d_orig, 1) ~= 0, 1)))
      disp('WARNING: Segmentation image is not binary!') ;
      exit ;
   end
   s_3d = s_3d_orig / max(max(max(s_3d_orig))) ;
   if (~isempty(find(mod(s_3d, 1) ~= 0, 1)))
      disp('WARNING: Segmentation image is not binary!') ;
      exit ;
   end
   clear s_3d_orig ;

   % Ensure that the nodule isn't on a border plane. Some of the
   % algorithms below require 2 "guard" planes on all sides.

   s_3d(1:2,       :,         :        ) = 0 ;
   s_3d(end-1:end, :,         :        ) = 0 ;
   s_3d(:,         1:2,       :        ) = 0 ;
   s_3d(:,         end-1:end, :        ) = 0 ;
   s_3d(:,         :,         1:2      ) = 0 ;
   s_3d(:,         :,         end-1:end) = 0 ;


   %**************************************************************************
   %**************************************************************************
   %**************************************************************************
   % CALCULATE 3-D METRICS
   %**************************************************************************
   %**************************************************************************
   %************************************************lac**************************

   disp('CALCULATING 3-D METRICS ...') ;

   %**************************************************************************
   % Calculate size metrics.
   %**************************************************************************

   disp('   Calculating size metrics ...') ;

   % Calculate the volume.

   disp('      Calculating volume ...') ;
   d3_nodule_voxels = sum(sum(sum(s_3d))) ;
   d3_nodule_volume = d3_nodule_voxels * dx_mm * dy_mm * dz_mm ;
   disp('      Finished calculating volume.') ;

   % Calculate a bounding box.

   disp('      Calculating bounding box ...') ;
   [iii_list, jjj_list, kkk_list] = ind2sub(size(s_3d), find(s_3d > 0)) ;
   d3_bounding_box = [min(iii_list) min(jjj_list) min(kkk_list) ; ...
                      max(iii_list) max(jjj_list) max(kkk_list)] ;
   d3_bounding_x_mm = (d3_bounding_box(2, 2) - d3_bounding_box(1, 2)) * dx_mm ;
   d3_bounding_y_mm = (d3_bounding_box(2, 1) - d3_bounding_box(1, 1)) * dy_mm ;
   d3_bounding_z_mm = (d3_bounding_box(2, 3) - d3_bounding_box(1, 3)) * dz_mm ;
   clear iii_list jjj_list kkk_list ;
   disp('      Finished calculating bounding box.') ;

   % Calculate the surface area.

   disp('      Calculating surface area ...') ;
   type_a_surfaces_3d = sum(sum(sum(abs(s_3d(:,:,1:(end-1)) - s_3d(:,:,2:end))))) ;
   type_b_surfaces_3d = sum(sum(sum(abs(s_3d(:,1:(end-1),:) - s_3d(:,2:end,:))))) ;
   type_c_surfaces_3d = sum(sum(sum(abs(s_3d(1:(end-1),:,:) - s_3d(2:end,:,:))))) ;

   d3_surface_area_mm2 = (type_a_surfaces_3d * dx_mm * dy_mm) ...
                       + (type_b_surfaces_3d * dx_mm * dz_mm) ...
                       + (type_c_surfaces_3d * dy_mm * dz_mm) ;
   disp('      Finished calculating surface area.') ;

   disp('   Finished calculating size metrics.') ;

   %**************************************************************************
   % Calculate shape metrics.
   %**************************************************************************

   disp('   Calculating shape metrics ...') ;

   % Calculate the sphericity.

   disp('      Calculating sphericity ...') ;
   d3_sphericity = (pi^(1/3)) * ((6 * d3_nodule_volume)^(2/3)) / d3_surface_area_mm2 ;
   disp('      Finished calculating sphericity.') ;

   % Calculate the compactness. Compactness is sphericity cubed.

   disp('      Calculating compactness ...') ;
   d3_compactness = (36 * pi * (d3_nodule_volume^2)) / (d3_surface_area_mm2^3) ;
   disp('      Finished calculating compactness.') ;

   % Calculate the principal rotational moments.
   % This calculation uses the binary segmentation image.
   % See http://en.wikipedia.org/wiki/Moment_of_inertia for a definition of the
   % tensor matrix that is calculated.
   % The tensor matrix is then diagonalized using a singular-value
   % decomposition (SVD) and normalized to have unit energy so that it is
   % independent of size.

   disp('      Calculating principal rotational moments ...') ;

   % Find the centroid. This will be a tuple.

   x = 1:d3_long_pix ;
   y = 1:d3_trans_pix ;
   z = 1:d3_z_pix ;
   [X, Y, Z] = meshgrid(x, y, z) ;
   x_c = sum(sum(sum(s_3d .* X))) / d3_nodule_voxels ;
   y_c = sum(sum(sum(s_3d .* Y))) / d3_nodule_voxels ;
   z_c = sum(sum(sum(s_3d .* Z))) / d3_nodule_voxels ;
   clear X Y Z ;

   x_shift = x - x_c ;
   y_shift = y - y_c ;
   z_shift = z - z_c ;
   [X_shift, Y_shift, Z_shift] = meshgrid(x_shift, y_shift, z_shift) ;
   clear x y z x_c y_c z_c x_shift y_shift z_shift ;

   x_sq_term_3d = sum(sum(sum((s_3d .* X_shift).^2))) / d3_nodule_voxels ;
   y_sq_term_3d = sum(sum(sum((s_3d .* Y_shift).^2))) / d3_nodule_voxels ;
   z_sq_term_3d = sum(sum(sum((s_3d .* Z_shift).^2))) / d3_nodule_voxels ;

   xy_term_3d   = sum(sum(sum((s_3d .* X_shift) .* (s_3d .* Y_shift)))) / d3_nodule_voxels ;
   xz_term_3d   = sum(sum(sum((s_3d .* X_shift) .* (s_3d .* Z_shift)))) / d3_nodule_voxels ;
   yz_term_3d   = sum(sum(sum((s_3d .* Y_shift) .* (s_3d .* Z_shift)))) / d3_nodule_voxels ;
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

   clear x_sq_term_3d y_sq_term_3d z_sq_term_3d xy_term_3d xz_term_3d yz_term_3d ;

   I_matrix_3d = [ I11_3d -I12_3d -I13_3d ; ...
                  -I21_3d  I22_3d -I23_3d ; ...
                  -I31_3d -I32_3d  I33_3d] ;
   clear I11_3d I12_3d I13_3d I21_3d I22_3d I23_3d I31_3d I32_3d I33_3d ;

   % Diagonalize it. The r_i_3d are the answers we seek.

   [U, S, V] = svd(I_matrix_3d) ;
   norm_fac = sqrt(sum(diag(S) .^ 2)) ;
   d3_r1 = S(1,1) / norm_fac ;
   d3_r2 = S(2,2) / norm_fac ;
   d3_r3 = S(3,3) / norm_fac ;
   clear norm_fac U S V I_matrix_3d ;

   d3_moment_ratio = d3_r1 / d3_r3 ;

   disp('      Finished calculating principal rotational moments.') ;

   disp('   Finished calculating shape metrics.') ;

   %**************************************************************************
   % Calculate attenuation value metrics.
   %**************************************************************************

   disp('   Calculating attenuation value metrics ...') ;
   index_inside_3d   = find(s_3d == 1) ;
   disp('      Calculating median of attenuation values ...') ;
   d3_median_atten   = median(g_3d(index_inside_3d)) ;
   disp('      Finished calculating median of attenuation values.') ;
   disp('      Calculating mean of attenuation values ...') ;
   d3_mean_atten     = mean(g_3d(index_inside_3d)) ;
   disp('      Finished calculating mean of attenuation values.') ;
   disp('      Calculating standard deviation of attenuation values ...') ;
   d3_stdev_atten    = std(g_3d(index_inside_3d))  ;
   disp('      Finished calculating standard deviation of attenuation values.') ;
   disp('      Calculating variance of attenuation values ...') ;
   d3_var_atten      = var(g_3d(index_inside_3d))  ;
   disp('      Finished calculating variance of attenuation values.') ;
   disp('      Calculating skewness of attenuation values ...') ;
   d3_skew_atten     = skewness(g_3d(index_inside_3d)) ;
   disp('      Finished calculating skewness of attenuation values.') ;
   disp('      Calculating kurtosis of attenuation values ...') ;
   d3_kurt_atten     = kurtosis(g_3d(index_inside_3d)) ;
   disp('      Finished calculating kurtosis of attenuation values.') ;
   disp('      Calculating entropy of attenuation values ...') ;
   d3_entropy_atten  = entropy(g_3d(index_inside_3d)) ;
   disp('      Finished calculating entropy of attenuation values.') ;

   disp('   Finished calculating attenuation value metrics.') ;

   %**************************************************************************
   % Calculate texture metrics.
   % Start by differencing the grayscale images with nearest neighbors.
   %**************************************************************************

   disp('   Calculating texture metrics ...') ;
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
   disp('      Calculating mean of difference images ...') ;
   d3_mean_diff     = mean(g_3d_diff(index_inside_3d)) ;
   disp('      Finished calculating mean of difference images.') ;
   disp('      Calculating standard deviation of difference images ...') ;
   d3_stdev_diff      = std(g_3d_diff(index_inside_3d)) ;
   disp('      Finished calculating standard deviation of difference images.') ;
   disp('      Calculating variance of difference images ...') ;
   d3_var_diff      = var(g_3d_diff(index_inside_3d))  ;
   disp('      Finished calculating variance of difference images.') ;
   disp('      Calculating skewness of difference images ...') ;
   d3_skew_diff = skewness(g_3d_diff(index_inside_3d)) ;
   disp('      Finished calculating skewness of difference images.') ;
   disp('      Calculating kurtosis of difference images ...') ;
   d3_kurt_diff = kurtosis(g_3d_diff(index_inside_3d)) ;
   disp('      Finished calculating kurtosis of difference images.') ;
   disp('      Calculating entropy of difference images ...') ;
   d3_entropy_diff  = entropy(g_3d_diff(index_inside_3d)) ;
   disp('      Finished calculating entropy of difference images.') ;
   disp('      Calculating correlation metric of difference images ...') ;
   disp('         PLACEHOLDER FOR 3D CORRELATION METRIC') ;
   d3_corr_diff = -9999.0 ;
   disp('      Finished calculating correlation metric of difference images.') ;
   disp('      Calculating lacunarity of segmentation images ...') ;
   s_3d_small = s_3d(d3_bounding_box(1,1):d3_bounding_box(2,1), ...
                     d3_bounding_box(1,2):d3_bounding_box(2,2), ...
                     d3_bounding_box(1,3):d3_bounding_box(2,3)) ;
   [d3_temp_lacunarity, d3_box_size_lacunarity] = calc_lacunarity_3d(s_3d_small) ;
   d3_lacunarity_seg = zeros(1, 10) ;
   d3_lacunarity_seg(1:numel(d3_temp_lacunarity)) = d3_temp_lacunarity ;
   d3_lacunarity_seg(numel(d3_temp_lacunarity):end) = d3_temp_lacunarity(end) ;
   disp('      Finished calculating lacunarity of segmentation images.') ;
   disp(['Entering ngtdm: ' datestr(now)]) ;
   [d3_ak_coarseness_1, d3_ak_contrast_1, d3_ak_busyness_1, ...
      d3_ak_complexity_1, d3_ak_texture_strength_1] = ngtdm(g_3d, s_3d, 1) ;
   disp(['Finished with ngtdm: ' datestr(now)]) ;
   disp(['Entering ngtdm: ' datestr(now)]) ;
   [d3_ak_coarseness_2, d3_ak_contrast_2, d3_ak_busyness_2, ...
      d3_ak_complexity_2, d3_ak_texture_strength_2] = ngtdm(g_3d, s_3d, 2) ;
   disp(['Finished with ngtdm: ' datestr(now)]) ;
   clear g_3d_diff q kern ;
   disp('   Finished calculating texture metrics.') ;

   %**************************************************************************
   % Calculate margin metrics.
   %**************************************************************************

   disp('   Calculating margin metrics ...') ;

   % Calculate the mean distance to surface.

   if (fast_only == 0)
      
      disp('      Calculating distance to surface ...') ;
      dist_img_name_3d = [filename_root '_3d.dst'] ;
      fid_3d = fopen(dist_img_name_3d, 'w') ;

      disp('         Finding border voxels ...') ;
      kern = ones(3,3,3) ;
      border_3d = sign((1 - s_3d) .* convn(s_3d, kern, 'same')) ;
      clear kern ;
      [iii_list, jjj_list, kkk_list] = ind2sub(size(s_3d), find(border_3d > 0)) ;
      disp('         Finished finding border voxels.') ;

      disp('         Allocating 3D image for distance to surface ...') ;
      dist_exterior_3d = single(zeros(d3_trans_pix, d3_long_pix)) ;
      disp('         Finished allocating 3D image for distance to surface.') ;

      d3_summed_dist = 0.0 ;

      for k = 1:d3_z_pix
         nnztp = length(find(s_3d(:, :, k) > 0)) ;
         disp_string = ['         Plane ' num2str(k) ' out of ' ...
            num2str(d3_z_pix) ', ' num2str(nnztp) ...
            ' voxels inside the nodule. '   datestr(now())] ;
         disp(disp_string) ;
         clear nnztp disp_string ;
         term_3 = ((k - kkk_list) * dz_mm).^2 ;
         for j = 1:d3_long_pix
            term_2 = ((j - jjj_list) * dx_mm).^2 ;
            for i = 1:d3_trans_pix
               if (s_3d(i,j,k) == 1)
                  term_1 = ((i - iii_list) * dy_mm).^2 ;
                  dist_sq_vec = term_1 + term_2 + term_3 ;
                  dist_exterior_3d(i, j) = sqrt(min(dist_sq_vec)) ;
                  d3_summed_dist = d3_summed_dist + dist_exterior_3d(i, j) ;
               end
            end
         end
         count = fwrite(fid_3d, dist_exterior_3d, 'single') ;
      end
      clear term_1 term_2 dist_sq_vec dist_exterior_3d iii_list jjj_list kkk_list i j k ;
      d3_mean_dist = d3_summed_dist / d3_nodule_voxels ;
      d3_norm_summed_dist = d3_summed_dist / (d3_nodule_volume^(1/3)) ;
      d3_norm_mean_dist = d3_mean_dist / (d3_nodule_volume^(1/3)) ;

      fclose(fid_3d) ;
      clear fid_3d ;
      disp('      Finished calculating distance to surface.') ;
      
   else
      
      disp('      Deliberately skipping calculating distance to surface ...') ;
      d3_summed_dist      = -9999.0 ;
      d3_mean_dist        = -9999.0 ;
      d3_norm_summed_dist = -9999.0 ;
      d3_norm_mean_dist   = -9999.0 ;
      disp('      Finished deliberately skipping calculating distance to surface.') ;

   end

   % Calculate the fractal dimension of the 3D volume.

   disp('      Calculating fractal dimension of the volume ...') ;
   [n_temp_vol_3d, r_temp_vol_3d] = boxcount(s_3d, 'slope') ;
   s_temp_vol_3d = -gradient(log(n_temp_vol_3d)) ./ gradient(log(r_temp_vol_3d)) ;
   d3_fractal_vol = median(s_temp_vol_3d) ;
   pause(5) ;
   close ;
   % clear n_temp_vol_3d r_temp_vol_3d s_temp_vol_3d ;
   disp('      Finished calculating fractal dimension of the volume.') ;

   % Calculate the fractal dimension of the border surface of the 3D volume.

   disp('      Calculating fractal dimension of the border surface of the volume ...') ;
   if (fast_only ~= 0)
      disp('         Finding border voxels ...') ;
      kern = ones(3,3,3) ;
      border_3d = sign((1 - s_3d) .* convn(s_3d, kern, 'same')) ;
      clear kern ;
      disp('         Finished finding border voxels.') ;
   end
   [n_temp_surf_3d, r_temp_surf_3d] = boxcount(border_3d, 'slope') ;
   s_temp_surf_3d = -gradient(log(n_temp_surf_3d)) ./ gradient(log(r_temp_surf_3d)) ;
   d3_fractal_surf = median(s_temp_surf_3d) ;
   pause(5) ;
   close ;
   % clear n_temp_surf_3d r_temp_surf_3d s_temp_surf_3d ;
   disp('      Finished calculating fractal dimension of the border surface of the volume.') ;

   % Calculate volume versus threshold.

   disp('      Calculating volume versus threshold ...') ;
   thresh_3d_vec = -1025:25:1000 ;
   d3_nodule_volume_vec = zeros(size(thresh_3d_vec)) ;
   for index = 1:length(thresh_3d_vec)
      d3_nodule_volume_vec(index) = length(find(g_3d(index_inside_3d) > thresh_3d_vec(index))) * dx_mm * dy_mm * dz_mm ;
   end
   clear index ;
   disp('      Finished calculating volume versus threshold.') ;

   % Calculate gradient margin metric.

   disp('      Calculating gradient margin metric ...') ;
   [iii, jjj, kkk] = ind2sub(size(s_3d), find(border_3d > 0)) ;
   temp_sum = 0 ;
   for mmm = 1:length(iii)
      temp_sum = temp_sum + ...
         (((g_3d(iii(mmm) + 1, jjj(mmm),     kkk(mmm))         ...
          - g_3d(iii(mmm),     jjj(mmm),     kkk(mmm)    ))^2) ...
        + ((g_3d(iii(mmm),     jjj(mmm),     kkk(mmm))         ...
          - g_3d(iii(mmm) - 1, jjj(mmm),     kkk(mmm)    ))^2) ...
        + ((g_3d(iii(mmm),     jjj(mmm) + 1, kkk(mmm))         ...
          - g_3d(iii(mmm),     jjj(mmm),     kkk(mmm)    ))^2) ...
        + ((g_3d(iii(mmm),     jjj(mmm),     kkk(mmm))         ...
          - g_3d(iii(mmm),     jjj(mmm) - 1, kkk(mmm)    ))^2) ...
        + ((g_3d(iii(mmm),     jjj(mmm),     kkk(mmm) + 1)     ...
          - g_3d(iii(mmm),     jjj(mmm),     kkk(mmm)    ))^2) ...
        + ((g_3d(iii(mmm),     jjj(mmm),     kkk(mmm))         ...
          - g_3d(iii(mmm),     jjj(mmm),     kkk(mmm) - 1))^2)) ;
   end
   d3_grad_margin = sqrt(temp_sum / (length(iii) * 6)) ;
   disp('      Finished calculating gradient margin metric.') ;

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
   d2_plane = 0 ;
   for i = 1:d3_z_pix
      test = sum(sum(s_3d(:,:,i))) ;
      if (test > max_inside_2d)
         max_inside_2d = test ;
         d2_plane = i ;
      end
   end
   s_2d = s_3d(:,:,d2_plane) ;
   g_2d = g_3d(:,:,d2_plane) ;
   clear max_inside_2d i test ;

   %**************************************************************************
   % Calculate size metrics ...') ;
   %**************************************************************************

   disp('   Calculating size metrics ...') ;

   % Calculate the area.

   disp('      Calculating area ...') ;
   d2_nodule_pixels = sum(sum(s_2d)) ;
   d2_nodule_area = d2_nodule_pixels * dx_mm * dy_mm ;
   disp('      Finished calculating area.') ;

   % Calculate a bounding box.

   disp('      Calculating bounding box ...') ;
   [iii_list, jjj_list] = ind2sub(size(s_2d), find(s_2d > 0)) ;
   d2_bounding_box = [min(iii_list) min(jjj_list) ; ...
                      max(iii_list) max(jjj_list)] ;
   d2_bounding_x_mm = (d2_bounding_box(2, 2) - d2_bounding_box(1, 2)) * dx_mm ;
   d2_bounding_y_mm = (d2_bounding_box(2, 1) - d2_bounding_box(1, 1)) * dy_mm ;
   clear iii_list jjj_list ;
   disp('      Finished calculating bounding box.') ;

   % Calculate the perimeter.

   disp('      Calculating perimeter ...') ;
   type_b_edges_2d = sum(sum(abs(s_2d(:,1:(end-1)) - s_2d(:,2:end)))) ;
   type_c_edges_2d = sum(sum(abs(s_2d(1:(end-1),:) - s_2d(2:end,:)))) ;

   d2_perim_mm = (type_b_edges_2d * dx_mm) ...
               + (type_c_edges_2d * dy_mm) ;
   disp('      Finished calculating perimeter.') ;

   disp('   Finished calculating size metrics.') ;

   %**************************************************************************
   % Calculate shape metrics.
   %**************************************************************************

   disp('   Calculating shape metrics ...') ;

   % Calculate the circularity (McNitt-Gray's definition).

   disp('      Calculating circularity ...') ;
   d2_circularity = (4 * pi * d2_nodule_area) / (d2_perim_mm^2) ;
   disp('      Finished calculating circularity.') ;

   % Calculate the principal rotational moments.
   % This calculation uses the binary segmentation image.
   % See http://en.wikipedia.org/Moment_of_inertia for a definition of the
   % tensor matrix that is calculated (modified for 2D).
   % The tensor matrix is then diagonalized using a singular-value
   % decomposition (SVD) and normalized to have unit energy so that it is
   % independent of size.

   disp('      Calculating principal rotational moments ...') ;

   % Find the centroid. This will be an ordered pair.

   x = 1:d3_long_pix ;
   y = 1:d3_trans_pix ;
   [X, Y] = meshgrid(x, y) ;
   x_c = sum(sum(s_2d .* X)) / d2_nodule_pixels ;
   y_c = sum(sum(s_2d .* Y)) / d2_nodule_pixels ;
   clear X Y ;

   x_shift = x - x_c ;
   y_shift = y - y_c ;
   [X_shift, Y_shift] = meshgrid(x_shift, y_shift) ;
   clear x y x_c y_c x_shift y_shift ;

   x_sq_term_2d = sum(sum((s_2d .* X_shift).^2)) / d2_nodule_pixels ;
   y_sq_term_2d = sum(sum((s_2d .* Y_shift).^2)) / d2_nodule_pixels ;

   xy_term_2d   = sum(sum((s_2d .* X_shift) .* (s_2d .* Y_shift))) / d2_nodule_pixels ;

   clear X_shift Y_shift Z_shift ;

   % Define the tensor as shown on the wikipedia web page, but modified for
   % 2D.

   I11_2d = y_sq_term_2d ;
   I22_2d = x_sq_term_2d ;

   I12_2d = xy_term_2d ;
   I21_2d = I12_2d ;

   clear x_sq_term_2d y_sq_term_2d xy_term_2d ;

   I_matrix_2d = [ I11_2d -I12_2d ; ...
                  -I21_2d  I22_2d] ;
   clear I11_2d I12_2d I21_2d I22_2d ;

   % Diagonalize it. The r_i_2d are the answers we seek.

   [U, S, V] = svd(I_matrix_2d) ;
   norm_fac = sqrt(sum(diag(S) .^ 2)) ;
   d2_r1 = S(1,1) / norm_fac ;
   d2_r2 = S(2,2) / norm_fac ;
   clear norm_fac U S V I_matrix_2d ;

   d2_moment_ratio = d2_r1 / d2_r2 ;

   disp('      Finished calculating principal rotational moments.') ;

   disp('   Finished calculating shape metrics.') ;

   %**************************************************************************
   % Calculate attenuation value metrics.
   %**************************************************************************

   disp('   Calculating attenuation value metrics ...') ;
   index_inside_2d  = find(s_2d == 1) ;
   disp('      Calculating median of attenuation values ...') ;
   d2_median_atten   = median(g_2d(index_inside_2d)) ;
   disp('      Finished calculating median of attenuation values.') ;
   disp('      Calculating mean of attenuation values ...') ;
   d2_mean_atten     = mean(g_2d(index_inside_2d)) ;
   disp('      Finished calculating mean of attenuation values.') ;
   disp('      Calculating standard deviation of attenuation values ...') ;
   d2_stdev_atten      = std(g_2d(index_inside_2d))  ;
   disp('      Finished calculating standard deviation of attenuation values.') ;
   disp('      Calculating variance of attenuation values ...') ;
   d2_var_atten      = var(g_2d(index_inside_2d))  ;
   disp('      Finished calculating variance of attenuation values.') ;
   disp('      Calculating skewness of attenuation values ...') ;
   d2_skew_atten = skewness(g_2d(index_inside_2d)) ;
   disp('      Finished calculating skewness of attenuation values.') ;
   disp('      Calculating kurtosis of attenuation values ...') ;
   d2_kurt_atten = kurtosis(g_2d(index_inside_2d)) ;
   disp('      Finished calculating kurtosis of attenuation values.') ;
   disp('      Calculating entropy of attenuation values ...') ;
   d2_entropy_atten  = entropy(g_2d(index_inside_2d)) ;
   disp('      Finished calculating entropy of attenuation values.') ;

   disp('   Finished calculating attenuation value metrics.') ;

   %**************************************************************************
   % Calculate texture metrics.
   % Start by differencing the grayscale images with nearest neighbors.
   %**************************************************************************

   disp('   Calculating texture metrics ...') ;
   kern = zeros(3, 3) ;
   q = -1/4 ;
   kern = [0 q 0 ; ...
           q 1 q ; ...
           0 q 0] ;
   g_2d_diff = convn(g_2d, kern, 'same') ;
   disp('      Calculating mean of difference images ...') ;
   d2_mean_diff     = mean(g_2d_diff(index_inside_2d)) ;
   disp('      Finished calculating mean of difference images.') ;
   disp('      Calculating standard deviation of difference images ...') ;
   d2_stdev_diff      = std(g_2d_diff(index_inside_2d)) ;
   disp('      Finished calculating standard deviation of difference images.') ;
   disp('      Calculating variance of difference images ...') ;
   d2_var_diff      = var(g_2d_diff(index_inside_2d)) ;
   disp('      Finished calculating variance of difference images.') ;
   disp('      Calculating skewness of difference images ...') ;
   d2_skew_diff = skewness(g_2d_diff(index_inside_2d)) ;
   disp('      Finished calculating skewness of difference images.') ;
   disp('      Calculating kurtosis of difference images ...') ;
   d2_kurt_diff = kurtosis(g_2d_diff(index_inside_2d)) ;
   disp('      Finished calculating kurtosis of difference images.') ;
   disp('      Calculating entropy of difference images ...') ;
   d2_entropy_diff  = entropy(g_2d_diff(index_inside_2d)) ;
   disp('      Finished calculating entropy of difference images.') ;
   disp('      Calculating correlation metric of difference images ...') ;
   disp('         PLACEHOLDER FOR 2D CORRELATION METRIC') ;
   d2_corr_diff = -9999.0 ;
   disp('      Finished calculating correlation metric of difference images.') ;
   disp('      Calculating lacunarity of segmentation images ...') ;
   s_2d_small = s_2d(d2_bounding_box(1,1):d2_bounding_box(2,1), ...
                     d2_bounding_box(1,2):d2_bounding_box(2,2)) ;
   [d2_temp_lacunarity, d2_box_size_lacunarity] = calc_lacunarity_2d(s_2d_small) ;
   d2_lacunarity_seg = zeros(1, 10) ;
   d2_lacunarity_seg(1:numel(d2_temp_lacunarity)) = d2_temp_lacunarity ;
   d2_lacunarity_seg(numel(d2_temp_lacunarity):end) = d2_temp_lacunarity(end) ;
   disp('      Finished calculating lacunarity of segmentation images.') ;
   disp(['Entering ngtdm: ' datestr(now)]) ;
   [d2_ak_coarseness_1, d2_ak_contrast_1, d2_ak_busyness_1, ...
      d2_ak_complexity_1, d2_ak_texture_strength_1] = ngtdm(g_2d, s_2d, 1) ;
   disp(['Finished with ngtdm: ' datestr(now)]) ;
   disp(['Entering ngtdm: ' datestr(now)]) ;
   [d2_ak_coarseness_2, d2_ak_contrast_2, d2_ak_busyness_2, ...
      d2_ak_complexity_2, d2_ak_texture_strength_2] = ngtdm(g_2d, s_2d, 2) ;
   disp(['Finished with ngtdm: ' datestr(now)]) ;
   clear g_2d_diff q kern ;
   disp('   Finished calculating texture metrics.') ;

   %**************************************************************************
   % Calculate margin metrics.
   %**************************************************************************

   disp('   Calculating margin metrics ...') ;

   % Calculate the mean distance to surface.

   disp('      Calculating distance to surface ...') ;
   dist_img_name_2d = [filename_root '_2d.dst'] ;
   fid_2d = fopen(dist_img_name_2d, 'w') ;

   disp('         Finding border voxels ...') ;
   kern = ones(3,3) ;
   border_2d = sign((1 - s_2d) .* convn(s_2d, kern, 'same')) ;
   clear kern ;
   [iii_list, jjj_list] = ind2sub(size(s_2d), find(border_2d > 0)) ;
   disp('         Finished finding border voxels.') ;

   disp('         Allocating 2D image for distance to surface ...') ;
   dist_exterior_2d = single(zeros(d3_trans_pix, d3_long_pix)) ;
   disp('         Finished allocating 2D image for distance to surface.') ;

   d2_summed_dist = 0.0 ;

   for j = 1:d3_long_pix
      term_2 = ((j - jjj_list) * dx_mm).^2 ;
      for i = 1:d3_trans_pix
         if (s_2d(i,j) == 1)
            term_1 = ((i - iii_list) * dy_mm).^2 ;
            dist_sq_vec = term_1 + term_2 ;
            dist_exterior_2d(i, j) = sqrt(min(dist_sq_vec)) ;
            d2_summed_dist = d2_summed_dist + dist_exterior_2d(i, j) ;
         end
      end
   end
   count = fwrite(fid_2d, dist_exterior_2d, 'single') ;
   clear count dist_sq_vec term_1 term_2 dist_exterior_2d iii_list jjj_list i j ;
   d2_mean_dist = d2_summed_dist / d2_nodule_pixels ;
   d2_norm_summed_dist = d2_summed_dist / (d2_nodule_area^(1/2)) ;
   d2_norm_mean_dist = d2_mean_dist / (d2_nodule_area^(1/2)) ;

   fclose(fid_2d) ;
   clear fid_2d ;
   disp('      Finished calculating distance to surface.') ;

   % Calculate the fractal dimension of the 2D area.

   disp('      Calculating fractal dimension of the area ...') ;
   [n_temp_area_2d, r_temp_area_2d] = boxcount(s_2d, 'slope') ;
   s_temp_area_2d = -gradient(log(n_temp_area_2d)) ./ gradient(log(r_temp_area_2d)) ;
   d2_fractal_area = median(s_temp_area_2d) ;
   pause(5) ;
   close ;
   % clear n_temp_area_2d r_temp_area_2d s_temp_area_2d ;
   disp('      Finished calculating fractal dimension of the area.') ;

   % Calculate the fractal dimension of the perimeter of the 2D area.

   disp('      Calculating fractal dimension of the perimeter of the area ...') ;
   [n_temp_perim_2d, r_temp_perim_2d] = boxcount(border_2d, 'slope') ;
   s_temp_perim_2d = -gradient(log(n_temp_perim_2d)) ./ gradient(log(r_temp_perim_2d)) ;
   d2_fractal_perim = median(s_temp_perim_2d) ;
   pause(5) ;
   close ;
   % clear n_temp_perim_2d r_temp_perim_2d s_temp_perim_2d ;
   disp('      Finished calculating fractal dimension of the perimeter of the area ...') ;

   % Calculate area versus threshold.

   disp('      Calculating area versus threshold ...') ;
   thresh_2d_vec = -1025:25:1000 ;
   d2_nodule_volume_vec = zeros(size(thresh_2d_vec)) ;
   for index = 1:length(thresh_2d_vec)
      d2_nodule_volume_vec(index) = length(find(g_2d(index_inside_2d) > thresh_2d_vec(index))) * dx_mm * dy_mm ;
   end
   clear index ;
   disp('      Finished calculating area versus threshold.') ;

   % Calculate gradient margin metric.

   disp('      Calculating gradient margin metric ...') ;
   [iii, jjj] = ind2sub(size(s_2d), find(border_2d > 0)) ;
   temp_sum = 0 ;
   for mmm = 1:length(iii)
      temp_sum = temp_sum + ...
         (((g_2d(iii(mmm) + 1, jjj(mmm)    )     - ...
            g_2d(iii(mmm),     jjj(mmm)    ))^2)   ...
        + ((g_2d(iii(mmm),     jjj(mmm)    )     - ...
            g_2d(iii(mmm) - 1, jjj(mmm)    ))^2)   ...
        + ((g_2d(iii(mmm),     jjj(mmm) + 1)       ...
          - g_2d(iii(mmm),     jjj(mmm)    ))^2)   ...
        + ((g_2d(iii(mmm),     jjj(mmm)    )       ...
          - g_2d(iii(mmm),     jjj(mmm) - 1))^2)) ;
   end
   d2_grad_margin = sqrt(temp_sum / (length(iii) * 4)) ;
   disp('      Finished calculating gradient margin metric.') ;

   disp('   Finished calculating margin metrics.') ;

   disp('FINISHED CALCULATING 2-D METRICS.') ;

   % Save the results in a Matlab ".mat" file. Clear the big arrays first.

   clear ans g_2d g_3d s_2d s_3d index_inside_2d index_inside_3d ;
   filename_mat  = [filename_root '.mat'] ;
   eval_string = ['save ''' filename_mat ''''] ;
   eval(eval_string) ;

   clear disp_string eval_string ;

   % Package the results in a nice dictionary for return.

   format long ;

   algo_stats = {};
   algo_stats{end+1} = {'Filename of grayscale image', grayscale_ana};
   algo_stats{end+1} = {'Filename of binary image', binary_ana};
   algo_stats{end+1} = {'Kernel' , kernel};
   algo_stats{end+1} = {'Delta d3_trans_pix (mm)', dy_mm};
   algo_stats{end+1} = {'Delta d3_long_pix (mm)',  num2str(dx_mm)};
   algo_stats{end+1} = {'Delta d3_z_pix (mm)', num2str(dz_mm)};
   algo_stats{end+1} = {'Height of input image (# rows)', num2str(d3_trans_pix)};
   algo_stats{end+1} = {'Width of input image (# columns)', num2str(d3_long_pix)};
   algo_stats{end+1} = {'Depth of input image (# planes)', num2str(d3_z_pix)};
   algo_stats{end+1} = {'Segmentation Number', segment_num};
   algo_stats{end+1} = {'Segmentation Threshold', gray_seg_thresh};
   algo_stats{end+1} = {'Lung nodule analysis (lna) version' , lna_version};


   disp(['****************************** STATISTICS FOR 3D ******************************']) ;
   metrics_3d = {};
   % disp(['********** Size Metrics **********']) ;
   metrics_3d{end+1} = {'Number of voxels inside', d3_nodule_voxels};
   metrics_3d{end+1} = {'Volume', d3_nodule_volume};
   metrics_3d{end+1} = {'Span of nodule in x direction (mm)', d3_bounding_x_mm};
   metrics_3d{end+1} = {'Span of nodule in y direction (mm)', d3_bounding_y_mm};
   metrics_3d{end+1} = {'Span of nodule in z direction (mm)', d3_bounding_z_mm};
   metrics_3d{end+1} = {'Surface Area', d3_surface_area_mm2};

   % disp(['********** Shape Metrics **********']) ; 
   metrics_3d{end+1} = {'Sphericity', d3_sphericity};
   metrics_3d{end+1} = {'Compactness', d3_compactness};
   metrics_3d{end+1} = {'Rotational moment, r_1', d3_r1};
   metrics_3d{end+1} = {'Rotational moment, r_2', d3_r2};
   metrics_3d{end+1} = {'Rotational moment, r_3', d3_r3};
   metrics_3d{end+1} = {'Ratio of largest to smallest rotational moment', d3_moment_ratio};

   % disp(['********** Attenuation Metrics **********']) ;
   metrics_3d{end+1} = {'Median of HU inside nodule' , d3_median_atten};
   metrics_3d{end+1} = {'Mean of HU inside nodule' , d3_mean_atten};
   metrics_3d{end+1} = {'Standard Deviation of HU inside nodule' , d3_stdev_atten};
   metrics_3d{end+1} = {'Variance of HU inside nodule' , d3_var_atten};
   metrics_3d{end+1} = {'Skewness of HU inside nodule' , d3_skew_atten};
   metrics_3d{end+1} = {'Kurtosis of HU inside nodule' , d3_kurt_atten};
   metrics_3d{end+1} = {'Entropy of HU inside nodule' , d3_entropy_atten};

   % disp(['********** Texture Metrics **********']) ;
   metrics_3d{end+1} = {'Mean of difference image of HU inside nodule' , d3_mean_diff};
   metrics_3d{end+1} = {'Standard Deviation of difference image of HU inside nodule' , d3_stdev_diff};
   metrics_3d{end+1} = {'Variance of difference image of HU inside nodule' , d3_var_diff};
   metrics_3d{end+1} = {'Skewness of difference image of HU inside nodule' , d3_skew_diff};
   metrics_3d{end+1} = {'Kurtosis of difference image of HU inside nodule' , d3_kurt_diff};
   metrics_3d{end+1} = {'Entropy of difference image of HU inside nodule' , d3_entropy_diff};
   metrics_3d{end+1} = {'Correlation metric of difference image of HU inside nodule' , d3_corr_diff};
   metrics_3d{end+1} = {'Lacunarity of segmentation image' , d3_lacunarity_seg};
   metrics_3d{end+1} = {'Box size for lacunarity calculation' , d3_box_size_lacunarity};
   metrics_3d{end+1} = {'Coarseness (d = 1)' , d3_ak_coarseness_1};
   metrics_3d{end+1} = {'Contrast (d = 1)' , d3_ak_contrast_1};
   metrics_3d{end+1} = {'Busyness (d = 1)' , d3_ak_busyness_1};
   metrics_3d{end+1} = {'Complexity (d = 1)' , d3_ak_complexity_1};
   metrics_3d{end+1} = {'Texture Strength (d = 1)' , d3_ak_texture_strength_1};
   metrics_3d{end+1} = {'Coarseness (d = 2)' , d3_ak_coarseness_2};
   metrics_3d{end+1} = {'Contrast (d = 2)' , d3_ak_contrast_2};
   metrics_3d{end+1} = {'Busyness (d = 2)' , d3_ak_busyness_2};
   metrics_3d{end+1} = {'Complexity (d = 2)' , d3_ak_complexity_2};
   metrics_3d{end+1} = {'Texture Strength (d = 2)' , d3_ak_texture_strength_2};

   % disp(['********** Margin Metrics **********']) ;
   metrics_3d{end+1} = {'Summed distance' , d3_summed_dist};
   metrics_3d{end+1} = {'Mean distance' , d3_mean_dist};
   metrics_3d{end+1} = {'Normalized summed distance' , d3_norm_summed_dist};
   metrics_3d{end+1} = {'Normalized mean distance' , d3_norm_mean_dist};
   metrics_3d{end+1} = {'Fractal dimension of volume' , d3_fractal_vol};
   metrics_3d{end+1} = {'Fractal dimension of surface' , d3_fractal_surf};
   metrics_3d{end+1} = {'cum_volume_thresh', thresh_3d_vec};
   metrics_3d{end+1} = {'cum_volume_value', d3_nodule_volume_vec};
   metrics_3d{end+1} = {'Gradient margin metric', d3_grad_margin};


   % disp(['****************************** STATISTICS FOR 2D ******************************']) ;
   metrics_2d = {};
   metrics_2d{end+1} = {'Transverse plane with maximal area' , d2_plane};
   % disp(['********** Size Metrics **********']) ;
   metrics_2d{end+1} = {'Number of pixels inside' , d2_nodule_pixels};
   metrics_2d{end+1} = {'Area' , d2_nodule_area};
   metrics_2d{end+1} = {'Span of nodule in x direction (mm)' , d2_bounding_x_mm};
   metrics_2d{end+1} = {'Span of nodule in y direction (mm)' , d2_bounding_y_mm};
   metrics_2d{end+1} = {'Perimeter' , d2_perim_mm};
   % disp(['********** Shape Metrics **********']) ;
   metrics_2d{end+1} = {'Circularity' , d2_circularity};
   metrics_2d{end+1} = {'Rotational moment, r_1' , d2_r1};
   metrics_2d{end+1} = {'Rotational moment, r_2' , d2_r2};
   metrics_2d{end+1} = {'Ratio of largest to smallest rotational moment' , d2_moment_ratio};
   % disp(['********** Attenuation Metrics **********']) ;
   metrics_2d{end+1} = {'Median of HU inside nodule' , d2_median_atten};
   metrics_2d{end+1} = {'Mean of HU inside nodule' , d2_mean_atten};
   metrics_2d{end+1} = {'Standard Deviation of HU inside nodule' , d2_stdev_atten};
   metrics_2d{end+1} = {'Variance of HU inside nodule' , d2_var_atten};
   metrics_2d{end+1} = {'Skewness of HU inside nodule' , d2_skew_atten};
   metrics_2d{end+1} = {'Kurtosis of HU inside nodule' , d2_kurt_atten};
   metrics_2d{end+1} = {'Entropy of HU inside nodule' , d2_entropy_atten};
   % disp(['********** Texture Metrics **********']) ;
   metrics_2d{end+1} = {'Mean of difference image of HU inside nodule' , d2_mean_diff};
   metrics_2d{end+1} = {'Standard Deviation of difference image of HU inside nodule' , d2_stdev_diff};
   metrics_2d{end+1} = {'Variance of difference image of HU inside nodule' , d2_var_diff};
   metrics_2d{end+1} = {'Skewness of difference image of HU inside nodule' , d2_skew_diff};
   metrics_2d{end+1} = {'Kurtosis of difference image of HU inside nodule' , d2_kurt_diff};
   metrics_2d{end+1} = {'Entropy of difference image of HU inside nodule' , d2_entropy_diff};
   metrics_2d{end+1} = {'Correlation metric of difference image of HU inside nodule' , d2_corr_diff};
   metrics_2d{end+1} = {'Lacunarity of segmentation image' , d2_lacunarity_seg};
   metrics_2d{end+1} = {'Box size for lacunarity calculation' , d2_box_size_lacunarity};
   metrics_2d{end+1} = {'Coarseness (d = 1)' , d2_ak_coarseness_1};
   metrics_2d{end+1} = {'Contrast (d = 1)' , d2_ak_contrast_1};
   metrics_2d{end+1} = {'Busyness (d = 1)' , d2_ak_busyness_1};
   metrics_2d{end+1} = {'Complexity (d = 1)' , d2_ak_complexity_1};
   metrics_2d{end+1} = {'Texture Strength (d = 1)' , d2_ak_texture_strength_1};
   metrics_2d{end+1} = {'Coarseness (d = 2)' , d2_ak_coarseness_2};
   metrics_2d{end+1} = {'Contrast (d = 2)' , d2_ak_contrast_2};
   metrics_2d{end+1} = {'Busyness (d = 2)' , d2_ak_busyness_2};
   metrics_2d{end+1} = {'Complexity (d = 2)' , d2_ak_complexity_2};
   metrics_2d{end+1} = {'Texture Strength (d = 2)' , d2_ak_texture_strength_2};
   % disp(['********** Margin Metrics **********']) ;
   metrics_2d{end+1} = {'Summed distance' , d2_summed_dist};
   metrics_2d{end+1} = {'Mean distance' , d2_mean_dist};
   metrics_2d{end+1} = {'Normalized summed distance' , d2_norm_summed_dist};
   metrics_2d{end+1} = {'Normalized mean distance' , d2_norm_mean_dist};
   metrics_2d{end+1} = {'Fractal dimension of area' , d2_fractal_area};
   metrics_2d{end+1} = {'Fractal dimension of perimeter' , d2_fractal_perim};
   metrics_2d{end+1} = {'cum_volume_thresh', thresh_2d_vec };
   metrics_2d{end+1} = {'cum volume_value', d2_nodule_volume_vec };
   metrics_2d{end+1} = {'Gradient margin metric', d2_grad_margin };

   diary off ;

   disp('Exiting lna01_4') ;
   toc ;

endfunction
