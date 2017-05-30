% Calculate the Neighboring Gray-Tone Difference Matrix (NGTDM)
%
% See:
%    Moses Amadasun and Robert T. King, "Textural Features Corresponding to
%    Texural Properties," IEEE Transactions on Systems, Man, and
%    Cybernetics. Vol. 19, No. 5, September/October 1989, pp. 1264-1274.
%
% This version is for 2D and 3D.
%
% NGTDM stands for "Neighborhood Gray-Tone Difference Matrix.
%
% g_nd is the input gray-scale image in two  or three dimensions.
% s_nd is the input binary-valued segmentation image in two or three
%      dimensions.
%
% d is the neighborhood size. I.e., if d=1 then the neighborhood is a
% region surrounding a pixel of size (2d+1) x (2d+1) = 3 x 3 or
% (2d+1) x (2d+1) x (2d+1) = 3 x 3 x 3.
% If d=2 then the neighborhod is of size 5 x 5 or 5 x 5 x 5.

function [ak_coarseness, ak_contrast, ak_busyness, ak_complexity, ...
   ak_texture_strength] = ngtdm(g_nd, s_nd, d)


   num_dim = length(size(s_nd)) ;

   if (num_dim == 2)
      
      disp('      Calculating NGTDM from Amadasun and King ...') ;

      g_min = min(min(g_nd)) ;
      g_max = max(max(g_nd)) ;
      g_2d = (g_nd - g_min) + 1 ;

      s = zeros(g_max - g_min + 1, 1) ;
      p = zeros(size(s)) ;
      num_this_gray = zeros(size(s)) ;

      sqrt_W = (2 * d) + 1 ;
      W = sqrt_W^2 ;
      Wm1 = W - 1 ;

      % Create a (row, col) list of pixels that are inside the binary object
      % and that have a "guard band" sufficiently big so that they can be
      % included in the calulations of the NGTDM.

      [iii_list, jjj_list] = ind2sub(size(s_nd), find(s_nd > 0)) ;
      iii_list_pruned = zeros(size(iii_list)) ;
      jjj_list_pruned = zeros(size(jjj_list)) ;

      count = 0 ;
      for m = 1:length(iii_list)
         if (sum(sum(s_nd(iii_list(m)-d:iii_list(m)+d, ...
                          jjj_list(m)-d:jjj_list(m)+d))) == W)
            count = count + 1 ;
            iii_list_pruned(count) = iii_list(m) ;
            jjj_list_pruned(count) = jjj_list(m) ;
         end
      end
      iii_list_pruned = iii_list_pruned(1:count) ;
      jjj_list_pruned = jjj_list_pruned(1:count) ;

      % Now loop through all of the included voxels and calculate tne NGTDM.

      for m = 1:length(iii_list_pruned)
         i = iii_list_pruned(m) ;
         j = jjj_list_pruned(m) ;
         index = g_2d(i, j) ;
         A = (sum(sum(g_2d(i-d:i+d,j-d:j+d))) - g_2d(i,j)) / Wm1 ;
         s(index) = s(index) + abs(index - A) ;
         num_this_gray(index) = num_this_gray(index) + 1 ;
      end

      p = num_this_gray / sum(num_this_gray) ;
      
      disp('      Finished calculating NGTDM from Amadasun and King.') ;

   elseif (num_dim == 3)
      
      disp('      Calculating NGTDM from Amadasun and King ...') ;

      g_min = min(min(min(g_nd))) ;
      g_max = max(max(max(g_nd))) ;
      g_3d = (g_nd - g_min) + 1 ;

      s = zeros(g_max - g_min + 1, 1) ;
      p = zeros(size(s)) ;
      num_this_gray = zeros(size(s)) ;

      cube_root_W = (2 * d) + 1 ;
      W = cube_root_W^3 ;
      Wm1 = W - 1 ;

      % Create a (row, col, plane) list of pixels that are inside the binary
      % object and that have a "guard band" sufficiently big so that they can
      % be included in the calculations of the NGTDM.

      [iii_list, jjj_list, kkk_list] = ind2sub(size(s_nd), find(s_nd > 0)) ;
      iii_list_pruned = zeros(size(iii_list)) ;
      jjj_list_pruned = zeros(size(jjj_list)) ;
      kkk_list_pruned = zeros(size(kkk_list)) ;

      count = 0 ;
      for m = 1:length(iii_list)
         if (sum(sum(sum(s_nd(iii_list(m)-d:iii_list(m)+d, ...
                              jjj_list(m)-d:jjj_list(m)+d, ...
                              kkk_list(m)-d:kkk_list(m)+d)))) == W)
            count = count + 1 ;
            iii_list_pruned(count) = iii_list(m) ;
            jjj_list_pruned(count) = jjj_list(m) ;
            kkk_list_pruned(count) = kkk_list(m) ;
         end
      end
      iii_list_pruned = iii_list_pruned(1:count) ;
      jjj_list_pruned = jjj_list_pruned(1:count) ;
      kkk_list_pruned = kkk_list_pruned(1:count) ;

      % Now loop through all of the included voxels and calculate the NGTDM.

      for m = 1:length(iii_list_pruned)
         i = iii_list_pruned(m) ;
         j = jjj_list_pruned(m) ;
         k = kkk_list_pruned(m) ;
         index = g_3d(i, j, k) ;
         A = (sum(sum(sum(g_3d(i-d:i+d,j-d:j+d,k-d:k+d)))) - g_3d(i,j,k)) / Wm1 ;
         s(index) = s(index) + abs(index - A) ;
         num_this_gray(index) = num_this_gray(index) + 1 ;
      end

      p = num_this_gray / sum(num_this_gray) ;
      
      disp('      Finished calculating NGTDM from Amadasun and King.') ;

   else

      disp('ERROR: The binary segmentation image is neither 2D nor 3D') ;
      exit ;

   end
   
   % Next, we'll compute the 5 texture measures, but first we'll compute
   % some variables that are common to more than 1 of the measures.
   
   epsilon = 10^(-7) ; % This was Amadasun and King's choice
   i_vec_pruned = find(p > 0) ;
   p_pruned = p(i_vec_pruned) ;
   s_pruned = s(i_vec_pruned) ;
   num_pix_or_vox = length(iii_list_pruned) ;
   N_g = length(i_vec_pruned) ;
   i_mat = meshgrid(i_vec_pruned) ;
   j_mat = i_mat' ;
   pi_mat = meshgrid(p_pruned) ;
   pj_mat = pi_mat' ;

   % Now compute the first texture metric, coarseness.

   disp('      Computing coarseness (Amadasun and King, 1989) ...') ;
   ak_coarseness = 1 / (epsilon + sum(p_pruned .* s_pruned)) ;
   disp('      Finished computing coarseness.') ;

   % Now compute the second texture metric, contrast.

   disp('      Computing contrast (Amadasun and King, 1989) ...') ;
   factor_1 = 1 / (N_g * (N_g - 1)) ;
%    factor_2 = sum(sum(p_pruned * p_pruned' * toeplitz((i_vec_pruned - 1).^2))) ;
   factor_2 = sum(sum(pi_mat .* pj_mat .* ((i_mat - j_mat).^2))) ;
   factor_3 = sum(s_pruned) / num_pix_or_vox ;
   ak_contrast = factor_1 * factor_2 * factor_3 ;
   disp('      Finished computing contrast.') ;

   % Now compute the third texture metric, busyness.

   disp('      Computing busyness (Amadasun and King, 1989) ...') ;
   numer = sum(p_pruned .* s_pruned) ;
   ip_mat = meshgrid(i_vec_pruned .* p_pruned) ;
   jp_mat = ip_mat' ;
   denom = sum(sum(abs(ip_mat - jp_mat))) ;
   ak_busyness = numer / denom ;
   disp('      Finished computing busyness.') ;

   % Now compute the fourth texture metric, complexity.

   disp('      Computing complexity (Amadasun and King, 1989) ...') ;
   psi_mat = meshgrid(p_pruned .* s_pruned) ;
   psj_mat = psi_mat' ;
   ak_complexity = sum(sum((abs(i_mat - j_mat) ./ (pi_mat + pj_mat)) ...
      .* (psi_mat + psj_mat))) / num_pix_or_vox ;
   disp('      Finished computing complexity.') ;

   % Now compute the fifth texture metric, texture_strength.

   disp('      Computing texture strength (Amadasun and King, 1989) ...') ;
   numer = sum(sum((pi_mat + pj_mat) .* ((i_mat - j_mat).^2))) ;
   denom = epsilon + sum(s_pruned) ;
   ak_texture_strength = numer / denom ;
   disp('      Finished computing texture strength.') ;

end