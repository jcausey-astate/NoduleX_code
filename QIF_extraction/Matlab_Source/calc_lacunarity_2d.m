function [lacunarity_vec, box_size_vec] = calc_lacunarity_2d(image_2d)

temp_a = size(image_2d) ;
num_rows = temp_a(1) ;
num_cols = temp_a(2) ;

temp_b = floor(log(min(temp_a)) / log(2)) ;
box_size_vec = 2.^(0:temp_b) ;
lacunarity_vec = zeros(size(box_size_vec)) ;

% Special case for box size equal to 1.

n = num_rows * num_cols ;
cum_sum    = sum(sum(image_2d   )) ;
cum_sum_sq = sum(sum(image_2d.^2)) ;
lacunarity_vec(1) = (n * cum_sum_sq)/  (cum_sum^2) ;

new_mat = image_2d ;

for q = 2:length(lacunarity_vec)
   
   disp_string = [num2str(q) ' ' datestr(now)] ;
   disp(disp_string) ;
   box_size = box_size_vec(q) ;
   old_mat = new_mat ;
   num_rows_new = num_rows - box_size + 1 ;
   num_cols_new = num_cols - box_size + 1 ;
   new_mat = zeros(num_rows_new, num_cols_new) ;
   n = num_rows_new * num_cols_new ;
   skip = 2^(q - 2) ;
   
   cum_sum = 0 ;
   cum_sum_sq = 0 ;
   
   for new_row = 1:num_rows_new
      
      for new_col = 1:num_cols_new
         
         t = sum(sum(old_mat([new_row new_row+skip], [new_col new_col+skip]))) ;
         cum_sum = cum_sum + t ;
         cum_sum_sq = cum_sum_sq + (t^2) ;
         
         new_mat(new_row, new_col) = t ;
         
      end
      
   end
   
   if (cum_sum ~= 0)
      lacunarity_vec(q) = (n * cum_sum_sq) / (cum_sum^2) ;
   else
      lacunarity_vec(q) = -9999 ;
   end

end