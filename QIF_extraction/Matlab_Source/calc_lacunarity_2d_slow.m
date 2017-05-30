function [lacunarity_vec, box_size_vec] = calc_lacunarity_2d_slow(image_2d)

temp_a = size(image_2d) ;
num_rows = temp_a(1) ;
num_cols = temp_a(2) ;

temp_b = floor(log(min(temp_a)) / log(2)) ;
box_size_vec = 2.^(0:temp_b) ;
lacunarity_vec = zeros(size(box_size_vec)) ;

for q = 1:length(lacunarity_vec)
   
   disp_string = [num2str(q) ' ' datestr(now)] ;
   disp(disp_string) ;
   
   box_size = box_size_vec(q) ;
   box_size_m1 = box_size - 1 ;
   cum_sum = 0 ;
   cum_sum_sq = 0 ;
   n = (num_rows - box_size + 1) * (num_cols - box_size + 1) ;
   
   for row = 1:(num_rows - box_size + 1)
      
      for col = 1:(num_cols - box_size + 1)
         
         t = sum(sum(image_2d(row:row + box_size_m1,col:col + box_size_m1))) ;
         cum_sum = cum_sum + t ;
         cum_sum_sq = cum_sum_sq + (t^2) ;
         
      end
      
   end
   
   if (cum_sum ~= 0)
      lacunarity_vec(q) = (n * cum_sum_sq) / (cum_sum^2) ;
   else
      lacunarity_vec(q) = -9999 ;
   end

end