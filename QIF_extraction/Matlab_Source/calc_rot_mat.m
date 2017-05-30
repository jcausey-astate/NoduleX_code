function [rot_mat] = calc_rot_mat(alpha, beta, gamma)

%calc_rot_mat.m
%  This function defines the 3D rotation matrix corresponding to Euler
%  angles of rotation alpha, beta, and gamma, which are the input
%  arguments.
   
   rx = [1           0           0           ; ...
         0           cos(alpha) -sin(alpha)  ; ...
         0           sin(alpha)  cos(alpha)] ;
      
   ry = [cos(beta)   0          -sin(beta)   ; ...
         0           1           0           ; ...
         sin(beta)   0           cos(beta) ] ;
      
   rz = [cos(gamma) -sin(gamma)  0           ; ...
         sin(gamma)  cos(gamma)  0           ; ...
         0           0           1         ] ;
      
   rot_mat = rx * ry * rz ;

end

