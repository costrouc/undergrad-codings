function [localmat] = Bfun(x_pts, isub) 
%
% This function represents the diffusion coefficient function.
%  
%  
%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global xpts nnds
global Global_r  Global_s  Global_u
global rad_bas_type  str_bas_type  vel_bas_type
%
%

% Test example: rtrue = 6x^2 + 4x + 2
%localmat = (3* x_pts + 1).' ;

localmat = (2+0.*x_pts).';


