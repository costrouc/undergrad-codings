function [localmat] = d1_ip_Dten0_ten0(isub, quad_rul, ...
   scal_fun, ten0a_type, ten0b_type) 
%
% This function computes, for subinterval isub, the integrals
% of the (scalar) derivative of ten0a basis functions times 
% the (scalar) ten0b basis functions multiplied by the scalar
% function scal_fun. The matrix of values is returned in localmat.
%  
%


%%%%%%%%%%%%%%%%%%%%%% Global Variables %%%%%%%%%%%%%%%%%%%
global xpts nnds
global Global_r  Global_s  Global_u
global rad_bas_type  str_bas_type  vel_bas_type

% Description of subinterval.
xleft = xpts(isub) ;
xright = xpts(isub + 1) ;
hsub = xright - xleft ;

% Evaluation of quadrature points and quadrature weights.
[quad_pts, quad_wghts] = feval(quad_rul) ;
nqpts = size(quad_pts,1) ;

% Evaluate Basis Functions and their Gradients at quad. points.
[ten0a, Gradten0a] = feval(ten0a_type, quad_pts) ;
nbas0a = size(ten0a,1) ;

[ten0b, Gradten0b] = feval(ten0b_type, quad_pts) ;
nbas0b = size(ten0b,1) ;

% Do appropriate scaling to get the true derivatives.
Gradten0b = Gradten0b / hsub ;

% Adjust points and weights to account for size of true triangle.
quad_pts = xleft + hsub* quad_pts ;
quad_wghts = hsub * quad_wghts ;

% Evaluate the scalar multiplier at the quadrature points.
sfun_vals = feval(scal_fun, quad_pts, isub) ;


% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   ten0a(:,iq) = quad_wghts(iq) * sfun_vals(iq) * ten0a(:,iq) ;  
end

mat1 = ten0a * Gradten0b.' ;

localmat = [ mat1 ] ; 


