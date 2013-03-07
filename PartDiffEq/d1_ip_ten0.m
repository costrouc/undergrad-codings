function [localmat] = d1_ip_ten0(isub, quad_rul, ...
   scal_fun, ten0_type) 
%
% This function computes, for subinterval isub, the integrals
% of the (scalar) ten0 basis functions multiplied by the
% scalar function  scal_fun.
% The vector of values is returned in localmat.
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
[ten0, Gradten0] = feval(ten0_type, quad_pts) ;
nbas0 = size(ten0,1) ;

% Adjust points and weights to account for size of true triangle.
quad_pts = xleft + hsub* quad_pts ;
quad_wghts = hsub * quad_wghts ;

% Evaluate the scalar multiplier at the quadrature points.
sfun_vals = feval(scal_fun, quad_pts, isub) ;

% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   ten0(:,iq) = quad_wghts(iq) * ten0(:,iq) ;  
end


% Note: localmat is a vector
localmat = ten0 * sfun_vals.' ; 


