function [localmat] = d1_ip_ten0_ten0_Dten0(isub, quad_rul, ...
   scal_fun, ten0a_type, ten0b_type) 
%
% This function computes, for subinterval isub, the integrals
% of the (scalar) derivative of ten0a basis functions times the 
% (scalar) ten0b basis functions the (scalar) ten0b basis functions 
% multiplied by the scalar function scal_fun. The matrix of values is 
% returned in localmat.
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

% Evaluate Non-linear function at known location 
[ten0c, Gradten0c] = feval(rad_bas_type, quad_pts) ;
nbas0c = size(ten0c,1) ;

% Determine weights for Non-Linear function
% Identify the global unknown coefficients
if strcmp(rad_bas_type, 'd1_CtsLin') == 1
     local_weights = Global_r(isub:isub + 1) ; 
elseif strcmp(rad_bas_type, 'd1_CtsQuad') == 1
     local_weights = Global_r(2*isub - 1:2*isub + 1) ;     
elseif strcmp(rad_bas_type, 'd1_CtsCub') == 1
     local_weights = Global_r(3*isub - 2:3*isub + 1) ;
end

% Do appropriate scaling to get the true derivatives.
Gradten0b = Gradten0b / hsub ;

% Adjust points and weights to account for size of true triangle.
quad_pts = xleft + hsub* quad_pts ;
quad_wghts = hsub * quad_wghts ;

% Evaluate the scalar multiplier at the quadrature points.
sfun_vals = feval(scal_fun, quad_pts, isub) ;

% Evaluate the known non-linear basis function at quadrature points.
ten0c = local_weights' * ten0c;

% Now to do the evaluations of the integrals.
for iq = 1:nqpts
   ten0a(:,iq) = quad_wghts(iq) * sfun_vals(iq) * ten0a(:,iq) * ten0c(:,iq);  
end

mat1 = ten0a * Gradten0b.' ;

localmat = [ mat1 ] ; 


