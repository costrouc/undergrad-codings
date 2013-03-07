function [CQuadVal, DerivCQuadVal] = d1_CtsQuad(quad_pts)
%
% This function computes the values of the continuous
% quadratic basis functions, and of its gradient, at
% the quadrature points quad_pts --- on the reference triangle.
%

nqpts = size(quad_pts,1) ;

CQuadVal(1,:) = 2*(quad_pts(1:nqpts)' - 0.5).* (quad_pts(1:nqpts)' - 1) ;
CQuadVal(2,:) = -4*quad_pts(1:nqpts)'.* (quad_pts(1:nqpts)' - 1) ;
CQuadVal(3,:) = 2*quad_pts(1:nqpts)'.* (quad_pts(1:nqpts)' - 0.5) ;

DerivCQuadVal(1,:) = 4*quad_pts(1:nqpts)' - 3 ;
DerivCQuadVal(2,:) = 4 - 8*quad_pts(1:nqpts)' ;
DerivCQuadVal(3,:) = 4*quad_pts(1:nqpts)' - 1 ;
