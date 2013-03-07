function [quad_pts, quad_wghts] = d1_quad_5
%
% This function contains weights and quadrature points
% which are exact for polynomials of degree five for evaluting
% the line integral from (0,0) to (1,0).
%

quad_pts = zeros(3,1) ;
quad_wghts = zeros(3,1) ;

quad_pts(1) = (1 - sqrt(3/5))/2 ;
quad_pts(2) = (1 + sqrt(3/5))/2 ;
quad_pts(3) = 1/2 ;


quad_wghts(1:2) = 5/18* ones(1,2) ;
quad_wghts(3) = 8/18 ;

