function [Amat, bvec] = ScaleM(Amat, bvec)

% this function does a row scale of the  matrix Amat
% and the associated right hand side vector bvec.

[Nukn,ncol] = size(Amat) ;

for i = 1:Nukn
    scale = norm(Amat(i,:), inf) ;
    Amat(i,:) = Amat(i,:) / scale ;
    bvec(i,:) = bvec(i,:) / scale ;
end
