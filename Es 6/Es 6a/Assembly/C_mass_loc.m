function [M_loc] = C_mass_loc(dphiq,w_1D,nln,BJ)
%% [M_loc] = C_mass_loc(dphiq,w_1D,nln,BJ)
%==========================================================================
% Build the local mass matrix for the term grad(u)grad(v)
%==========================================================================
%    called in C_matrix1D.m
%
%    INPUT:
%          dphiq       : (array real) vecotr containing the evaluation of the
%                         i-th basis function on the quadrature nodes
%          w_1D        : (array real) quadrature weights
%          nln         : (integer) number of local unknowns
%          BJ          : (array real) Jacobian of the map 
%
%    OUTPUT:
%          M_loc       :  (array real) Local stiffness matrix


M_loc = zeros(nln,nln);

%% General implementation -- to be used with general finite element spaces
for i=1:nln
    for j=1:nln
        for k=1:length(w_1D)
            Binv = 1./BJ;    % inverse
            Jdet = BJ;       % determinant 
            M_loc(i,j) = M_loc(i,j) + (Jdet.*w_1D(k)) .*  dphiq(1,k,i).* dphiq(1,k,j);
        end
    end
end