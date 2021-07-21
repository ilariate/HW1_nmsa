function [M,A,f]=C_matrix1D(Dati,femregion)
%% [M,A,f] = C_matrix1D(Dati,femregion)
%==========================================================================
% Assembly of the stiffness matrix A, the mass matrix M and rhs f
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%          femregion   : (struct)  see C_create_femregion.m
%
%    OUTPUT:
%          A           : (sparse(ndof,ndof) real) stiffnes matrix
%          M           : (sparse(ndof,ndof) real) mass matrix
%          f           : (sparse(ndof,1) real) rhs vector


addpath FESpace
addpath Assembly

fprintf('============================================================\n')
fprintf('Assembling matrices and right hand side ... \n');
fprintf('============================================================\n')


% connectivity infos
ndof         = femregion.ndof; % degrees of freedom
nln          = femregion.nln;  % local degrees of freedom
ne           = femregion.ne;   % number of elements
connectivity = femregion.connectivity; % connectivity matrix


% shape functions
[basis] = C_shape_basis(Dati);

% quadrature nodes and weights for integrals
[nodes_1D, w_1D] = C_quadrature(Dati);

% evaluation of shape bases 
[dphiq,Grad] = C_evalshape(basis,nodes_1D);


% Assembly begin ...
A = sparse(ndof,ndof);  % Global Stiffness matrix
M = sparse(ndof,ndof);
f = sparse(ndof,1);     % Global Load vector

for ie = 1 : ne
     
    % Local to global map --> To be used in the assembly phase
    iglo = connectivity(1:nln,ie);
  
    
    [BJ, pphys_1D] = C_get_Jacobian(femregion.coord(iglo,:), nodes_1D);
    % BJ        = Jacobian of the elemental map 
    % pphys_2D = vertex coordinates in the physical domain 
   
    %=============================================================%
    % STIFFNESS MATRIX
    %=============================================================%
    
    % Local stiffness matrix 
    [A_loc] = C_lap_loc(Grad,w_1D,nln,BJ);

    % Assembly phase for stiffness matrix
    A(iglo,iglo) = A(iglo,iglo) + Dati.mu*A_loc; 
    
    %=============================================================%
    % MASS MATRIX
    %=============================================================%
    
    % Local stiffness matrix 
    [M_loc] = C_mass_loc(dphiq,w_1D,nln,BJ);

    % Assembly phase for stiffness matrix
    M(iglo,iglo) = M(iglo,iglo) + M_loc; 
    
    %==============================================
    % FORCING TERM --RHS
    %==============================================

    % Local load vector
    [load] = C_loc_rhs1D(Dati.force,dphiq,BJ,w_1D,pphys_1D,nln);    

    % Assembly phase for the load vector
    f(iglo) = f(iglo) + load;

end
