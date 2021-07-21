%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test 
%                'Domain',            % set the domain [x1,x2]
%                'exact_sol',         % set the exact solution
%                'force',             % set the forcing term
%                'grad_exact_1',      % set the first componenet of the gradient of the exact solution
%                'fem',               % set finite element space
%                'nqn_1D',            % number of quadrature nodes for integrals over lines
%                'refinement_vector', % set the level of refinement for the grid
%                'visual_graph',      % if you want to display the graphical results ['Y','N']
%                'print_out',         % if you want to print out the results ['Y','N']
%                'plot_errors'        % you want to print the computed errors ['Y','N']
% 
%========================================================================================================

function [Dati]=C_dati(test)

if test=='Test1'
Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,1],...                          
               ... % Domain bounds       
               'mu',               1, ...
               ... % Diffusive term ...
               'omega',               1, ...  
               'n',               1, ...
               'exact_sol',        'sin(2*pi*x)',...      
               ... % Definition of exact solution
               'force',            'sin(2*pi*x) - 4*pi^2*sin(2*pi*x)',...  
               ... % Forcing term
               'grad_exact',     '2*pi*cos(2*pi*x)',...    
               ... % Definition of exact gradient (x comp) 
               'fem',              'P1',...
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % Triangular Structured mesh
               'refinement_vector', [2,3,4,5],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );

end



