%=======================================================================================================
% This contain all the information for running main
% TEMPLATE OF THE STRUCT DATI
%=======================================================================================================
%
%  DATI= struct( 'name',              % set the name of the test 
%                'Domain',            % set the domain [x1,x2]
%                'c2'                 % c^2 wave speed
%                'T'                  % final time
%                'dt'                 % time step 
%                'u0',                % Initial condition u               
%                'v0',                % Initial condition  du/dt        
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
               'c2',               1, ...
               ... % Diffusive term ...
               'bc',               'DD', ...
               ... % boundary conditions ...             
               'T',               4, ...
               ... % Final time ...
               'dt',            0.01, ...
               ... % Time step
               'u0',      'sin(2*pi*x)', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'sin(2*pi*x).*cos(t)',...      
               ... % Definition of exact solution
               'grad_exact',     '2*pi*cos(2*pi*x).*cos(t)',...    
               ... % du/dx 
               'force',           '(-sin(2*pi*x) + 4*pi^2*sin(2*pi*x)).*cos(t).*(t>0)',...  
               ... % Forcing term
               'g1',           '1',...  
               ... % boundary condition in x=0 
               'g2',           '2',...  
               ... % boundary condition in x=L
               'fem',              'P1',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [2,3,4,5],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );
elseif test=='Test2'
       Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [0,2],...                          
               ... % Domain bounds       
               'c2',               1, ...
               ... % Diffusive term ...
               'bc',               'PP', ...
               ... % boundary conditions ...             
               'T',               2.5, ...
               ... % Final time ...
               'dt',            0.005, ...
               ... % Time step
               'u0',      '0.*x', ...
               ... % Initial condition u               
               'v0',      '0.*x', ...
               ... % Initial condition  du/dt        
               'exact_sol',       '0.*x.*t',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...    
               ... % du/dx 
               'force',           '100*exp(-100.*(x-0.25).^2).*(t<=0.25)',...  
               ... % Forcing term
               'g1',           '0.*t',...  
               ... % boundary condition in x=0 
               'g2',           '0.*t',...  
               ... % boundary condition in x=L
               'fem',              'P1',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [2,3,4,5],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );
           
elseif test=='Test3'
       Dati = struct( 'name',             test,...
               ... % Test name
               'domain',           [-1,1],...                          
               ... % Domain bounds       
               'c2',               1, ...
               ... % Diffusive term ...
               'bc',               'PP', ...
               ... % boundary conditions ...             
               'T',               5, ...
               ... % Final time ...
               'dt',            0.005, ...
               ... % Time step
               'u0',      'cos(4*pi*x)', ...
               ... % Initial condition u               
               'v0',      '-4*pi*sin(4*pi*x)', ...
               ... % Initial condition  du/dt        
               'exact_sol',       'cos(4*pi*x+4*pi*t)',...      
               ... % Definition of exact solution
               'grad_exact',     '0.*x.*t',...    
               ... % du/dx 
               'force',           '0.*x.*t',...  
               ... % Forcing term
               'g1',           '0.*t',...  
               ... % boundary condition in x=0 
               'g2',           '0.*t',...  
               ... % boundary condition in x=L
               'fem',              'P1',...         
               ... % P1-fe
               'nqn_1D',            2,...           
               ... % Number of quad. points per element
               'MeshType',         'TS', ...        
               ... % uniform regular mesh
               'refinement_vector', [2,3,4,5],...   
               ... % Refinement levels for  the error analysis
               'visual_graph',      'Y',...
               ... % Visualization of the solution
               'plot_errors',       'Y' ...
               ...% Compute Errors 
               );           

end



