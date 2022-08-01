%  Simple finite element analysis code
%
%  ME4291 FEA
%  Elasticity 2D Triangle Elements

%  Load node coordinates and element connectivity
%  Always load 1 row mesh model only!! (We will reduce to 1D beam element)
nodes = load('Input_node_coord_10.txt')
elems = load('Input_elem_connect_4.txt')

%  Determine the number of nodes
node_size = size(nodes);
num_of_nodes = node_size(1)/2 % reduced to 1D 
nodes = nodes(1:num_of_nodes) % trim nodes (only x axis)

%  Allocate space for [K] matrix and [F] vector
K = zeros(2*num_of_nodes, 2*num_of_nodes);
F = zeros(2*num_of_nodes, 1);

%  Determine the number of elements
elem_size = size(elems);
num_of_elem = elem_size(1);

u=0;
Nmax = 100;
N_iter = 0;
tol = 0.001;
k_mode = 1;
percentage_diff = 1;
n = 0.9; %non linear iteration power 

while(percentage_diff > tol) && (N_iter < Nmax)  
    
    %initialising first set of displacements using linear E 
    if k_mode == 1
        
        %  Loop through every element
        for j = 1:num_of_elem
            %  Obtain global node numbers
            node_1 = elems(j,1); %1, 2, ..., 4
            node_2 = elems(j,2); %2, 3, ..., 5

            %   Form stiffness beam matrix for uniform beam
            Ke = Sub_uniform_beam_stiffness( nodes(node_1), nodes(node_2), j, u, n, k_mode, elems); % NOTE: only 1 row mesh model should be used for beam model
        
            % Assemble into global matrix
            global_dof = [2*node_1-1, 2*node_1, 2*node_2-1, 2*node_2]; % 1, 2, 3, 4/ 3, 4, 5, 6/ ...
            K( global_dof, global_dof ) = K( global_dof, global_dof ) + Ke;

        end
        %  End of for-loop for each element 

        %  Apply force boundary condition (single downward nodal force at coord(5.0,2.0))
        force_data = load('Input_bc_nodalforce.txt'); % coordinate of nodal force to be applied
        force_data_size = size(force_data);
        num_of_force = force_data_size(1);

        % Loop through all coords
        for j = 1:num_of_force
           % find node to be 
           for i = 1:num_of_nodes
            if force_data(j,1) == nodes(i)
                Q = [-10e3; 0]; % constand nodal force, nom moment (Fy, M)
            else
                Q = [0;0]; % Fy, M
            end

            % Assemble Q into F
            global_dof = [2*i-1, 2*i];
            F( global_dof ) = F( global_dof ) + Q;
            end 
        end


        %  Apply disp boundary condition
        F_copy = F;
        K_copy = K;

        disp_data = load('Input_bc_disp_1.txt');
        disp_nodes = disp_data(:,1);
        disp_dof = disp_data(:,2); % disp or rotation
        disp_values = disp_data(:,3);
        num_of_disp = length(disp_values)/2; % reducing to 1D

        for j = 1:num_of_disp
            % Determine global dof
            if (disp_dof(j) == 1)
                global_dof = 2*disp_nodes(j) - 1; % displacement
            else 
                global_dof = 2*disp_nodes(j) ;  % rotation
            end   
            % the convenient way to solve for BC
            F(global_dof) = disp_values(j);
            K(global_dof,:) = zeros(1,2*num_of_nodes);
            K(global_dof,global_dof) = 1;    
        end

        %  Solve for the displacements
        u = K\F;
        %initialising variable for previous u 
        u_copy = u;
        

        % Calculate reaction force
        R = K_copy * u - F_copy;
        
        %initiate iteration as displacement changes
        k_mode = 2; 
    
    % iterating with changes in displacement 
    elseif k_mode ==2
    
        %  Loop through every element
        for j = 1:num_of_elem
            %  Obtain global node numbers
            node_1 = elems(j,1); %1, 2, ..., 4
            node_2 = elems(j,2); %2, 3, ..., 5

            %   Form stiffness beam matrix for uniform beam
            Ke = Sub_uniform_beam_stiffness( nodes(node_1), nodes(node_2), j, u, n, k_mode, elems); % NOTE: only 1 row mesh model should be used for beam model

            % Assemble into global matrix
            global_dof = [2*node_1-1, 2*node_1, 2*node_2-1, 2*node_2]; % 1, 2, 3, 4/ 3, 4, 5, 6/ ...
            K( global_dof, global_dof ) = K( global_dof, global_dof ) + Ke;
        end
        %  End of for-loop for each element 

        %  Apply force boundary condition (single downward nodal force at coord(5.0,2.0))
        force_data = load('Input_bc_nodalforce.txt'); % coordinate of nodal force to be applied
        force_data_size = size(force_data);
        num_of_force = force_data_size(1);

        % Loop through all coords
        for j = 1:num_of_force
           % find node to be 
           for i = 1:num_of_nodes
            if force_data(j,1) == nodes(i)
                Q = [-10e3; 0]; % constand nodal force, nom moment (Fy, M)
            else
                Q = [0;0]; % Fy, M
            end

            % Assemble Q into F
            global_dof = [2*i-1, 2*i];
            F( global_dof ) = F( global_dof ) + Q;
            end 
        end


        %  Apply disp boundary condition
        F_copy = F;
        K_copy = K;

        disp_data = load('Input_bc_disp_1.txt');
        disp_nodes = disp_data(:,1);
        disp_dof = disp_data(:,2); % disp or rotation
        disp_values = disp_data(:,3);
        num_of_disp = length(disp_values)/2; % reducing to 1D

        for j = 1:num_of_disp
            % Determine global dof
            if (disp_dof(j) == 1)
                global_dof = 2*disp_nodes(j) - 1; % displacement
            else 
                global_dof = 2*disp_nodes(j) ;  % rotation
            end   
            % the convenient way to solve for BC
            F(global_dof) = disp_values(j);
            K(global_dof,:) = zeros(1,2*num_of_nodes);
            K(global_dof,global_dof) = 1;    
        end

        %  Solve for the displacements
        u = K\F;
        
        % Checking for deviation
        u_diff = u - u_copy;
        percentage_diff = norm(u_diff)/norm(u);                
        
        %initialising variable for previous u 
        u_copy = u;
        

        % Calculate reaction force
        R = K_copy * u - F_copy;
        
        % updating N_iter
        N_iter = N_iter + 1
        per_diff = percentage_diff *100
    end
    
end

u
R





% displacement in y due to poisson ratio. bottom is set to bc (disp) = 0,
% therefore the top part of the elements will experience the y negative
% displacement 

% can keep copies of the matrix before setting to 0 if we want to find
% the force vector
