function E = Input_2D_E(u,n,k_mode, node_1, node_2)

if k_mode == 1
    E = 70e+9;
elseif k_mode == 2
    
    n1 = node_1 *2;
    n2 = node_2 *2;

    n1x = n1-1;
    n2x = n2-1;
    
    uy = u([n1x, n2x]) % ydisplacement
 
    strainy = uy/2;
     strain_avg = mean(strainy,'all');
 
    %assumption that overall length is not changed massively 
    E = norm(70e+9*strain_avg^(n-1)); % E is now not a single variable, but a vector the same size as the number of nodes

%  Returns the Young's modulus for subsequent displacement 
end