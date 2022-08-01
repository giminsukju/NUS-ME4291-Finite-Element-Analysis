function Ke = Sub_uniform_beam_stiffness(node_1_x_coord, node_2_x_coord, j, u, n, k_mode, elems)
    b = 0.01; %  thickness of the beam, assumbed to be constant
    h0 = 2.0; % height of the beam element
    L = node_2_x_coord - node_1_x_coord; % length of the beam element

    % Eqn 540 from the lecture note
    NiNj_integrated = [12,   6*L,    -12,    6*L;k
                       6*L,  4*L*L,  -6*L,   2*L*L;
                       -12,  -6*L,   12,     -6*L;
                       6*L,  2*L*L,  -6*L,   4*L*L];

    %Load Node Numbers
    node_1 = elems(j,1);
    node_2 = elems(j,2);
%    node_3 = elems(j,3);
%    node_4 = elems(j,4);

    E = Input_2D_E(u,n,k_mode, node_1, node_2); % to use only node_1 and node_2

    Ke = b*h0/(12*L^3)*E*NiNj_integrated;
end
