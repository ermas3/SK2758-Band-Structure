function band_structure = calculate_band_structure(A, L, M, N, E_g, delta, E_v, P, n_points, k_x_mesh)
% CALCULATE_BAND_STRUCTURE calculates the band structure of a crystal using the given Hamiltonian parameters
%
% INPUTS:
% A: constant
% L: constant
% M: constant
% N: constant
% E_g: band gap energy
% delta: spin-orbit splitting energy
% E_v: valence band energy
% P: spin-orbit interaction constant
% n_points: number of points in the k_x_mesh
% k_x_mesh: a vector of k_x values
%
% OUTPUTS:
% band_structure: an array of size (8, n_points) that contains the energy eigenvalues for each point in k_x_mesh

    % Define k-vector mesh
    k_y = 0;
    k_z = 0;

    band_structure = zeros(8, n_points);

    for i = 1:n_points
        k_x = k_x_mesh(i);
        k = sqrt(k_x^2 + k_y^2 + k_z^2);

        % Define submatrices and hamiltonian
        G_1 = [            E_v + E_g, 1i*P*k_x, 1i*P*k_y, 1i*P*k_z;            -1i*P*k_x, E_v - delta/3, 0, 0;            -1i*P*k_y, 0, E_v - delta/3, 0;            -1i*P*k_z, 0, 0, E_v - delta/3];
        G_2 = [            A*k^2, 0, 0, 0;            0, L*k_x^2 + M*(k_y^2 + k_z^2), N*k_x*k_y, N*k_x*k_z;            0, N*k_x*k_y, L*k_y^2 + M*(k_x^2 + k_z^2), N*k_y*k_z;            0, N*k_x*k_z, N*k_y*k_z, L*k_z^2 + M*(k_x^2 + k_y^2)];

        G_so = -delta/3 * [0, 0, 0, 0;            0, 0, 1i, 0;            0, -1i, 0, 0;            0, 0, 0, 0];

        GAMMA = -delta/3 * [0, 0, 0, 0;            0, 0, 0, -1;            0, 0, 0, 1i;            0, 1, -1i, 0];

        G = G_1 + G_2 + G_so;

        H = [G, GAMMA;            -conj(GAMMA), conj(G)];

        % Solve eigenvalue problem
        energy_values = eig(H);

        band_structure(:, i) = energy_values;
    end
end