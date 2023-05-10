function bands = band_structure(k_x_min, k_x_max, n_points)
    % Define physical constants
    m_0 = 5.6667e-12;
    h_bar = 6.582e-16;

    % Define superconductor constants for GaAs
    a = 5.6532e-10;
    gamma_1 = 6.79;
    gamma_2 = 1.924;
    gamma_3 = 2.782;
    m_e = 0.0665;
    E_p = 28.8;
    E_g = 1.424;
    E_v = 0; 
    delta = 0.341;
    P = sqrt(h_bar^2*E_p / (2*m_0));

    % Define Hamiltonian parameters
    A = h_bar^2 / (2*m_e*m_0) - P^2 * (E_g + 2/3 * delta) / (E_g * (E_g + delta));
    L = - h_bar^2 * (1 + gamma_1 + 4*gamma_2) / (2*m_0) + P^2/E_g;
    M = - h_bar^2 * (1 + gamma_1 - 2*gamma_2) / (2*m_0);
    N = - 3*h_bar^2 * gamma_3 / m_0 + P^2 / E_g;
    
    % Define k-vector mesh
    k_x_mesh = linspace(k_x_min, k_x_max, n_points);
    k_y = 0;
    k_z = 0;

    % Solve for eigenenergies at every k-value
    bands = zeros(n_points, 8);
    for i = 1:n_points
        k_x = k_x_mesh(i);
        k = sqrt(k_x^2 + k_y^2 + k_z^2);

        % Define submatrices and hamiltonian
        G_1 = [E_v + E_g, 1j*P*k_x, 1j*P*k_y, 1j*P*k_z;
            -1j*P*k_x, E_v - delta/3, 0, 0;
            -1j*P*k_y, 0, E_v - delta/3, 0;
            -1j*P*k_z, 0, 0, E_v - delta/3];

        G_2 = [
            A*k^2, 0, 0, 0;
            0, L*k_x^2 + M*(k_y^2 + k_z^2), N*k_x*k_y, N*k_x*k_z;
            0, N*k_x*k_y, L*k_y^2 + M*(k_x^2 + k_z^2), N*k_y*k_z;
            0, N*k_x*k_z, N*k_y*k_z, L*k_z^2 + M*(k_x^2 + k_y^2)];

        G_so = -delta/3 * [0, 0, 0, 0;
            0, 0, 1j, 0;
            0, -1j, 0, 0;
            0, 0, 0, 0];

        GAMMA = -delta/3 * [0, 0, 0, 0;
            0, 0, 0, -1;
            0, 0, 0, 1j;
            0, 1, -1j, 0];

        G = G_1 + G_2 + G_so;

        H = [G, GAMMA;
            -conj(GAMMA), conj(G)];

        % Solve eigenvalue problem
        bands(i, :) = eig(H);
    end
end

