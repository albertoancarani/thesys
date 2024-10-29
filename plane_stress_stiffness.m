function [Ke,B,D] = plane_stress_stiffness()
    E = 1;    % Modulo di Young
    nu = 0.3;     % Coefficiente di Poisson
    % coordinate: coordinate dei 4 nodi dell'elemento (4x2 matrix)
    coordinate=[-1, -1;  % Nodo 1
                 1, -1;  % Nodo 2
                 1,  1;  % Nodo 3
                -1,  1]; % Nodo 4
    % Matrice costitutiva D
    D = (E / (1 - nu^2)) * [1, nu, 0;
                            nu, 1, 0;
                            0,  0, (1 - nu) / 2];

    % Punti e pesi per la quadratura di Gauss 
    gauss_points = [-1 1] / sqrt(3);
    weights = [1, 1];
    
    Ke = zeros(8, 8);  
    
    for xi = 1:2
        for eta = 1:2
            %Calcolo le funzioni di forma nei punti di Gauss
            [B, detJ] = B_matrix(coordinate, gauss_points(xi), gauss_points(eta));
            
            % Contributo della matrice di rigidezza nel punto di Gauss
            % considerato
            Ke = Ke + B' * D * B * detJ * weights(xi) * weights(eta);
        end
    end
end

%%This is good! Verified from Top88.m
