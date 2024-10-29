function [B, detJ] = B_matrix(coordinate, xi, eta)
    % Questa function mi calcola lo Jacobiano e la matrice B  matrix dati 
    % xi, eta e le coordinate dell'elemento
    % coordinate: matrice 4x2  con le coordinate nodali dell'elemento
    % xi, eta: Gauss point in natural coordinates
    
%% Derivate delle funzioni di forma rispetto a xi e eta
naturalDerivatives = 1/4 * [-(1-eta), -(1-xi); 
                            1-eta, -(1+xi); 
                            1+eta,  1+xi; 
                            -(1+eta), 1-xi];

%% Matrice Jacobiana
 
    % J = naturalDerivatives' * coordinate; % versione dinamica
    J = [1/2 0; 0 1/2];     %YG
    detJ = det(J);
    invJ = J\eye(size(J));

        %% DEBUG 
        % disp('Determinante della Jacobiana:');
        % disp(detJ);
%% Derivate rispetto x e y

derivativesXY = zeros(4, 2);
for i = 1:4
    derivativesXY(i, :) = invJ * naturalDerivatives(i, :)';
end

%% B matrix
B = zeros(3, 8);  

for i = 1:4
    B(1, 2*i-1) = derivativesXY(i, 1);  % Derivata rispetto a x
    B(2, 2*i)   = derivativesXY(i, 2);  % Derivata rispetto a y
    B(3, 2*i-1) = derivativesXY(i, 2);  % Derivata rispetto a y (per gamma_xy)
    B(3, 2*i)   = derivativesXY(i, 1);  % Derivata rispetto a x (per gamma_xy)
end
% disp('Matrice B:');
% disp(B);
    
end

%%
 %shape = 1/4*[(1-xi)*(1-eta); (1+xi)*(1-eta);
 % (1+xi)*(1+eta); (1-xi)*(1+eta)];
 
 % Funzioni di forma in coordinate naturali
    % dN_dxi = 0.25 * [-1 + eta, 1 - eta, 1 + eta, -1 - eta];
    % dN_deta = 0.25 * [-1 + xi, -1 - xi, 1 + xi, 1 - xi];
    % 
    % naturalDerivatives = 1/4*[-(1-eta),-(1-xi); 1-eta,-(1+xi);1+eta ,1+xi;-(1+eta), 1-xi];

%% Matrice Jacobiana
    % % J = [dN_dxi; dN_deta] * coordinate;
    % J = [1/2 0; 0 1/2];                 %YG
    % % invJ = inv(J);                      % This not a inverse; wrong! 
    % detJ = det(J);
    % invJ = J\eye(size(J));
    % 
    
%% Funzioni di forma nelle coordinate x e y
    % dN_dx = [invJ(1, 1) * dN_dxi , invJ(1, 2) * dN_deta];
    % dN_dy = [invJ(2, 1) * dN_dxi , invJ(2, 2) * dN_deta];
    % 

    %% Matrice B
    % B = zeros(3, 8);
    % for i = 1:4
    %     B(1, 2*i-1) = dN_dx(i);
    %     B(2, 2*i) = dN_dy(i);
    %     B(3, 2*i-1) = dN_dy(i);
    %     B(3, 2*i) = dN_dx(i);
    % end

