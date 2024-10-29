function [pnorm, pnorm_sen, MISES] = Stress_2D_Sensitivity_Comp(x, nelx, nely, pl, q, p)

      % cantilever beam
      
    %% Parametri del materiale

    [KE, B, D] = plane_stress_stiffness();
    
    % MATERIAL PROPERTIES
    E0 = 1;         % Modulo di Young del materiale solido
    Emin = 1e-9;    % Modulo di Young del materiale "vuoto"
    nele = nelx * nely;             % Numero totale di elementi      
    ndof = 2*(nelx+1)*(nely+1);     % NUMERO TOTALE DI GRADI DI LIBERTÀ

   %% MATRICE DI CONNETTIVITA edofMat

    nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
    edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
    edofMat = repmat(edofVec,1,8)+ ... 
        repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
   
    %% CARICHI && VINCOLI YG


    % Cantilever case boundary condtion 
    % fixeddofs = 1:1:2*(nely+1);    
    % CARICO                                      
    % loaddofs = 2*(nelx +1)*(nely+1);
    % F = sparse(loaddofs,1,-1,ndof,1);  % Forza verso il basso (negativa)

    %% CARICHI && VINCOLI 

    leftBoundaryNodes = nodenrs(:,1);  % Nodi lungo il lato sinistro

    %  Vincolo completo a sinistra (bloccato in x: x=0, y=free )
    fixed_dofs_x = 2*leftBoundaryNodes - 1;  % DOF in x (movimento orizzontale bloccato)
    fixed_dofs_y = 2*leftBoundaryNodes;      % DOF in y (movimento verticale bloccato)

    % Unisci i DOF bloccati in un unico vettore
    fixeddofs = [fixed_dofs_x; fixed_dofs_y];
  
    bottom_right_node = nodenrs(end, end);  % Nodo inferiore destro
    F = sparse(2*bottom_right_node, 1, -1, 2*(1+nely)*(1+nelx), 1);
    
    %% MATRICE DI RIGIDEZZA GLOBALE
    % iK = reshape(kron(edofMat, ones(8, 1))', 8*8*nele, 1);
    % jK = reshape(kron(edofMat, ones(1, 8))', 8*8*nele, 1);
    % sK = reshape(KE(:)*(Emin + x(:)'.^pl*(E0-Emin)), 8*8*nele, 1);
    % K = sparse(iK, jK, sK); 
    % K = (K + K') / 2; %posso calcolarla così perchè è simmetrica

    % From top88.m
    iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
    jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

    sK = reshape(KE(:)*(Emin+x(:)'.^pl*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;

    %% RISOLUZIONE DEL SISTEMA

    % GRADI DI LIBERTÀ LIBERI
    freedofs = setdiff(1:ndof, fixeddofs); % Gradi di libertà liberi

    U = zeros(ndof, 1);
    % K(freedofs, freedofs) = K(freedofs, freedofs) + 1e-6 * speye(length(freedofs));
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);   %3D U(freedofs,:) = K(freedofs, freedofs) \ F(freedofs,:);
    %% STRESS DI VON MISES

    MISES = zeros(nele, 1);  % Stress di Von Mises
    S = zeros(nele, 3);  % Matrice di stress rilassato

    for i = 1:nele
        temp = x(i)^q * (D * B * U(edofMat(i, :)))';
        S(i, :) = temp;
        % Calcolo dello stress di Von Mises
        MISES(i) = real(sqrt(temp(1)^2 - temp(1)*temp(2) + temp(2)^2 + 3*temp(3)^2)); %% QUESTA FORMULA è DA CONTROLLARE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % disp(MISES(i))
        % MISES(i) = sqrt(abs(temp(1).^2 - temp(1).*temp(2) + temp(2).^2 + 3.*temp(3).^2));
    end

    %% P-NORM  
    
    pnorm = (sum(MISES.^p))^(1/p);   
     % dpn_dvms = (sum(MISES.^p)).^((1/p) - 1) * .p .* MISES.^(p - 1);  % !!!!!!!!!!!!!!SENSIBILITA 2D!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     % con questa formula però dpn_dvms non è scalare... genera problemi nel calcolo di T1
   
 %% SENSIBILITA 
  
    dpn_dvms=(sum(MISES.^p))^(1/p-1);
    DvmDs = zeros(nele, 3); % Derivata di von Mises rispetto allo stress

    for i = 1:nele
        % Derivate dello stress di von Mises in base alla componenente
        % dello stress
        DvmDs(i, 1) = 1 / (2 * MISES(i)) * (2 * S(i, 1) - S(i, 2));
        DvmDs(i, 2) = 1 / (2 * MISES(i)) * (2 * S(i, 2) - S(i, 1));
        DvmDs(i, 3) = 3 / MISES(i) * S(i, 3);
        % LA MATEMATICA DIETRO è CORRETTA
    end

    %% T1

    %T1
    beta = zeros(nele, 1);

    for i = 1:nele
        % u = U(edofMat(i, :), :);  % Prende gli spostamenti nodali per l'elemento i
        u=reshape(U(edofMat(i,:),:)',[],1);
        beta(i) = q * (x(i)^(q-1)) * MISES(i)^(p-1) * DvmDs(i, :) * D * B * u;
    end


    % T1 = sum(beta);
    T1 = dpn_dvms * beta;

    %% T2

    index_matrix = edofMat';
    gama = zeros(ndof, 1);


    for i = 1:nele
        index = index_matrix(:, i);
        gama(index) = gama(index) + x(i)^q * dpn_dvms * (B' * D' * (DvmDs(i, :)' * MISES(i).^(p-1)));                                         
    end

    lamda = zeros(ndof, 1);
    lamda(freedofs, :) = K(freedofs, freedofs) \ gama(freedofs, :);

    T2 = zeros(nele, 1);
    for i = 1:nele
        index = index_matrix(:, i);
        T2(i) = -lamda(index)' * pl * x(i)^(pl-1) * KE * U(index);
    end

    %% SENSIBILITA DELLA P-NORM

    pnorm_sen = T1 + T2;


    %% CONTROLLO PER DEBUG
    
    % fprintf('T1')
    % disp(T1(1:10));
    % fprintf('T2')
    % disp(T2(1:10));
    % fprintf('beta')
    % disp(beta(1:10));  % Controlla i primi valori di beta
    % fprintf('DvmDs')
    % disp(DvmDs(1:10, :));  % Verifica i valori delle prime righe di DvmDs
    % fprintf('MISES')
    % disp(MISES(1:10));     % Verifica i valori dei primi elementi di MISES
end