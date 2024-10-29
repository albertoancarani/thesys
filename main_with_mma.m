clc
clear
close('all');
warning('off', 'all')
% dbclear all
%% GENERAZIONE GRIGLIA

% Parametri della griglia 
nelx = 200;  % Numero di elementi lungo l'asse x
nely = 100;   % Numero di elementi lungo l'asse y

% Variabile di design (inizialmente con densità uniforme 0.3)
x = 0.3 * ones(nely, nelx);

%% CREAZIONE FILTRO 
% rmin =  0.04*nelx;
rmin =  4;
[Hs, H] = prepare_filter_2D(rmin, nelx, nely);

%% INIZIALIZZAZIONE VARIABILI

% Parametri dell'ottimizzazione
m = 1;                                                                      %YG: # of constraint
epsimin = 1e-7;
n = length(x(:));
xval = x(:);
xold1 = xval;
xold2 = xval;
xlb = 1e-3 * ones(n, 1);  % Limite inferiore del design
xub = ones(n, 1);         % Limite superiore del design
xmin = xlb;
xmax = xub;
low = xlb;
upp = xub;
% c = [1e4]';
% d = [0]';
% a0 = 0;
% a = [0]';
c = 1000*ones(m,1);
d = ones(m,1);
a0 = 0;
a = zeros(m, 1);
raa0 = 0.0001;
raa = 0.0001;
raa0eps = 1e-7;
raaeps = 1e-7;
outeriter = 0;
maxoutit = 200;  % Numero massimo di iterazioni
kkttol = 0;
x_his = zeros(nelx * nely, maxoutit);  % Storico delle variabili di design

% Prima iterazione
if outeriter < 0.5
    [f0val, df0dx, fval, dfdx] = stress_minimize_2D(xval, Hs, H, nelx, nely);  % Versione 2D
    innerit = 0;
    outvector1 = [outeriter, innerit, xval'];
    outvector2 = [f0val, fval'];
end

% Inizializzazione del ciclo di ottimizzazione
kktnorm = kkttol + 1;
outit = 0;
move = 0.2;                                                                  %YG: move-limit
%% Ciclo principale di ottimizzazione
while outit < maxoutit
   
    % Debug all'inizio del ciclo per verificare le variabili chiave
    % disp(['Iterazione: ', num2str(outit)]);
    % disp(['Valore di f0val: ', num2str(f0val)]);
    % disp(['Norma KKT: ', num2str(kktnorm)]);
    % disp(['Valore medio di xval: ', num2str(mean(xval))]);

    outit = outit + 1;
    outeriter = outeriter + 1;
    % 
    % % Impose move limits by modifying lower and upper bounds passed to MMA
    % mlb = max(low, xval - move);
    % mub = min(upp, xval + move);
    % 
    % %%%% Solve MMA subproblem for current design x
    % [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
    %     mmasub(1,n,outeriter,xval,mlb,mub,xold1, ...
    %            xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);


    % Aggiornamento dei limiti inferiori, superiori e parametri raa
    [low, upp, raa0, raa] = ...
        asymp(outeriter, n, xval, xold1, xold2, xmin, xmax, low, upp, ...
              raa0, raa, raa0eps, raaeps, df0dx, dfdx);
    % 
    % % Sottoproblema MMA risolto con gcmmasub (aggiornamento del design)
    [xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, f0app, fapp] = ...
        gcmmasub(m, n, outeriter, epsimin, xval, xmin, xmax, low, upp, ...
                 raa0, raa, f0val, df0dx, fval, dfdx, a0, a, c, d);

    % Aggiornamento delle variabili di design e calcolo dei nuovi valori
    xold2 = xold1;
    xold1 = xval;
    xval = xmma;

    % Calcolo dello stress e delle sensibilità per il nuovo design
    [f0val, df0dx, fval, dfdx] = stress_minimize_2D(xval, Hs, H, nelx, nely);

    % Stampa dei risultati dell'iterazione corrente
    fprintf(' It.:%5i      P-norm Stress.:%11.4f   Vol.:%7.3f \n', ...
            outit, f0val, mean(xval(:)));

    % Calcolo del vettore residuo delle condizioni KKT
    [residu, kktnorm, residumax] = ...
        kktcheck(m, n, xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, ...
                 xmin, xmax, df0dx, fval, dfdx, a0, a, c, d);

    % Aggiornamento dello storico delle variabili di design
    outvector1 = [outeriter, innerit, xval'];
    outvector2 = [f0val, fval'];
    x_his(:, outit) = xmma;
end



