function [f0val, df0dx, fval, dfdx] = stress_minimize_2D(x, Hs, H, nelx, nely)
% Parametri
pl = 3;       % Parametro di penalizzazione
q = 0.5;      % Parametro per sensibilità
p = 10;       % Parametro per norma p

% Applica il filtro (nel caso 2D non c'è più la terza dimensione)
x(:) = (H * x(:)) ./ Hs;

% Calcola la norma p dello stress di Von Mises e la sua sensibilità
[pnorm, pnorm_sen, MISES] = Stress_2D_Sensitivity_Comp(x, nelx, nely, pl, q, p);

% Visualizzazione della distribuzione del materiale
figure(1)
x_plot = reshape(x, nely, nelx);
contourf(flipud(x_plot), [0.5, 0.5]);
colormap([0, 0, 0]);
set(gcf, 'color', 'w');
axis equal;
axis off;
title('Material Layout');
drawnow;



% Visualizzazione dello stress di von Mises
figure(2)
imagesc(reshape(MISES .* (0.5 * sign(x - 0.5) + 0.5), nely, nelx));
axis equal;
axis off;
colormap('jet');
title('Von-Mises Stress');
colorbar;
drawnow;

% Calcolo delle sensibilità e del gradiente rispetto al design
dv = ones(nely, nelx) / (nelx * nely);  % Derivata del volume
sen(:) = H * (pnorm_sen(:) ./ Hs);      % Sensibilità normalizzata
dv(:) = H * (dv(:) ./ Hs);              % Sensibilità del volume normalizzata

% Vincolo sul volume e sensibilità
fval = [mean(x(:)) - 0.3];  % Vincolo sul volume (0.3 è il target di volume)
dfdx = [dv(:)'];            % Sensibilità del vincolo di volume

% Sensibilità dell'obiettivo
df0dx = sen';

% Valore dell'obiettivo (norma p dello stress)
f0val = pnorm;

%%CONTROLLO PER DEBUG
    % disp('df0dx (first few elements):');
    % disp(df0dx(1:10));  % Mostra i primi 10 elementi di lamda
    % disp('dfdx (first few elements):');
    % disp(dfdx(1:10));  % Mostra i primi 10 elementi di lamda
end

