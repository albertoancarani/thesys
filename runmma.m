function history = runmma(x0,obj,nonlcon)

global OPT GEOM FE LAM

%%

% Initialize history object
history.x = [];
history.fval = [];
history.fconsval = [];

% Initialize lower and upper bounds vectors
if OPT.options.dv_scaling
    lb_point = zeros(FE.dim,1);
    ub_point = ones(FE.dim,1);
    lb_radius = 0;
    % Consider case when max_bar_radius and min_bar_radius are
    % the same (when bars are of fixed radius)
    if GEOM.max_bar_radius - GEOM.min_bar_radius < 1e-12
        ub_radius = 0;
    else
        ub_radius = 1;
    end
else
    lb_point = FE.coord_min;
    ub_point = FE.coord_max;
    lb_radius = GEOM.min_bar_radius;
    ub_radius = GEOM.max_bar_radius;
end
lb_size = zeros(FE.n_mat,1);
ub_size = ones (FE.n_mat,1);

lb_bar = [lb_point;lb_point;lb_radius;lb_size];
ub_bar = [ub_point;ub_point;ub_radius;ub_size];

lb = zeros(size(OPT.dv)); 
ub = zeros(size(OPT.dv)); 
lb(OPT.bar_dv) = repmat(lb_bar,1,GEOM.n_bar);
ub(OPT.bar_dv) = repmat(ub_bar,1,GEOM.n_bar);


%
ncons = OPT.functions.n_func - 1;  % Number of optimization constraints
ndv = OPT.n_dv; % Number of design variables

% Initialize vectors that store current and previous two design iterates
x = x0;
xold1 = x0; 
xold2 = x0;

OPT.xold1 = x0;
OPT.xold2 = x0;
OPT.x = x0;
% Initialize move limits 
ml_step = OPT.options.move_limit * abs(ub - lb);  % Compute move limits once

% Initialize lower and upper asymptotes
low = lb;
upp = ub;
OPT.low = low;
OPT.upp = upp;

% These are the MMA constants (Svanberg, 1998 DACAMM Course)
c = 1000*ones(ncons,1);
d = ones(ncons,1);
a0 = 1;
a = zeros(ncons, 1);

% Evaluate the initial design and print values to screen 
iter = 0;
[f0val , df0dx] = obj(x);
[fval, ~, dfdx, ~] = nonlcon(x);
dfdx = dfdx';
fprintf('It. %i, Obj= %-12.5e, ConsViol = %-12.5e\n', ...
    iter, f0val, max(max(fval, zeros(ncons,1))));
OPT.iter = iter;
OPT.f0old = f0val;
OPT.f0val = f0val;
OPT.df0dx = df0dx;
OPT.fval = fval;
OPT.dfdx = dfdx;


%%%
% Save initial design to history
history.fval = [history.fval; f0val];
history.fconsval = [history.fconsval; fval];
history.x = [history.x x(:)];

%%%
% Plot initial design 
plotfun(iter);
          
%%%% Initialize stopping values
kktnorm = OPT.options.kkt_tol*10;
dv_step_change = 10*OPT.options.step_tol;
if isfield(OPT.options,"objective_tol")
    obj_step_change = 10*abs(OPT.options.objective_tol);
else
    OPT.options.objective_tol = -1; % value is absolute change
    obj_step_change = 1;
end
%%%% Continuation method
bufferSize   = 10;
bufferIndex = 1;
buffer = zeros(bufferSize,1);
%
% ******* MAIN MMA LOOP STARTS *******
%
LAM.OPT(:) = {OPT};
LAM.GEOM(:) = {GEOM};   
LAM.FE(:) = {FE};
L = 1:LAM.num_ply;

while kktnorm > OPT.options.kkt_tol && iter < OPT.options.max_iter && ...
        dv_step_change > OPT.options.step_tol && ...
        obj_step_change > OPT.options.objective_tol

    iter = iter+1;
    LAM.iter = iter;

    for l = L
        LAM.count = l;
        GEOM =  LAM.GEOM{l}; FE = LAM.FE{l}; OPT = LAM.OPT{l};
% 
        x = OPT.x; xold1 = OPT.xold1; xold2 = OPT.xold2;
        f0val = OPT.f0val; fval = OPT.fval;
        df0dx = OPT.df0dx; dfdx = OPT.dfdx;
        low = OPT.low; upp = OPT.upp;
        % Impose move limits by modifying lower and upper bounds passed to MMA
        mlb = max(lb, x - ml_step);
        mub = min(ub, x + ml_step);
    
        %%%% Solve MMA subproblem for current design x
        [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
        mmasub(ncons,ndv,iter,x,mlb,mub,xold1, ...
               xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
    
        %%%% Updated design vectors of previous and current iterations
        xold2 = xold1;
        xold1 = x;
        x  = xmma;

        OPT.xold1 = xold1; OPT.xold2 = xold2; OPT.x = x;
        OPT.low = low; OPT.upp = upp;
        
        update_and_project(x);
        LAM.GEOM{l} = GEOM; LAM.FE{l} = FE; LAM.OPT{l} = OPT;
    end

    for l = L
        % Update function values and gradients
        LAM.count = l;
        GEOM =  LAM.GEOM{l}; FE = LAM.FE{l}; OPT = LAM.OPT{l};
        
%         laminate_analysis_v2();
        evaluate_relevant_functions();
        x = OPT.x; xold1 = OPT.xold1; f0val = OPT.f0val; 
        f0old = f0val;
        [f0val , df0dx] = obj(x);
        [fval, ~, dfdx, ~] = nonlcon(x);
        dfdx = dfdx';

        OPT.f0old = f0old; OPT.f0val = f0val; OPT.df0dx = df0dx;
        OPT.fval = fval; OPT.dfdx = dfdx;
        
        % Compute change in design variables
        % Check only after first iteration
        if iter > 1
            dv_step_change = norm(x - xold1);
        
            % if all(ml_step < 0.05)
            %     ml_step  = ml_step + 2.0e-04;
            % end
            
            % if all(bufferIndex > bufferSize) && all(ml_step > 0.005)
            %     if all(buffer < dv_step_change) 
            %         ml_step = ml_step*exp(-0.2);
            %     end
            %     buffer = circshift(buffer,[1,0]);
            %     buffer(1) = dv_step_change;
            % else
            %     buffer = circshift(buffer,[1,0]);
            %     buffer(1) = dv_step_change;
            %     bufferIndex = bufferIndex + 1;
            % end

            if dv_step_change < OPT.options.step_tol
                L = setxor(l,L);
                fprintf('Design step convergence tolerance satisfied.\n');
            end
            obj_step_change = norm(f0val - f0old);
            if obj_step_change < OPT.options.objective_tol
                L = setxor(l,L);
                fprintf('Ojective function convergence tolerance satisfied.\n');
            end
        end
        if iter == OPT.options.max_iter
            fprintf('Reached maximum number of iterations.\n');
        end    
        
        % Compute norm of KKT residual vector
        [residu,kktnorm,residumax] = ...
        kktcheck(ncons,ndv,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
               lb,ub,df0dx,fval,dfdx,a0,a,c,d);

        % Produce output to screen
        fprintf('It. %i, Obj= %-12.5e, ConsViol = %-12.5e, KKT-norm = %-12.5e, DV norm change = %-12.5e\n', ...
            iter, f0val, max(max(fval, zeros(ncons,1))), kktnorm, dv_step_change);
        
        % Save design to .mat file
        [folder, baseFileName, ~] = fileparts(GEOM.initial_design.path);
        mat_filename = fullfile(folder, strcat(baseFileName, '.mat'));
    %     save(mat_filename, 'GEOM');
        
        
        % Update history
        history.fval = [history.fval; f0val];
        history.fconsval = [history.fconsval; fval];
        history.x = [history.x x(:)];
        
        % Plot current design
    
        plotfun(iter);
        LAM.GEOM{l} = GEOM; LAM.FE{l} = FE; LAM.OPT{l} = OPT;
    end
   
%
end

% Write vtk for final iteration if requested
if strcmp(OPT.options.write_to_vtk, 'all') || ...
        strcmp(OPT.options.write_to_vtk, 'last')
    writevtk(OPT.options.vtk_output_path, 'dens', iter);
end 

% ============================================


    function plotfun(iter)
        % Note that this function has a slightly different format than its
        % equivalent for fmincon.
        % here.
        
        if OPT.options.plot == true   
            if mod(iter,1) == 0
                
                if iter>0
                    if ~(isfield(OPT.options.plotting.design,'axes'))
                       figure(OPT.options.plotting.design.figure);       
                       OPT.options.plotting.design.axes = cell(1,LAM.num_ply);
                       for il = L
                           OPT.options.plotting.design.axes{il} = subplot(1,LAM.num_ply,il);
%                            OPT.options.plotting.design.axes{2} = subplot(1,LAM.num_ply,2);
                       end
                    end
                    OPT.options.plotting.design.axis = OPT.options.plotting.design.axes{l};
                    
                    plot_design(OPT.options.plotting.design.axis)
                    title(OPT.options.plotting.design.axis,sprintf('design, iteration = %i',iter))
                else
                    plot_design(OPT.options.plotting.design.axis)
                    title(OPT.options.plotting.design.axis,sprintf('design, iteration = %i',iter))
                end
                
                for m = 1:FE.n_mat
                    axis(OPT.options.plotting.design.axis,'equal')
                    xlim(OPT.options.plotting.design.axis,[FE.coord_min(1), FE.coord_max(1)])
                    ylim(OPT.options.plotting.design.axis,[FE.coord_min(2), FE.coord_max(2)])
                    if FE.dim == 2
                        view(OPT.options.plotting.design.axis,2)
                    else
                        zlim(OPT.options.plotting.design.axis,[FE.coord_min(3), FE.coord_max(3)])
                        view(OPT.options.plotting.design.axis,[50,22])
                    end
                end
         
                if iter>0
                    plot_density(l)
                else
                    plot_density()               
                end
                
                for m = 1:FE.n_mat
                    axis(OPT.options.plotting.density.axes{m},'equal')
                    xlim(OPT.options.plotting.density.axes{m},[FE.coord_min(1), FE.coord_max(1)])
                    ylim(OPT.options.plotting.density.axes{m},[FE.coord_min(2), FE.coord_max(2)])
                    if FE.dim == 2
                        view(OPT.options.plotting.density.axes{m},2)
                    else
                        zlim(OPT.options.plotting.density.axes{m},[FE.coord_min(3), FE.coord_max(3)])
                        view(OPT.options.plotting.density.axes{m},[50,22])
                    end
                end
               drawnow; 
            end
        end
    end
end