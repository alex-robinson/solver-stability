clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User options

% Choose experiment (weak, strong): 
experiment = 'weak';

% Choose solver (diva, hybrid, l1l2):
solver     = 'diva';

% Calculate and save the sample (true), or simply load a saved version (false)?
run_sample = true; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional parameters
alpha    = 1e-3;
epsh     = 0.1;
nsteps   = 100;

% Sampling parameters
lx_bound = [5 45e3];
lt_bound = [1e-3 50]; 
n_x      = 60;
n_t      = 20;

% Define constants depending on experiment
if (experiment=="weak") 
    L       = 1e5;
    mu      = 1e5;
    H0      = 1000;
    beta_sl = 1e3;
else % experiment=='strong'
    L       = 4e5;
    mu      = 4e5;
    H0      = 500;
    beta_sl = 30;
end

% Calculate additional constants
rhog = 910*9.8;
eta  = beta_sl*H0/mu; 

% Define a vector of dx values for analytical solutions 
dx_anl = 10.^(linspace(log10(lx_bound(1)),log10(lx_bound(2)),100));
    
if (solver=="diva")
    
    % Define specific solver function
    stability_solver = str2func('stability_test1D_diva');
    
    % Calculate analytical solution
    anl = analytical_diva(dx_anl,alpha,mu,H0,beta_sl);
    
elseif (solver=="l1l2")
    
    % Define specific solver function
    stability_solver = str2func('stability_test1D_l1l2');
    
    % Calculate analytical solution
    anl = analytical_l1l2(dx_anl,alpha,mu,H0,beta_sl);
    
else % Hybrid
    
    % Define specific solver function
    stability_solver = str2func('stability_test1D_hybrid');
    
    % Calculate analytical solution
    anl = analytical_hybrid(dx_anl,alpha,mu,H0,beta_sl);
    
end

if (run_sample)
    % Sample dx/dt space and save results to a .mat file 
    
    % Calculate dx values to test
    dx      = 10.^(linspace(log10(lx_bound(1)),log10(lx_bound(2)),n_x));

    % Generate a table to hold the results
    clear res;
    res.dx  = zeros(n_x*n_t,1);
    res.dt  = zeros(n_x*n_t,1);
    res.fac = zeros(n_x*n_t,1);

    ktot = 0; 

    for nx = 1:length(dx)
    %for nx = 1:1

        % Generate 5 dt values
        dt      = zeros(n_t);
        dt(1:5) = 10.^(linspace(log10(lt_bound(1)),log10(lt_bound(2)),5));

        factor  = zeros(n_t);

        for nt = 1:5

            [fac u0 Lm] = stability_solver(dx(nx),dt(nt),nsteps,...
                                        beta_sl,L,mu,H0,alpha,epsh);

            factor(nt) = fac;

            % Set convergence to high value for now for printing 
            conv = 1e5;

            % Print results
            disp(['dx, dt, fac, conv = ' num2str(dx(nx)) ', ' ...
                    num2str(dt(nt)) ', ' num2str(fac) ', ' num2str(conv)]);

            % Save results to output table
            ktot = ktot + 1;
            res.dx(ktot)  = dx(nx);
            res.dt(ktot)  = dt(nt);
            res.fac(ktot) = factor(nt);

        end

        if (min(res.fac(1:ktot)) < 1 & max(res.fac(1:ktot)) > 1)
            % Initial sampling covers stable and unstable range,
            % so now sample further to refine estimate

            for nt = 6:n_t

                kk = find(dt > 0 & factor < 1);
                [val idx] = max(dt(kk));
                qmax = kk(idx);
                kk = find(dt > dt(qmax));
                [val idx] = min(dt(kk));
                qmax1 = kk(idx);

                if (length(qmax1)==1)
                    dt(nt) = 0.5*(dt(qmax)+dt(qmax1));
                else
                    dt(nt) = dt(qmax);
                end

                [fac u0 Lm] = stability_solver(dx(nx),dt(nt),nsteps,...
                                        beta_sl,L,mu,H0,alpha,epsh);

                factor(nt) = fac;

                conv = abs(dt(nt)-dt(nt-1))/dt(nt-1);

                % Print results
                disp(['dx, dt, fac, conv = ' num2str(dx(nx)) ', ' ...
                        num2str(dt(nt)) ', ' num2str(fac) ', ' num2str(conv)]);

                % Save results to output table
                ktot = ktot + 1;
                res.dx(ktot)  = dx(nx);
                res.dt(ktot)  = dt(nt);
                res.fac(ktot) = factor(nt);

                % Check for convergence
                if ( conv < 1e-3) 
                    break;
                end
            end   

        end

        % Print empty line to start next loop 
        disp('  ')

    end

    % Select non-zero indices 
    kk = find(res.dt > 0);
    res.dx  = res.dx(kk);
    res.dt  = res.dt(kk);
    res.fac = res.fac(kk);

    % Save output table
    save(['sampled_' solver '_' experiment '.mat'], '-struct', 'res')
    
else
    % Load sample from saved .mat file
    
    % Load output table
    load(['sampled_' solver '_' experiment '.mat'])
    clear res;
    res.dx  = dx;
    res.dt  = dt;
    res.fac = fac; 
    
end

% Plot to be sure everything looks ok
plot(res.dx,res.dt,'.')
set(gca,'xscale','Log');
set(gca,'yscale','Log');

% Add analytical solution
hold on;
plot(anl.dx,anl.dt_adv,'--','color','black','linewidth',1)
plot(anl.dx,anl.dt_dyn,'-', 'color','black','linewidth',2)

