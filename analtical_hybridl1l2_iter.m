function [anl] = analtical_hybridl1l2_iter(dx,alpha,mu,H0,beta_sl,solver)
    % Calculate analytical solution for a given value of dx
    
    % Define constants
    rhog = 910*9.8;
    F2   = H0/(3*mu);
    B    = beta_sl; 
    D    = 4*mu*H0/dx^2;
    Q    = B/D;
    r    = 1 + Q/2 - sqrt((1+Q/2)^2-1);
    a0   = 1/(2+Q-2*r);
    
    eta  = beta_sl*H0/mu;
    
    % Define a range of dt values to try
    dts = 10.^(linspace(log10(1e-4),log10(50),1000));
    
    % Calculate the limit
    if (solver == "hybrid")
        val = abs(1 - dts*H0*a0*(rhog/mu)*((1-r)/(1+r)) - (4*rhog*H0^3)/(3*mu)*dts/dx^2);
    else  % solver=="l1l2"
        val = abs( 1 - dts*H0*a0*(rhog/mu)*((1-r)/(1+r))*(1+eta/3) );
    end
    
    % Determine maximum viable dt value
    ii = find(val < 1);
    if (length(ii) > 0)
        dt = max(dts(ii));
    else
        dt = NaN;
    end
    
    % Store analytical maximum stable timestep 
    anl = dt; 
    
end