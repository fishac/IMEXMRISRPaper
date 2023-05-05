function lienard_exact()
    clearvars
    
    dydt = @(t,y) dydt_func(t,y);
    tspan = 0:(25/10):25;
    y0 = [1.45; 0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Lienard/Lienard_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Lienard/Lienard_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:25;
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Lienard/Lienard_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Lienard/Lienard_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,25,11);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Lienard/Lienard_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Lienard/Lienard_fixed_truesol_11.csv";
    writematrix(y', filename);
end

function dydt = dydt_func(t,y)
    dydt = [0; 0];
    dydt(1) = y(2);
    dydt(2) = -y(1) - 8.53*(y(1)^2-1)*y(2) + 1.2*sin(pi/5 * t);
end