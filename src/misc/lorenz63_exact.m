function lorenz63_exact()
    sigma = 10.0;
    beta = 8.0/3.0;
    rho = 28.0;
    dydt = @(t,y) [sigma*(y(2)-y(1)); ...
        y(1)*(rho-y(3))-y(2); ...
        y(1)*y(2)-beta*y(3)];
    
    tspan = 0:(2/10):2;
    y0 = [1.0; 1.0; 1.0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Lorenz63/Lorenz63_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Lorenz63/Lorenz63_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:2;
    y0 = [1.0; 1.0; 1.0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Lorenz63/Lorenz63_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Lorenz63/Lorenz63_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,2,11);
    y0 = [1.0; 1.0; 1.0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Lorenz63/Lorenz63_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Lorenz63/Lorenz63_fixed_truesol_11.csv";
    writematrix(y', filename);
end