function brusselator_exact()
    a = 1;
    b = 3.5;
    epsilon = 0.01;
    dydt = @(t,y) [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2); ...
        y(3)*y(1) - y(1)*y(1)*y(2); ...
        (b-y(3))/epsilon - y(1)*y(3)];
    
    tspan = 0:(2/10):2;
    y0 = [1.2; 3.1; 3];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/brusselator/brusselator_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Brusselator/brusselator_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:2;
    y0 = [1.2; 3.1; 3];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/brusselator/brusselator_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Brusselator/brusselator_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,2,11);
    y0 = [1.2; 3.1; 3];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/brusselator/brusselator_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Brusselator/brusselator_fixed_truesol_11.csv";
    writematrix(y', filename);
end