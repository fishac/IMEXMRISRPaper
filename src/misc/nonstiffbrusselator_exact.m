function nonstiffbrusselator_exact()
    a = 1.0;
    b = 3.0;
    dydt = @(t,y) [a - y(1)*(b+1.0) + y(1)*y(1)*y(2); ...
        b*y(1) - y(1)*y(1)*y(2)];
    
    tspan = 0:(2/10):2;
    y0 = [1.0; 1.0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/nonstiffbrusselator/nonstiffbrusselator_t.csv";
    writematrix(t, filename);
    filename = "./../resources/nonstiffbrusselator/nonstiffbrusselator_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:2;
    y0 = [1.0; 1.0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/nonstiffbrusselator/nonstiffbrusselator_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/nonstiffbrusselator/nonstiffbrusselator_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,2,11);
    y0 = [1.0; 1.0];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/nonstiffbrusselator/nonstiffbrusselator_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/nonstiffbrusselator/nonstiffbrusselator_fixed_truesol_11.csv";
    writematrix(y', filename);
end