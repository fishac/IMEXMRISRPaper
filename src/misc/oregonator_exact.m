function oregonator_exact()
    clearvars
    
    dydt = @(t,y) dydt_func(t,y);

    tspan = 0:(360/100):360;
    y0 = [1; 2; 3];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Oregonator/Oregonator_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Oregonator/Oregonator_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:360;
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Oregonator/Oregonator_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Oregonator/Oregonator_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,360,11);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Oregonator/Oregonator_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Oregonator/Oregonator_fixed_truesol_11.csv";
    writematrix(y', filename);
end

function dydt = dydt_func(t,y)
    dydt = [0; 0; 0];
    dydt(1) = 77.27*(y(2)+y(1)*(1-(8.375e-6)*y(1)-y(2)));
    dydt(2) = (y(3)-(1+y(1))*y(2))/77.27;
    dydt(3) = 0.161*(y(1)-y(3));
end