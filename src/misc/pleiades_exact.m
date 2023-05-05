function pleiades_exact()
    clearvars
    
    global g softening_length n d masses;
    g = 1;
    softening_length = 0;
    n = 7;
    d = 2;
    masses = 1:7;
    
    dydt = @(t,y) dydt_func(t,y);
    tspan = 0:(3/10):3;
    positions = [3; 3; 3; -3; -1; 2; -3; 0; 2; 0; -2; -4; 2; 4];
    velocities = [0; 0; 0; 0; 0; 0; 0; -1.25; 0; 1; 1.75; 0; -1.5; 0];
    y0 = [positions; velocities];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Pleiades/Pleiades_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Pleiades/Pleiades_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:3;
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Pleiades/Pleiades_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Pleiades/Pleiades_fixed_truesol_11.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,3,11);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/Pleiades/Pleiades_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Pleiades/Pleiades_fixed_truesol_11.csv";
    writematrix(y', filename);
end

function dydt = dydt_func(t,y)
    global g n d softening_length masses;
    posLength = n * d;

    dydt = [y(posLength + 1:end); zeros(posLength, 1)];

    for i = 1:n
        iStartIdx = d * (i - 1) + 1;
        iEndIdx = d * i;

        posI = y(iStartIdx:iEndIdx);

        for j = (i + 1):n
            jStartIdx = d * (j - 1) + 1;
            jEndIdx = d * j;
            deltaPos = y(jStartIdx:jEndIdx) - posI;

            deltaAccel = g * deltaPos / (sum(deltaPos.^2) + softening_length^2)^(1.5);

            dydt(posLength + (iStartIdx:iEndIdx)) = dydt(posLength + (iStartIdx:iEndIdx)) + masses(j) * deltaAccel;
            dydt(posLength + (jStartIdx:jEndIdx)) = dydt(posLength + (jStartIdx:jEndIdx)) - masses(i) * deltaAccel;
        end
    end
end