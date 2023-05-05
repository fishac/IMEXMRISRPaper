function fourbody3d_exact()
    clearvars
    
    global g softening_length n d masses;
    g = 1;
    softening_length = 0;
    n = 4;
    d = 3;
    masses = [8 10 12 14];
    
    dydt = @(t,y) dydt_func(t,y);
    
    tspan = 0:(15/10):15;
    positions = [ ...
                0; 0; 0; ...
                4; 3; 1; ...
                3; -4; -2; ...
                -3; 4; 5];
    velocities = zeros(12, 1);
    y0 = [positions; velocities];
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/FourBody3d/FourBody3d_t.csv";
    writematrix(t, filename);
    filename = "./../resources/FourBody3d/FourBody3d_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:15;
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/FourBody3d/FourBody3d_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/FourBody3d/FourBody3d_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,15,11);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/FourBody3d/FourBody3d_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/FourBody3d/FourBody3d_fixed_truesol_11.csv";
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