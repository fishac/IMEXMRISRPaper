function brusselatorpde801_exact()  
    tspan = 0:(3/10):3;
    n = 801;
    y0 = get_y0(n);
    yprime = @(t,y) calc_dydt(t,y);
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/BrusselatorPDE801/BrusselatorPDE801_t.csv";
    writematrix(t, filename);
    filename = "./../resources/BrusselatorPDE801/BrusselatorPDE801_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:3;
    y0 = get_y0(n);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/BrusselatorPDE801/BrusselatorPDE801_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/BrusselatorPDE801/BrusselatorPDE801_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,3,11);
    y0 = get_y0(n);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/BrusselatorPDE801/BrusselatorPDE801_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/BrusselatorPDE801/BrusselatorPDE801_fixed_truesol_11.csv";
    writematrix(y', filename);
end

function y0 = get_y0(n)
    y0 = zeros(3*n,1);
    u_start = 1;
    v_start = n+1;
    w_start = 2*n+1;
    x = linspace(0,1,n);
    a = 0.6;
    b = 2;
    epsilon = 0.01;
    for i=0:(n-1)
        y0(u_start+i) = a + 0.1*sin(pi*x(i+1));
        y0(v_start+i) = b/a + 0.1*sin(pi*x(i+1));
        y0(w_start+i) = b + 0.1*sin(pi*x(i+1));
    end
end


function f = calc_dydt(t,y)
    a = 0.6;
    b = 2.0;
    epsilon = 0.01;
    alpha = 0.01;
    rho = 0.001;
    n = 801;
    f = zeros(3*n,1);
    u_start = 1;
    v_start = n+1;
    w_start = 2*n+1;
    dx = 1/(n-1);
    for i=1:(n-2)
    	f(u_start+i) = alpha*(1.0*y(u_start+i-1) -2.0*y(u_start+i) + 1.0*y(u_start+i+1))/(dx*dx) - rho*(-y(u_start+i-1)+y(u_start+i+1))/(2.0*dx) + (a-(y(w_start+i)+1.0)*y(u_start+i) + y(u_start+i)*y(u_start+i)*y(v_start+i));
		f(v_start+i) = alpha*(1.0*y(v_start+i-1) -2.0*y(v_start+i) + 1.0*y(v_start+i+1))/(dx*dx) - rho*(-y(v_start+i-1)+y(v_start+i+1))/(2.0*dx) + (y(w_start+i)*y(u_start+i) - y(u_start+i)*y(u_start+i)*y(v_start+i));
		f(w_start+i) = alpha*(1.0*y(w_start+i-1) -2.0*y(w_start+i) + 1.0*y(w_start+i+1))/(dx*dx) - rho*(-y(w_start+i-1)+y(w_start+i+1))/(2.0*dx) + ((b - y(w_start+i))/epsilon - y(w_start+i)*y(u_start+i));
    end
end