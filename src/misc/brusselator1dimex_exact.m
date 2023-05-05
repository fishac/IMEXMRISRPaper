function brusselator1dimex_exact()  
    tspan = 0:(3/10):3;
    n = 101;
    y0 = get_y0(n);
    yprime = @(t,y) calc_dydt(t,y);
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/Brusselator1DIMEX/Brusselator1DIMEX_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Brusselator1DIMEX/Brusselator1DIMEX_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:3;
    y0 = get_y0(n);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/Brusselator1DIMEX/Brusselator1DIMEX_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Brusselator1DIMEX/Brusselator1DIMEX_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,3,11);
    y0 = get_y0(n);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/Brusselator1DIMEX/Brusselator1DIMEX_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Brusselator1DIMEX/Brusselator1DIMEX_fixed_truesol_11.csv";
    writematrix(y', filename);
end

function y0 = get_y0(n)
    y0 = zeros(3*n,1);
    u_start = 1;
    v_start = n+1;
    w_start = 2*n+1;
    x = linspace(0,1,n);
    for i=0:(n-1)
        y0(u_start+i) = 1.2 + 0.1*sin(pi*x(i+1));
        y0(v_start+i) = 3.1 + 0.1*sin(pi*x(i+1));
        y0(w_start+i) = 3.0 + 0.1*sin(pi*x(i+1));
    end
end

function doft = d(t)
    doft = 0.00006 + 0.00005*cos(pi*t);
end

function roft = r(t)
    roft = 0.6 + 0.5*cos(4*pi*t);
end

function soft = s(t)
    soft = 0.00006 + 0.00005*cos(pi*t);
end

function f = calc_dydt(t,y)
    a = 1;
    b = 3.5;
    epsilon = 0.001;
    n = 101;
    f = zeros(3*n,1);
    u_start = 1;
    v_start = n+1;
    w_start = 2*n+1;
    dx = 1/(n-1);
    for i=1:(n-2)
    	f(u_start+i) = d(t)*(1.0*y(u_start+i-1) -2.0*y(u_start+i) + 1.0*y(u_start+i+1))/(dx*dx) + s(t)*(-1.0*y(u_start+i-1) + 1.0*y(u_start+i+1))/(dx) + r(t)*(a-(y(w_start+i)+1.0)*y(u_start+i) + y(u_start+i)*y(u_start+i)*y(v_start+i));
		f(v_start+i) = d(t)*(1.0*y(v_start+i-1) -2.0*y(v_start+i) + 1.0*y(v_start+i+1))/(dx*dx) + s(t)*(-1.0*y(v_start+i-1) + 1.0*y(v_start+i+1))/(dx) + r(t)*(y(w_start+i)*y(u_start+i) - y(u_start+i)*y(u_start+i)*y(v_start+i));
		f(w_start+i) = d(t)*(1.0*y(w_start+i-1) -2.0*y(w_start+i) + 1.0*y(w_start+i+1))/(dx*dx) + s(t)*(-1.0*y(w_start+i-1) + 1.0*y(w_start+i+1))/(dx) + r(t)*((b - y(w_start+i))/epsilon - y(w_start+i)*y(u_start+i));
    end
end