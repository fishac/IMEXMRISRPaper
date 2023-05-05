function brusselator1d_exact()  
    tspan = 0:(2/10):2;
    n = 100;
    y0 = get_y0(n);
    yprime = @(t,y) calc_dydt(t,y);
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/Bicoupling1D/bicoupling1d_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Bicoupling1D/bicoupling1d_truesol.csv";
    writematrix(y', filename);
    
    H = 2^-5;
    tspan = 0:H:2;
    y0 = zeros(3*n,1);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/Bicoupling1D/bicoupling1d_fixed_t.csv";
    writematrix(t, filename);
    filename = "./../resources/Bicoupling1D/bicoupling1d_fixed_truesol.csv";
    writematrix(y', filename);
    
    tspan = linspace(0,2,11);
    y0 = zeros(3*n,1);
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(yprime,tspan,y0,opts);
    
    filename = "./../resources/Bicoupling1D/bicoupling1d_fixed_t_11.csv";
    writematrix(t, filename);
    filename = "./../resources/Bicoupling1D/bicoupling1d_fixed_truesol_11.csv";
    writematrix(y', filename);
end

function y0 = get_y0(n)
    y0 = zeros(3*n,1);
    u_start = 1;
    v_start = n+1;
    w_start = 2*n+1;
    x = linspace(0,1,n);
	a = 1.0;
	b = 20.0;
	w = 100.0;
	l = 5.0;
	p = 0.001;
    for i=0:(n-1)
        y0(u_start+i) = 1.0+a + 0.1*sin(pi*x(i+1));
        y0(v_start+i) = b + 0.1*sin(pi*x(i+1));
        y0(w_start+i) = a*l+b*w + 0.1*sin(pi*x(i+1));
    end
end

function doft = d(t)
    doft = 0.006 + 0.005*cos(pi*t);
end

function roft = r(t)
    roft = 0.6 + 0.5*cos(4*pi*t);
end

function f = calc_dydt(t,y)
	a = 1.0;
	b = 20.0;
	w = 100.0;
	l = 5.0;
	p = 0.001;
    n = 100;
    f = zeros(3*n,1);
    u_start = 1;
    v_start = n+1;
    w_start = 2*n+1;
    dx = 1/(n-1);
    for i=1:(n-2)
    	f(u_start+i) = d(t)*(1.0*y(u_start+i-1) -2.0*y(u_start+i) + 1.0*y(u_start+i+1))/(dx*dx) + r(t)*(w*y(v_start+i)-y(w_start+i)-p*t);
		f(v_start+i) = d(t)*(1.0*y(v_start+i-1) -2.0*y(v_start+i) + 1.0*y(v_start+i+1))/(dx*dx) + r(t)*(-w*y(u_start+i));
		f(w_start+i) = d(t)*(1.0*y(w_start+i-1) -2.0*y(w_start+i) + 1.0*y(w_start+i+1))/(dx*dx) + r(t)*(-l*y(w_start+i)-l*p*t-p*(y(u_start+i)-a*y(w_start+i)/(a*l + b*w)-a*p*t/(a*l + b*w)).^2-p*(y(v_start+i)-b*y(w_start+i)/(a*l + b*w)-b*p*t/(a*l + b*w)).^2);
    end
end