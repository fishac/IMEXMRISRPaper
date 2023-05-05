function grayscott_exact()
    clearvars
    
    global n dspace dspace2 n_interior f_param k_param;
    n=29;
    dspace=1/n;
    dspace2=dspace*dspace;
    n_interior=n-1;
    f_param=0.018;
    k_param=0.052;
    
    dydt = @(t,y) dydt_func(t,y);
    H = 2^-5;
    tspan = 0:H:2;
    y0 = IC();
    
    opts = odeset('RelTol',2.5e-14,'AbsTol',1e-14);
    [t,y] = ode15s(dydt,tspan,y0,opts);
    
    filename = "./../resources/GrayScott/GrayScott_t.csv";
    writematrix(t, filename);
    filename = "./../resources/GrayScott/GrayScott_truesol.csv";
    writematrix(y', filename);
end

function y0 = IC()
    global n_interior dspace;
    y0 = zeros(2*(n_interior+2)^2,1);
    for i = 1:n_interior
        for j=1:n_interior
            if i*dspace >= 0.25 && i*dspace <= 0.75 && j*dspace >= 0.25 && j*dspace <= 0.75 
                y0(u_index_map(i+1,j+1)) = 0.5;
                y0(v_index_map(i+1,j+1)) = 1.0;
            end
        end
    end
end

function u_idx = u_index_map(i, j) 
    global n_interior;
    u_idx = (i-1)*(n_interior+2)+j;
end

function v_idx = v_index_map(i, j) 
    global n_interior;
    v_idx = (n_interior+2)^2+(i-1)*(n_interior+2)+j;
end

function ep_u = epsilon_u(i, j, u) 
    global dspace;
    ep_u = 0.0625*exp(-u/100.0)*sin(pi*i*dspace)*sin(pi*j*dspace);
end

function ep_v = epsilon_v(i, j, v) 
    global dspace;
    ep_v = 0.0312*exp(-v/100.0)*sin(pi*i*dspace)*sin(pi*j*dspace);
end

function dydt = dydt_func(t,y)
    global n_interior dspace2 f_param k_param;
    dydt = zeros(2*(n_interior+2)^2,1);
    for i=1:n_interior 
        for j=1:n_interior 
            u_idx = u_index_map(i+1,j+1);
            u_ip1_idx = u_index_map(i+1+1,j+1);
            u_im1_idx = u_index_map(i-1+1,j+1);
            u_jp1_idx = u_index_map(i+1,j+1+1);
            u_jm1_idx = u_index_map(i+1,j-1+1);

            v_idx = v_index_map(i+1,j+1);
            v_ip1_idx = v_index_map(i+1+1,j+1);
            v_im1_idx = v_index_map(i-1+1,j+1);
            v_jp1_idx = v_index_map(i+1,j+1+1);
            v_jm1_idx = v_index_map(i+1,j-1+1);

            u = y(u_idx);
            u_ip1 = y(u_ip1_idx);
            u_im1 = y(u_im1_idx);
            u_jp1 = y(u_jp1_idx);
            u_jm1 = y(u_jm1_idx);

            v = y(v_idx);
            v_ip1 = y(v_ip1_idx);
            v_im1 = y(v_im1_idx);
            v_jp1 = y(v_jp1_idx);
            v_jm1 = y(v_jm1_idx);
            
            offset_u_factor = (epsilon_u(i+1,j,u)-epsilon_u(i-1,j,u)+epsilon_u(i,j+1,u)-epsilon_u(i,j-1,u))/(4.0*dspace2) + epsilon_u(i,j,u)/dspace2;
			offset_v_factor = (epsilon_v(i+1,j,v)-epsilon_v(i-1,j,v)+epsilon_v(i,j+1,v)-epsilon_v(i,j-1,v))/(4.0*dspace2) + epsilon_v(i,j,v)/dspace2;

            dydt(u_idx) = u*(-4.0*epsilon_u(i,j,u)/dspace2 - f_param)+f_param ...
            + u_ip1*offset_u_factor ...
            + -u_im1*offset_u_factor ...
            + u_jp1*offset_u_factor ...
            + -u_jm1*offset_u_factor ...
            + -u*v*v;

            dydt(v_idx) = v*(-4.0*epsilon_v(i,j,v)/dspace2 - (f_param+k_param)) ...
            + v_ip1*offset_v_factor ...
            + -v_im1*offset_v_factor ...
            + v_jp1*offset_v_factor ...
            + -v_jm1*offset_v_factor ...
            + u*v*v;
        end
    end
end