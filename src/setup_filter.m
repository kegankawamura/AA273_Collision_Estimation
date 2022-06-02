% script that creates the filter and contains the internal functions
function [filter] = setup_filter(robot_params,map,mu_robot,sigma_robot)
    has_prior = false;
    if nargin==4
        has_prior = true;
    end

    % number of particles
    N = 1000;


    dt = robot_params.dt;
    m = robot_params.m;
    I = robot_params.I;

    % made up numbers
    Q_robot = dt*blkdiag(0.5*eye(2), 0.5, 0.2*eye(2), 0.2);
    Q_omega_bounce = 0.3;

    Q_ekf = blkdiag(0.5*eye(2),1,eye(2)*3);
    % R is augmented with the force pseudomeasurements
    R_ekf = 5*blkdiag(0.5*eye(3),1.5*eye(2),1);
    sigma_0_ekf = 10*eye(5);

    pf_params = robot_params;
    pf_params.Map = map;
    pf_params.Q_robot = Q_robot;
    pf_params.Q_omega_bounce = Q_omega_bounce;

    pf_system = struct('params',pf_params,...
        'transition',@trans_model ...
        );

    ekf_system = struct('params',pf_params,...
        'dyn',@dyn_ekf,...
        'meas',@meas_ekf,...
        'A',@dyn_jac_ekf,...
        'C',@meas_jac_ekf,...
        'Qmodel',Q_ekf,...
        'Rmodel',R_ekf ...
        );
    ekf_system = EKFSystem(ekf_system);

    filter = CIOFilter(pf_system,ekf_system,robot_params);
    if has_prior
        initialize_filter_prior(filter,mu_robot,sigma_robot,0.1*ones(5,1),sigma_0_ekf,map,N,robot_params);
    else
        initialize_filter_uniform(filter,0.1*ones(5,1),sigma_0_ekf,map,N,robot_params);
    end
end

% initialize filter with known mean + std of pose
function [particles] = initialize_filter_prior(filter,mu_robot,sigma_robot,mu_ekf,sigma_ekf,map,N,robot_params)

    bounding_box = [min(map.corners);max(map.corners)];
    span = diff(bounding_box);
    P = mvnrnd(mu_robot,sigma_robot,N)';
    %P = [rand(N,2).*span + bounding_box(1,:)]';

    for i=1:N
        while ~map.is_inside(P(1:2,i))
            P(:,i) = mvnrnd(mu_robot,sigma_robot,1)';
        end
    end
    particles = struct('X',P,'W',ones(1,N)*1/N);
    filter.initialize(particles,mu_ekf,sigma_ekf);
end

% initialize filter with uniform particles distributed in the map
function [particles] = initialize_filter_uniform(filter,mu,sigma,map,N,robot_params)
    bounding_box = [min(map.corners);max(map.corners)];
    span = diff(bounding_box);
    P = [rand(N,2).*span + bounding_box(1,:)]';

    for i=1:N
        while ~map.is_inside(P(:,i))
            P(:,i) = [rand(1,2).*span + bounding_box(1,:)]';
        end
    end
    V = [2*(0.5-rand(2,N))] * robot_params.v_max;
    om = [2*(0.5-rand(1,N))] * robot_params.om_max;
    th = [2*(0.5-rand(1,N))] * pi;

    X = [P;th;V;om];
    particles = struct('X',X,'W',ones(1,N)*1/N);

    filter.initialize(particles,mu,sigma);
end


% simplified collision model
function [X_t1s] = trans_model(X_ts,u_t,params)
    dt = params.dt;

    X_t1s = X_ts;
    for i=1:size(X_ts,2)
        X_t1s(:,i) = TruthSim.collision_model(X_ts(:,i),params);
    end
    % pass through point mass dynamics and gaussian noise
    %w = mvnrnd(zeros(size(X_ts(:,1))),params.Q_robot,size(X_ts(1,:),2))';
    %X_t1s = quad_dyn(X_ts,u_t,params) + w;

    %% check for collision
    %for i=1:size(X_ts,2)
    %    X_prev = X_ts(:,i); X_next = X_t1s(:,i);

    %    [hitsWall, wall_normal, collision_loc] = params.Map.hit_wall(X_prev(1:2),X_next(1:2),params.R);

    %    if hitsWall
    %        if X_prev(4) ~= 0
    %            pre_dt = (collision_loc(1) - X_prev(1))/X_prev(4);
    %            post_dt = params.dt - pre_dt;
    %        else
    %            pre_dt = (collision_loc(2) - X_prev(2))/X_prev(5);
    %            post_dt = params.dt - pre_dt;
    %        end
    %        v_perp = wall_normal * (wall_normal'*X_prev(4:5));
    %        v_paral  = X_prev(4:5)-v_perp;

    %        % assumes perfect collision
    %        post_v_perp = -1*v_perp;
    %        %acc = (post_v_perp-v_perp)/dt;
    %        X_next(4:5) = post_v_perp + v_paral;
    %        X_next(1:2) = collision_loc + post_dt*X_next(4:5) + wall_normal*1e-6;

    %        % noise in angular velocity due to bounce
    %        X_next(6) = X_next(6) + randn()*params.Q_omega_bounce;
    %    else
    %        %acc = (X_next(4:5)-X_prev(4:5))/dt;
    %    end
    %    X_t1s(:,i) = X_next;
    %end
end

% point mass linear dyn
function [X_t1] = quad_dyn(X_t,u_t,params)
    dt = params.dt;
    F = u_t(1:2); M = u_t(3);
    X_t1(1:3,:) = X_t(1:3,:) + dt*X_t(4:6,:);
    X_t1(4:5,:) = X_t(4:5,:) + dt*F/params.m;
    X_t1(6,:) = X_t(6,:) + dt*M/params.I;
end

% ekf params will have 
%   params.num - number of ekf filter if in array (0 otherwise)
%   params.Xprev - states from last filter step
%   params.X   - states of each of the particles
%   need both to create a finite difference acceleration
% ekf state is 
% [ b_ax
%   b_ay
%   b_w
%   alpha_a
%   alpha_w ]
% measurement model is
% y = [f_b;om]
%
% f_b = R(th)*a + b_a + v_a1 , a = (V-Vprev)/dt
% om  = om + b_w + v_w1
%

% dynamics are a first order markov OU process
function [X_t1] = dyn_ekf(X_t,u_t,params)
    X_t1 = X_t;
    X_t1(1:2,:) = X_t(4,:).*X_t(1:2,:);
    X_t1(3,:) = X_t(5,:).*X_t(3,:);
end

function [Y_t] = meas_ekf(X_t,u_t,params)
    b_a = X_t(1:2); b_w = X_t(3); 
    n = params.num;
    % robot pose
    P = params.X(:,n);
    P_prev = params.X_prev(:,n);
    a = (P(4:5)-P_prev(4:5))/params.dt;
    dcm = [cos(P(3)) sin(P(3));
           -sin(P(3)) cos(P(3))];
    f_b = dcm*a + b_a;
    om = P(6) + b_w;
    F_hat = params.K_F*(params.m*f_b - u_t(1:2) - params.F_hat)*params.dt;
    M_hat = params.K_M*(params.I*om - (u_t(3) - params.M_hat)*params.dt);
    Y_t = [f_b;om;F_hat;M_hat];
end

function [A] = dyn_jac_ekf(X_t,u_t,params)
    b_ax = X_t(1,:); b_ay = X_t(2,:); b_w = X_t(3,:); alpha_a = X_t(4,:); alpha_w = X_t(5,:);

    db_ax = [alpha_a,0,0,b_ax,0];
    db_ay = [0,alpha_a,0,b_ay,0];
    db_w = [0,0,alpha_w,0,b_w];
    dalpha_a = [0,0,0,1,0];
    dalpha_w = [0,0,0,0,1];

    A = [db_ax;db_ay;db_w;dalpha_a;dalpha_w];
end

% this seems wrong
function [C] = meas_jac_ekf(X_t,u_t,params)
    df_b = [1 0 0 0 0;
            0 1 0 0 0];
    dom  = [0 0 1 0 0];
    dF_hat = [params.K_F*params.m 0 0 0 0;
              0 params.K_F*params.m 0 0 0];
    dM_hat = [0 0 params.K_M*params.I 0 0];
    % ignore augmented measurement
    dF_hat = zeros(2,5);
    dM_hat = zeros(1,5);
    C = [df_b;dom;dF_hat;dM_hat];
end
