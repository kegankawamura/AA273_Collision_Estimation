% state:
% 1:2 [x,y,      2d position 
% 3:4  vx,vy,    2d velocity
% 5    th,       heading
% 6    om,       angular rate ccw
% 7    b_a       accelerometer bias
% 8    b_w       gyro bias
% 9:10 F_hatx,F_haty    estimated force
% 11   M_hat        estimated moment

% accelerometer takes in true state X and model map/ + system model 
% returns specific force, gyro angular rate, and updated accelerometer bias
% measurement model

function [X_1] = dyn(X,U_in,U_ex,model)
    m = model.m;
    I_z = model.I_z;
    alpha_a = model.alpha_a;
    alpha_b = model.alpha_b;
    K_F = model.K_F;
    K_M = model.K_M;
    % last measurement from imu
    f_b = model.f_b;

    dt = model.dt;

    p = [X(1:2);0];
    v = [X(3:4);0];
    th = X(5);
    om = [0;0;X(6)];
    b_a = X(7);
    b_w = X(8);
    F_hat = [X(9:10);0];
    M_hat = [0;0;X(11)];

    F_in = [U_in(1:2);0];
    F_ex = [U_ex(1:2);0];
    M_in = [0;0;U_in(3)];
    M_ex = [0;0;U_ex(3)];

    v_dot = 1/m*( cross(om,v)+(F_in+F_ex) );
    % simplified since 2d;
    w_dot = 1/I_z*(M_in+M_ex);

    % external wrench estimation
    F_hat_1 = K_F * (m*f_b - F_in - F_hat)*dt;
    M_hat_1 = K_M * (I_z*om -(M_in + M_hat)*dt );


    % euler integration 
    p_1 = p + v*dt;
    v_1 = v + v_dot*dt;
    th_1 = th + om(3)*dt;
    om_1 = om + om_dot*dt;
    b_a_1 = alpha_a*b_a ;
    b_w_1 = alpha_b*b_a ;

    %model.accel_I = (v_1-v)/dt;

    X_1 = [p_1(1:2);
           v_1(1:2);
           th_1;
           om_1(3);
           b_a_1;
           b_w_1;
           F_hat_1(1:2);
           M_hat_1(3);
           ];
       % need collision dynamics
end

function [f_b,w_b, b_a_1, b_w_1] = accelerometer(X,model)
    p = X(1:2);
    v = X(3:4);
    th = X(5);
    om = X(6);
    b_a = X(7);
    b_w = X(8);
    % updated from true (noisy) dynamics
    a = model.accel_I;
    alpha_a = model.alpha_a;
    alpha_b = model.alpha_b;

    R_th = [cos(th) sin(th);
        -sin(th) cos(th)];

    f_b = R_th*a + b_a;% + v_a;
    w_b = om + b_w;% + v_w;
    % these are actually part of the dynamics model
    b_a_1 = alpha_a*b_a ;% + v_ba;
    b_w_1 = alpha_b*b_a ;% + v_bw;
end

function [f_b,w_b, b_a_1, b_w_1] = accelerometer_sample(X,model)
    [f_b,w_b,b_a_1, b_w_1] = accelerometer(X,model);
    Q_a = model.Q_a; % noise in accelerometer
    Q_ba = model.Q_ba; % noise in accel bias
    Q_w = model.Q_w; % noise in gyro
    Q_bw = model.Q_bw; % noise in gyro bias
    v_a = sample_normal([0;0],Q_a,1);
    v_ba = sample_normal(0,Q_ba,1);
    v_w = sample_normal(0,Q_w,1);
    v_bw = sample_normal(0,Q_bw,1);

    f_b = f_b + v_a;
    w_b = w_b + v_w;
    b_a_1 = b_a_1 + v_ba;
    b_w_1 = b_w_1 + v_bw;
end

% collision pseudomeasurement
function [v_parallel] = collision_meas(X,model)
    if ~model.collision
        v_parallel = [0;0];
        return
    end
    v = X(3:4);
    F = X(9:10);
    v_parallel = v - (v'*F)/(F'*F)*F;
end

function X = sample_normal(mu,sigma,N)
    dim = size(sigma,1);
    std_sample = randn(dim,N);
    X = mu+sqrtm(sigma)*std_sample;
end
