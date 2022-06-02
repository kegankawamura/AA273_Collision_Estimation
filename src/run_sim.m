clear;
% Final Time
tf = 60;
load_robot_params;

% Integration Timestep
dt = robot_params.dt;

iter = tf/dt;

% Define map coordinates
map_coord = [   0,  0;
               3,  2;
               5, 5;
                3, 4];

% output system
%define_system

% Intialize Truth and Discrete sims
truth = TruthSim(dt, iter, map_coord);
discrete = CPUSim(dt, iter);

% create filter
mu_0 = zeros(6,1);
mu_0(1:2) = [4;4];
mu_0(3:4) = [2;0];
mu_0(5:6) = [pi/4;0];
sigma = eye(6);
filter = setup_filter(robot_params,truth.Map,mu_0,sigma);
discrete.filter = filter;

inter = RobotInterface();
truth.Interface = inter;
discrete.Interface = inter;


% Set custom initial conditions
truth = truth.set_initial_state([4, 4, pi/4, 2, 0, 0, 0.5,0.5,0.2]);
Q = zeros(9);
Q = blkdiag( zeros(3), dt*0.01*eye(2), dt*0.001, 0.3*eye(3) );
%Q(4:5, 4:5) = dt * 0.01 * eye(2);
%Q(6,6) = dt * 0.001;
truth = truth.setProcessNoise(zeros(9,1), Q);

% state:
% 1:2 [x,y,      2d position 
% 3    th,       heading
% 4:5  vx,vy,    2d velocity
% 6    om,       angular rate ccw
% 7:8  b_ax,b_ay accelerometer bias
% 9    b_w       gyro bias
mu_0 = zeros(9,1);
%mu_0 = [5;5;2;0;pi/4;0;0;0;0;0;0];
mu_0(1:2) = [5;5];
mu_0(3:4) = [2;0];
mu_0(5:6) = [pi/4;0];

sigma_0 = zeros(11);
sigma_0(3:4,3:4) = 0.1*dt*eye(2);
sigma_0(6,6) = 0.01*dt;
sigma_0(7:8,7:8) = 0.001*dt*eye(2);
sigma_0(9:11,9:11) = 0.1*dt*eye(3);

%discrete = discrete.initialize_filter(mu_0, sigma_0);


figure(1); clf; hold on;
    ax = gca;
% Simulate
tic;
for i = 1:iter
    if ~mod(i,10)
        t_elapse = toc;
        fprintf('completed %i / %i in %.1f sec \n',i,iter,t_elapse);
        fprintf('true position is [%.2f;%.2f]\n',truth.x(1),truth.x(2));
        tic;
    end
    t = dt * i;
    truth = truth.propagate_one_timestep(t);
    discrete = discrete.run_one_timestep(t);

    if numel(ax.Children)<3
        scatter(truth.x(1),truth.x(2),50,'o','LineWidth',2);
    else
        set(ax.Children(1),'XData',truth.x(1));
        set(ax.Children(1),'YData',truth.x(2));
        drawnow;
    end
end


% Collect Data
truth_state = truth.state_hist;
est_state = discrete.state_hist;
est_error = discrete.error_hist;

meas_wo_noise = truth.meas_hist;
measurements = discrete.meas_hist;

T = [0:dt:iter*dt];


% Plot output
show_map_coord = [map_coord; map_coord(1,:)];
close all
figure
hold on
plot(show_map_coord(:,1), show_map_coord(:,2))
plot(truth_state(1,:), truth_state(2,:))
hold off

figure
subplot(3,1,1)
plot(T, truth_state(1,:))
title('Pose: x, y, theta')
subplot(3,1,2)
plot(T, truth_state(2,:))
subplot(3,1,3)
plot(T, truth_state(3,:))

figure
subplot(3,1,1)
plot(T, truth_state(4,:))
title('Pose Velocity: vx, vy, omega')
subplot(3,1,2)
plot(T, truth_state(5,:))
subplot(3,1,3)
plot(T, truth_state(6,:))

figure
subplot(3,1,1)
plot(T, meas_wo_noise(1,:))
title('Measurements w/o Noise: ax, ay, omega')
subplot(3,1,2)
plot(T, meas_wo_noise(2,:))
subplot(3,1,3)
plot(T, meas_wo_noise(3,:))

