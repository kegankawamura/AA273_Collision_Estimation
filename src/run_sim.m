
% Final Time
tf = 60;

% Integration Timestep
dt = 0.1;

iter = tf/dt;

% Define map coordinates
map_coord = [   -1,  -1;
               10,  0;
               10, 10;
                0, 10];

% Intialize Truth and Discrete sims
truth = TruthSim(dt, iter, map_coord);
discrete = CPUSim(dt, iter);

inter = RobotInterface();
truth.Interface = inter;
discrete.Interface = inter;

% Set custom initial conditions
truth = truth.set_initial_state([5, 5, pi/4, 2, 0, 0]);
Q = zeros(6);
Q(4:5, 4:5) = dt * 0.01 * eye(2);
Q(6,6) = dt * 0.001;
truth = truth.setProcessNoise(zeros(6,1), Q);

% Simulate
for i = 1:iter
    t = dt * i;
    truth = truth.propagate_one_timestep(t);
    discrete = discrete.run_one_timestep(t);
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

%%

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

