
tf = 100;

dt = 0.1;

iter = tf/dt;

map_coord = [   -1,  -1;
               10,  0;
               10, 10;
                0, 10];

truth = TruthSim(dt, iter, map_coord);
discrete = CPUSim(interface, dt, iter);

interface = RobotInterface();
truth.Interface = interface;
discrete.Interface = interface;

for i = 1:iter
    t = dt * i;
    truth.propagate_one_timestep(t)
    discrete.run_one_timestep(t)
end

truth_state = truth.state_hist;
est_state = discrete.state_hist;
est_error = discrete.error_hist;

meas_wo_noise = truth.meas_hist;
measurements = discrete.meas_hist;