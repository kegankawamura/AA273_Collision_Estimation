classdef TruthSim
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        StateDim
        MeasDim
        dt
        iter
        state_hist
        meas_hist
        Interface
        Map
        x
        % px, py, theta, vx, vy, omega, bax,bay,bw
        % theta is relative to inertial x-axis
        % omega is strictly in inertial z-axis
        w_mean
        Q
        m
        I
        R
        % accelerometer and gyro true bias rates, alpha<1
        alpha_b
        alpha_w
        F_int
        F_ext
        M_int
        M_ext
        mu
    end

    methods
        function obj = TruthSim(dt,num_iter, map_coord)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.StateDim = 9; %6;
            obj.MeasDim = 3;
            obj.dt = dt;
            obj.iter = 0;
            obj.state_hist = zeros(obj.StateDim, num_iter+1);
            obj.meas_hist = zeros(obj.MeasDim, num_iter+1);
            obj.Interface = NaN;
            obj.Map = ObstacleMap(map_coord);
            obj = obj.initialize_properties();
        end

        function obj = set_initial_state(obj, x)
            obj.x = reshape(x, obj.StateDim, 1);
            obj.state_hist(:,1) = obj.x;
            obj.iter = 1;
        end

        function obj = setProcessNoise(obj, w_mean, Q)
            obj.w_mean = w_mean;
            obj.Q = Q;
        end

        function obj = propagate_one_timestep(obj, t)
            obj = obj.compute_input_wrench();
            obj = obj.compute_external_wrench();
            u = zeros(3,1);
            u(1:2) = obj.F_int + obj.F_ext;
            u(3) = obj.M_int + obj.M_ext;

            k1 = obj.transition(t, obj.x, u);
            %k2 = obj.transition(t + obj.dt/2, obj.x + obj.dt/2 * k1, u);
            %k3 = obj.transition(t + obj.dt/2, obj.x + obj.dt/2 * k2, u);
            %k4 = obj.transition(t + obj.dt,   obj.x + obj.dt * k3,   u);

            % 4th order Runge-Kutta causes robot to clip out of map without
            % more sophisticated logic
            %del_x = 1/6 * (k1 + 2*k2 + 2*k3 + k4);
            del_x = k1;

            x_next = obj.x +  obj.dt * del_x;
            x_next = x_next + mvnrnd(obj.w_mean, obj.Q, 1)';

            acc = (x_next(4:5) - obj.x(4:5))/obj.dt;
            x = obj.x;

            num_tries = 0; tries_thresh = 10;
            while ~obj.Map.is_valid_loc(x_next(1:2),obj.R) && num_tries < tries_thresh
                num_tries = num_tries+1;
                [hitsWall, wall_normal, collision_loc] = obj.Map.hit_wall(x(1:2), x_next(1:2), obj.R);
                % consider including obj.x, or outputting equation of wall
                % intersection can be used to find position of robot against
                % wall.

                if hitsWall
                    disp('collide!');

                    % parallel velocity is retained
                    %
                    % perpendicular velocity goes to zero plus some noise in
                    % the opposite direction
                    % potential distance from wall to simulate "bounce"
                    %
                    % compute force (acceleration) felt by robot

                    if x(4) ~= 0
                        pre_dt = (collision_loc(1) - x(1))/x(4);
                        post_dt = obj.dt - pre_dt;
                    else
                        pre_dt = (collision_loc(2) - x(2))/x(5);
                        post_dt = obj.dt - pre_dt;
                    end

                    %disp(wall_normal)
                    %disp(obj.x(4:5))
                    %disp(dot(wall_normal, obj.x(4:5)))
                    v_perp = wall_normal * dot(wall_normal, x(4:5));
                    v_paral = x(4:5) - v_perp;

                    v_ball_contact = x(6)*obj.R*[wall_normal(2); -wall_normal(1)];

                    f_wall = -obj.mu * (v_ball_contact + v_paral);

                    post_v_paral = v_paral + obj.dt * f_wall/obj.m;

                    post_omega = x(6) - obj.dt * obj.R*(wall_normal(1)*f_wall(2) - wall_normal(2)*f_wall(1))/obj.I;

                    bounce_coeff = randn/30; % std dev is approx. 0.001

                    while abs(bounce_coeff) > 1
                        bounce_coeff = randn/30;
                    end
                    post_v_perp = -1*(1 - abs(bounce_coeff))*v_perp;

                    acc_perp = (post_v_perp-v_perp)/obj.dt; % momentum change in perp direction
                    acc_paral = f_wall/obj.m;
                    ang_acc = - obj.R*(wall_normal(1)*f_wall(2) - wall_normal(2)*f_wall(1))/obj.I;

                    acc = acc_perp + acc_paral;

                    %disp(x_next(1:2))
                    %disp(collision_loc)
                    %disp(post_v_perp)
                    %disp(v_paral)
                    x_next(4:5) = post_v_perp + post_v_paral;
                    x_next(6) = post_omega;
                    x_next(1:2) = collision_loc + post_dt * x_next(4:5);
                    x = x_next;
                else
                    acc = (x_next(4:5) - obj.x(4:5))/obj.dt;
                end
            end

            obj.iter = obj.iter + 1;
            obj.state_hist(:,obj.iter) = x_next;

            b_a = x_next(7:8);
            b_w = x_next(9);
            th = x_next(3);
            dcm = [cos(th), sin(th);
                -sin(th), cos(th)];
            acc_meas = dcm*acc + b_a;
            gyro_meas = x_next(6) + b_w;

            obj.Interface.setMeasAccel(acc_meas);
            obj.Interface.setMeasAngRate(gyro_meas);
            obj.meas_hist(:,obj.iter) = [acc_meas; gyro_meas];

            obj.x = x_next;

        end


        % biases are first order markov processes (Ornstein-Uhlenbeck to be more precise) 
        function dx = transition( obj, t, x, u)
            state_shape = size(x);
            dx = zeros(state_shape);
            dx(1:3) = x(4:6);

            %dx(4) = x(6)*x(5) + u(1)/obj.m;
            %dx(5) = -x(6)*x(4) + u(2)/obj.m;
            dx(4:5) = u(1:2)/obj.m;
            dx(6) = u(3)/obj.I;
            dx(7:8) = x(7:8)*(obj.alpha_b-1);
            dx(9) = x(9)*(obj.alpha_w-1);
        end

        function obj = compute_input_wrench(obj)

            obj.F_int = obj.Interface.getControlForce();
            obj.M_int = obj.Interface.getControlMoment();

        end

        function obj = compute_external_wrench(obj)
            % empty
        end


        function obj = initialize_properties(obj)
            obj.x = zeros(obj.StateDim,1);
            obj.w_mean = zeros(obj.StateDim,1);
            obj.Q = zeros(obj.StateDim);
            %obj.Q(obj.StateDim/2+1:obj.StateDim, obj.StateDim/2+1:obj.StateDim) = obj.dt * 0.1 * eye(obj.StateDim/2);
            obj.Q = obj.dt*blkdiag( zeros(3), 0.1*eye(3), 1*eye(2), 0.5);

            load_robot_params;
            obj.m = robot_params.m;
            obj.I = robot_params.I;
            obj.R = robot_params.R;
            obj.alpha_b = robot_params.alpha_b;
            obj.alpha_w = robot_params.alpha_w;
            obj.F_int = zeros(2,1);
            obj.F_ext = zeros(2,1);
            obj.M_int = 0;
            obj.M_ext = 0;

            obj.mu = 0.1;
            %% propeller properties
            %obj.Cq = 0.008;
            %obj.D = 0.2286;
            %obj.rho_air = 1.18;
        end

    end
    methods (Static)
        function X_t1 = point_model(X_t,u,params)
            dt = params.dt;
            F = u(1:2); M = u(3);
            X_t1(1:3,:) = X_t(1:3,:) + dt*X_t(4:6,:);
            X_t1(4:5,:) = X_t(4:5,:) + dt*F/params.m;
            X_t1(6,:) = X_t(6,:) + dt*M/params.I;
        end
        % collision model assuming zero control input used for the filter
        function x_next = collision_model(x,params)
            u = zeros(3,1);

            x_orig = x;
            is_valid = false;
            while ~is_valid
                w = mvnrnd(zeros(size(x)), params.Q_robot, 1)';
                x = x_orig+w;
                is_valid = params.Map.is_valid_loc(x(1:2),params.R);
            end

            x_next = TruthSim.point_model(x,u,params);
            %x_next = x_next + mvnrnd(zeros(size(x)), params.Q_robot, 1)';
            num_tries = 0; tries_thresh = 10;
            while ~params.Map.is_valid_loc(x_next(1:2),params.R) && num_tries < tries_thresh
                [hitsWall, wall_normal, collision_loc] = params.Map.hit_wall(x(1:2), x_next(1:2), params.R);
                if hitsWall
                    if x(4) ~= 0
                        pre_dt = (collision_loc(1) - x(1))/x(4);
                        post_dt = params.dt - pre_dt;
                    else
                        pre_dt = (collision_loc(2) - x(2))/x(5);
                        post_dt = params.dt - pre_dt;
                    end

                    v_perp = wall_normal * dot(wall_normal, x(4:5));
                    v_paral = x(4:5) - v_perp;

                    v_ball_contact = x(6)*params.R*[wall_normal(2); -wall_normal(1)];

                    f_wall = -params.mu * (v_ball_contact + v_paral);

                    post_v_paral = v_paral + params.dt * f_wall/params.m;

                    post_omega = x(6) - params.dt * params.R*(wall_normal(1)*f_wall(2) - wall_normal(2)*f_wall(1))/params.I;

                    bounce_coeff = randn/5; % greater than in truth

                    while abs(bounce_coeff) > 1
                        bounce_coeff = randn/5;
                    end
                    post_v_perp = -1*(1 - abs(bounce_coeff))*v_perp;

                    acc_perp = (post_v_perp-v_perp)/params.dt; % momentum change in perp direction
                    acc_paral = f_wall/params.m;
                    ang_acc = - params.R*(wall_normal(1)*f_wall(2) - wall_normal(2)*f_wall(1))/params.I;

                    %acc = acc_perp + acc_paral;

                    x_next(4:5) = post_v_perp + post_v_paral;
                    x_next(6) = post_omega;
                    x_next(1:2) = collision_loc + post_dt * x_next(4:5);

                else
                    %acc = (x_next(4:5) - x(4:5))/params.dt;
                end

                num_tries = num_tries+1;
                if num_tries>=tries_thresh
                    x_next = x;
                end
            end
        end
    end
end
