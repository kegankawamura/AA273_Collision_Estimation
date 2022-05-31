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
        % px, py, theta, vx, vy, omega
        % theta is relative to inertial x-axis
        % omega is strictly in inertial z-axis
        w_mean
        Q
        m
        I
        F_int
        F_ext
        M_int
        M_ext
    end

    methods
        function obj = TruthSim(dt,num_iter, map_coord)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.StateDim = 6;
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
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
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
    
            [hitsWall, wall_normal, collision_loc] = obj.Map.hit_wall(obj.x(1:2), x_next(1:2));
            % consider including obj.x, or outputting equation of wall
            % intersection can be used to find position of robot against
            % wall.

            if hitsWall

                % parallel velocity is retained
                %
                % perpendicular velocity goes to zero plus some noise in
                % the opposite direction
                % potential distance from wall to simulate "bounce"
                %
                % compute force (acceleration) felt by robot

                if obj.x(4) ~= 0
                    pre_dt = (collision_loc(1) - obj.x(1))/obj.x(4);
                    post_dt = obj.dt - pre_dt;
                else
                    pre_dt = (collision_loc(2) - obj.x(2))/obj.x(5);
                    post_dt = obj.dt - pre_dt;
                end

                %disp(wall_normal)
                %disp(obj.x(4:5))
                %disp(dot(wall_normal, obj.x(4:5)))
                v_perp = wall_normal * dot(wall_normal, obj.x(4:5));
                v_paral = obj.x(4:5) - v_perp;

                bounce_coeff = randn/30; % std dev is approx. 0.001

                post_v_perp = -1*abs(bounce_coeff)*v_perp;

                acc_avg = (post_v_perp-v_perp)/obj.dt; % momentum change in perp direction


                %disp(x_next(1:2))
                %disp(collision_loc)
                %disp(post_v_perp)
                %disp(v_paral)
                x_next(4:5) = post_v_perp + v_paral;
                x_next(1:2) = collision_loc + post_dt * x_next(4:5) + wall_normal*1e-6; % add bias to prevent clipping through wall
                
                obj.iter = obj.iter + 1;
                obj.state_hist(:,obj.iter) = x_next;

                obj.Interface.setMeasAccel(acc_avg);
                obj.Interface.setMeasAngRate(x_next(6));
                obj.meas_hist(:,obj.iter) = [acc_avg; x_next(6)];

            else
                obj.iter = obj.iter + 1;
                obj.state_hist(:,obj.iter) = x_next;
    
                acc_noisy = (x_next(4:5) - obj.x(4:5))/obj.dt;
    
                obj.Interface.setMeasAccel(acc_noisy);
                obj.Interface.setMeasAngRate(x_next(6));.
                
                obj.meas_hist(:,obj.iter) = [acc_noisy; x_next(6)];

            end

            obj.x = x_next;

        end


        function dx = transition( obj, t, x, u)
            state_shape = size(x);
            dx = zeros(state_shape);
            dx(1) = x(4);
            dx(2) = x(5);
            dx(3) = x(6);
    
            dx(4) = x(6)*x(5) + u(1)/obj.m;
            dx(5) = -x(6)*x(4) + u(2)/obj.m;
            dx(6) = u(3)/obj.I;
        end

        function obj = compute_input_wrench(obj)

            obj.F_int = obj.Interface.getControlForce();
            obj.M_int = obj.Interface.getControlMoment();

        end
    
        function obj = compute_external_wrench(obj)
            
        end


        function obj = initialize_properties(obj)
            obj.x = zeros(obj.StateDim,1);
            obj.w_mean = zeros(obj.StateDim,1);
            obj.Q = zeros(obj.StateDim);
            obj.Q(obj.StateDim/2+1:obj.StateDim, obj.StateDim/2+1:obj.StateDim) = obj.dt * 0.1 * eye(obj.StateDim/2);
            obj.m = 4.036;
            obj.I = 0.09;
            obj.F_int = zeros(2,1);
            obj.F_ext = zeros(2,1);
            obj.M_int = 0;
            obj.M_ext = 0;

            %% propeller properties
            %obj.Cq = 0.008;
            %obj.D = 0.2286;
            %obj.rho_air = 1.18;
        end

    end
end