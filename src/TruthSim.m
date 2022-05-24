classdef TruthSim
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        StateDim
        MeasDim
        dt
        iter
        state_hist
        Interface
        x
        w_mean
        Q
        m
        I
        F_int
        F_ext
        M_int
        M_ext
        Cq
        D
        rho_air
    end

    methods
        function obj = TruthSim(dt,num_iter)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.StateDim = 6;
            obj.MeasDim = 6;
            obj.dt = dt;
            obj.state_hist = zeros(obj.StateDim, num_iter)
            obj.Interface = None;
            obj.initialize_properties(obj)
        end

        function obj = propagate_one_timestep(obj, t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            obj.compute_input_wrench();
            obj.computer_external_wrench();
            u = zeros(3,1);
            u(1:2) = obj.F_int + obj.F_ext;
            u(3) = obj.M_int + obj.M_ext;
            
            k1 = obj.transition(t, obj.x, u);
            k2 = obj.transition(t + obj.dt/2, obj.x + obj.dt/2 * k1, u);
            k3 = obj.transition(t + obj.dt/2, obj.x + obj.dt/2 * k2, u);
            k4 = obj.transition(t + obj.dt,   obj.x + obj.dt * k3,   u);
    
            x_next = obj.x + 1/6 * obj.dt * (k1 + 2*k2 + 2*k3 + k4);
            x_next = x_next + mvnrnd(obj.w_mean, obj.Q, 1)';
    
            obj.iter = obj.iter + 1;
            obj.state_hist(:,obj.iter) = x_next;
            obj.x = x_next;
        end


        function dx = transition( obj, t, x, u)
            state_shape = size(x);
            dx = zeros(state_shape);
            dx(1) = x(4);
            dx(2) = x(5);
            dx(3) = x(6);
    
            dx(4) = x(6)*x(5) + u(1)/obj.m;
            dx(5) = -x(6)*x(4) + u(2)/self.m;
            dx(6) = u(3)/obj.I;
        end

        function obj = compute_input_wrench(obj)
    
            u = obj.Interface.get_control_force_and_moment();
            obj.F_int(1) = u(1);
            obj.F_int(2) = u(2);
    
            obj.M_int = u(3);
        end
    
        function obj = compute_external_wrench(obj)
            
        end


        function obj = initialize_properties(obj)
            obj.x = zeros(obj.StateDim,1);
            obj.w_mean = zeros(obj.StateDim,1);
            obj.Q = obj.dt * 0.1 * eye(obj.StateDim);
            obj.m = 4.036;
            obj.I = 0.09;
            obj.F_int = zeros(2,1);
            obj.F_ext = zeros(2,1);
            obj.M_int = 0;
            obj.M_ext = 0;

            % propeller properties
            obj.Cq = 0.008;
            obj.D = 0.2286;
            obj.rho_air = 1.18;
        end

    end
end