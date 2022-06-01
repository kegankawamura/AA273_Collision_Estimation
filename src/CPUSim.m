classdef CPUSim
    properties
        StateDim
        MeasDim
        dt
        iter
        state_hist
        error_hist
        meas_hist
        Interface
        PF
        EKF
        x
        w_mean
        Q
        v_mean
        R
        m = 4.036;
        I = 0.09;
        F_int = [0;0];
        F_ext
        M_int = 0;
        M_ext
    end

    methods
        function obj = CPUSim(dt,num_iter, system)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.StateDim = 11;
            obj.MeasDim = 3;
            obj.dt = dt;
            obj.state_hist = zeros(obj.StateDim, num_iter);
            obj.error_hist = zeros(obj.StateDim, obj.StateDim, num_iter);
            obj.meas_hist = zeros(obj.MeasDim, num_iter);
            obj.PF = ParticleFilter(system);
            obj.EKF = EKF(system);
            obj.Interface = NaN;
            obj = obj.initialize_properties(system);
        end

        function obj = initialize_filter(obj, mu_0, sigma_0)
            mu_PF = [mu_0(1:6); mu(9:11)];
            mu_EKF = mu_0(7:8);
            
            sigma_PF = [sigma_0(1:6, 1:6), sigma_0(1:6, 9:11);
                        sigma_0(9:11, 1:6), sigma_0(9:11, 9:11);];
            sigma_EKF = sigma_0(7:8, 7:8);

            obj.PF.initialize(mu_PF, sigma_PF);
            obj.EKF.initialize(mu_EKF, sigma_EKF);
            
        end

        function obj = run_one_timestep(obj, t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            % Collect measurements
            y_t = zeros(3,1);
            y_t(1:2) = obj.Interface.getMeasAccel();
            y_t(3) = obj.Interface.getMeasAngRate();

            u_t = [obj.F_int; obj.M_int];

            %disp('Run State Estimation...')
            obj.PF.predict(u_t);
            obj.EKF.predict(u_t);

            obj.PF.update(u_t, y_t);
            obj.EKF.update(u_t, y_t);

            %disp('Output Control...')

            % compute controls u_tpp
            u_tpp = zeros(3,1);
            
            obj.Interface.setControlForce(u_tpp(1:2));
            obj.Interface.setControlMoment(u_tpp(3));

            obj.F_int = u_tpp(1:2);
            obj.M_int = u_tpp(3);

        end
      
        function obj = initialize_properties(obj)
            obj.x = zeros(obj.StateDim,1);
            obj.w_mean = zeros(obj.StateDim,1);
            obj.Q = obj.dt * 0.1 * eye(obj.StateDim);
            obj.v_mean = zeros(obj.MeasDim,1);
            obj.R = 0.1 * eye(obj.MeasDim);
            obj.m = 4.036;
            obj.I = 0.09;
            obj.F_int = zeros(2,1);
            obj.F_ext = zeros(2,1);
            obj.M_int = 0;
            obj.M_ext = 0;
        end
    
    end
end
