classdef CIOFilter < handle
    % specific filter for the CIO localization problem
    %   particle filter that maintains pose
    %   array of ekfs that are each associated with a particle filter
    %   force estimate which is passed as an additional measurement to each filter
    properties 
        N;
        % struct of known robot parameters
        params;
        % a normal struct
        pf_system;
        % instantiation of class
        ekf_system;
        % particle filter of pose
        pf;
        % ekfs associated with each particle 
        ekfs = EKF();
        F_hat_hist;
        M_hat_hist;

    end
    methods
        function obj = CIOFilter(pf_system,ekf_system,robot_params)
            obj.pf_system = pf_system;
            obj.ekf_system = ekf_system;
            obj.pf = ParticleFilter(pf_system);
            obj.params = robot_params;
        end
        % initializes with prior particles of robot pose and mean+sigma of bias states  and drift rates
        % particles are struct particles.X, particles.W
        function [] = initialize(obj,particles,mu,sigma)
            obj.pf.initialize_particles(particles);
            obj.N = numel(particles.W);
            obj.ekfs(1:obj.N) = copy(EKF(obj.ekf_system));
            obj.ekfs.initialize(mu,sigma);
            keyboard
        end

        function [] = predict(obj,u_t)
            obj.ekfs(1).system.params.X_prev = obj.pf.particles.X;
            obj.pf.predict(u_t);
            obj.ekfs.predict(u_t);
            % give ekfs the list of particle states 
            % which they should all share supposedly
            obj.ekfs(1).system.params.X = obj.pf.particles.X;
            obj.ekfs(1).system.params.F_hat = obj.F_hat_hist(:,end);
            obj.ekfs(1).system.params.M_hat = obj.M_hat_hist(:,end);

        end

        % replaces the particle filter update step
        function [] = update(obj,u_t,y_t)
            f_b = y_t(1:2);
            w_b = y_t(3);
            % calculate estimated F and M
            K_F = obj.params.K_F;
            K_M = obj.params.K_M;
            m = obj.params.m;
            I_z = obj.params.I;
            dt = obj.params.dt;
            F_hat_prev = obj.F_hat_hist(:,end);
            M_hat_prev = obj.M_hat_hist(:,end);
            F = K_F * (m*f_b - F_in - F_hat_prev)*dt;
            M = K_M * (I_z*w_b -(M_in + M_hat_prev)*dt);

            % augment with pseudomeasurement
            y_t_aug = [y_t;F;M];

            % updating ekfs need a way to index ekfs
            nums = num2cell([1:obj.N]);
            [obj.ekfs.num] = deal(nums{:});
            % obj.ekfs.F_hat_prev = F_hat_prev;
            obj.ekfs.update(u_t,y_t_aug);

            % particle filter update a-la fast slam
            % idea: use F and M and the wall normals from collisions as additional weighting, since force estimates dont work super well in ekf
            W_bar = prod(reshape([obj.particles.ekf.p_yj],[],obj.N));
            W_hat = W_bar.*obj.particles.W;
            W_t = W_hat/sum(W_hat);
            obj.particles.W = W_t;
            % need to create a resample condition
            %resample = std(W_t)>0.03/sqrt(obj.N);
            resample = true;
            if resample
                [idx_resample] = obj.pf.importance_resample();
                obj.ekfs = obj.ekfs(idx_resample);
            end

            obj.F_hat_hist(:,end+1) = F;
            obj.M_hat_hist(:,end+1) = M;
        end
    end
end
