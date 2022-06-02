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
        F_hat_hist = [0;0];
        M_hat_hist = [0];
        ekf_mu_0; ekf_sigma_0;
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
            obj.ekfs(obj.N) = EKF(obj.ekf_system);
            [obj.ekfs.system] = deal(obj.ekf_system);
            obj.ekfs.initialize(mu,sigma);
            obj.ekf_mu_0 = mu;
            obj.ekf_sigma_0 = sigma;
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
            F = K_F * (m*f_b - u_t(1:2) - F_hat_prev)*dt;
            M = K_M * (I_z*w_b -(u_t(3) + M_hat_prev)*dt);

            % augment with pseudomeasurement
            y_t_aug = [y_t;F;M];

            % updating ekfs need a way to index ekfs
            nums = num2cell([1:obj.N]);
            [obj.ekfs.num] = deal(nums{:});
            % obj.ekfs.F_hat_prev = F_hat_prev;
            obj.ekfs.update(u_t,y_t_aug);

            % particle filter update a-la fast slam
            % idea: use F and M and the wall normals from collisions as additional weighting, since force estimates dont work super well in ekf
            W_bar = [obj.ekfs.p_yj];
            alive = (W_bar>1e-10);
            if (sum(alive)>=0.25*obj.N)
                obj.ekf_mu_0 = mean([obj.ekfs(alive).mu_current].*W_bar(alive),2);
                obj.ekf_sigma_0 = mean(...
                    reshape([obj.ekfs(alive).sigma_current],5,5,[])... 
                    .*reshape(W_bar(alive),1,1,[]),3);
            end

            if sum(W_bar==0)>0.75*obj.N
                fprintf('%i particles were killed, resetting\n',obj.N-sum(alive))
                % too many ekfs are implausible, wat do
                % try to reset them to the prior
                for i = 1:obj.N
                    if ~alive(i)
                        obj.ekfs(i).mu_current = obj.ekf_mu_0;
                        obj.ekfs(i).sigma_current = obj.ekf_sigma_0;
                        W_bar(i) = 1/obj.N;
                        obj.pf.particles.X(:,i) = obj.pf.Mu;
                    end
                end
                keyboard
            end
            W_hat = W_bar.*obj.pf.particles.W;
            W_t = W_hat/sum(W_hat);
            obj.pf.particles.W = W_t;
            % need to create a resample condition
            %resample = std(W_t)>0.03/sqrt(obj.N);
            resample = true;
            if resample
                [idx_resample] = obj.pf.importance_resample();
                % if only it were this easy
                %obj.ekfs = obj.ekfs(idx_resample);
                ekfs_copy = copy(obj.ekfs);
                for i=1:numel(obj.ekfs)
                    obj.ekfs(i) = copy(ekfs_copy(idx_resample(i)));
                end
            end
            ax = gca;
            if numel(ax.Children)<2
                obj.pf.plot_particles(1,2,1);
                obj.pf_system.params.Map.plot_map();
            else
                set(ax.Children(3),'XData',obj.pf.particles.X(1,:));
                set(ax.Children(3),'YData',obj.pf.particles.X(2,:));
                drawnow;
            end

                if any(obj.pf.particles.X(1:2)>10) || any(obj.pf.particles.X(1:2)<0) 
                    keyboard
                end

            obj.F_hat_hist(:,end+1) = F;
            obj.M_hat_hist(:,end+1) = M;
        end
    end
end
