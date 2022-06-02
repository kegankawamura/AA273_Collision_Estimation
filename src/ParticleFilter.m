classdef ParticleFilter < handle
    properties
        type = 'PF'
        N = 2000;
        system
        % state of filter
        particles
        last_step = '';

        % results of filter, after update step
        Mu
        Sigma
        resample_count = 0;
    end
    methods
        function obj = ParticleFilter(system)
            obj.system = system;
        end

        % ekf
        
        function [] = initialize(obj,mu_0, sigma_0)
            % generate prior particles from normal prior
            N = obj.N;
            X = ParticleFilter.sample_normal(mu_0,sigma_0,N);
            W = 1/N * ones(1,N);
            particles.X = X;
            particles.W = W;
            obj.particles = particles;
            obj.Mu = mu_0;
            obj.Sigma = sigma_0;
        end

        function [] = initialize_particles(obj,particles)
            obj.N = numel(particles.W);
            obj.particles = particles;
            W = particles.W; X = particles.X;
            obj.Mu = sum(W.*X,2);
            obj.Sigma = (W.*(X-obj.Mu))*(X-obj.Mu)';
        end

        function [] = predict(obj,u_t)
            system = obj.system;
            params = system.params;
            particles = obj.particles;
            X_t = particles.X;

            % system.transition defines a stochastic transition model (e.g. simplified collision model) 
            X_t1 = system.transition(X_t,u_t,params);
            %W_noise = ParticleFilter.sample_normal(...
            %    zeros(size(X_t(:,1))),...
            %    system.Qmodel,obj.N);
            %X_t1 = system.dyn(X_t,u_t,params) + W_noise;

            obj.particles.X = X_t1;
            obj.last_step = 'predict';
        end

        % unused function
        function [] = update(obj,u_t,y_t)
            system = obj.system;
            params = system.params;
            
            particles = obj.particles;
            X_t = obj.particles.X;

            Y_hat = system.meas(X_t,u_t,params);
            W_hat = ParticleFilter.prob_normal(...
                Y_hat,y_t,system.Rmodel).*obj.particles.W;
            W_t = W_hat/sum(W_hat);
            obj.particles.W = W_t;

            % condition for resampling, tune
            resample = std(W_t)>0.03/sqrt(obj.N);
            if resample
                obj.importance_resample();
            end

            obj.last_step = 'update';
            % compute empirical mean and cov
            W = obj.particles.W; X = obj.particles.X;
            mu = sum(W.*X,2);
            sigma = (W.*(X-mu))*(X-mu)';
            obj.Mu(:,end+1) = mu;
            obj.Sigma(:,:,end+1) = sigma;
        end

        function [idx] = importance_resample(obj);
            X = obj.particles.X;
            W = obj.particles.W;
            pdf = cumsum(W);
            r = rand(1,obj.N);
            idx = sum(r>pdf',1)+1;
            X_new = X(:,idx);
            obj.particles.X = X_new;
            obj.particles.W = 1/obj.N*ones(1,obj.N);
            obj.resample_count = obj.resample_count+1;
        end

        function plt = plot_histogram(obj,n,figno)
            % creates a copy of the particle filter with resampled particles 
            pf_copy = ParticleFilter(obj.system);
            pf_copy.particles = obj.particles;
            pf_copy.importance_resample();
            x = pf_copy.particles.X(n,:);
            figure(figno);
            plt = histogram(x);
            end

        function plt = plot_particles(obj,n1,n2,figno)
            x1 = obj.particles.X(n1,:);
            x2 = obj.particles.X(n2,:);
            w  = obj.particles.W;
            figure(figno);
            plt = scatter(x1,x2,w*obj.N*50,'o','filled');
            plt.MarkerFaceAlpha = 0.15;
            axis equal;
        end
            
    end
    methods (Static)
        function X = sample_normal(mu,sigma,N)
            dim = size(sigma,1);
            std_sample = randn(dim,N);
            X = mu+sqrtm(sigma)*std_sample;
        end
        % normal pdf of drawing x from N(mu,sigma)
        function P = prob_normal(x,mu,sigma)
            n = numel(mu);
            eta = 1./sqrt( (2*pi)^n * det(sigma) );
            P = eta.*exp( -0.5* dot((x-mu), sigma^-1 * (x-mu) ,1) );
        end
    end
end
