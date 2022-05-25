classdef EKF < handle
    properties
        type = 'EKF'
        system
        % state of filter
        mu_current;
        sigma_current;
        last_step = '';

        % results of filter, after update step
        Mu
        Sigma
    end
    methods
        function obj = EKF(system)
            if nargin ==0
                return
            end
            obj.system = system;
        end

        function [] = initialize(obj,mu_0, sigma_0)
            obj.Mu = mu_0;
            obj.Sigma = sigma_0;
            obj.mu_current = mu_0;
            obj.sigma_current = sigma_0;

            system = obj.system;
        end
        function [mu_p,sigma_p] = predict(obj,u_t)
            system = obj.system;
            mu_0 = obj.mu_current;
            sigma_0 = obj.sigma_current;
            mu_p = system.dyn(mu_0,u_t,system.params);
            A = system.A(mu_0,u_t,system.params);
            sigma_p= A*sigma_0*A' + system.Qmodel;

            obj.mu_current = mu_p;
            obj.sigma_current = sigma_p;
            obj.last_step = 'predict';
        end

        function [mu_1,sigma_1] = update(obj,u_t,y_t)
            system = obj.system;
            params = system.params;
            mu_p = obj.mu_current;
            sigma_p = obj.sigma_current;
            C = system.C(mu_p,u_t,params);
            A = system.A(mu_p,u_t,params);
            Q = system.Qmodel;
            R = system.Rmodel;

            K_t = sigma_p*C'*(C*sigma_p*C'+R)^-1;

            %innov = y_t-system.meas(mu_p,u_t,params)
            %if norm(innov)>2
            %    %keyboard
            %end
            mu_1 = mu_p + K_t*( y_t-system.meas(mu_p,u_t,params) );
            sigma_1 = sigma_p - K_t*C*sigma_p;

            obj.mu_current = mu_1;
            obj.sigma_current = sigma_1;
            obj.last_step = 'update';

            obj.Mu(:,end+1) = mu_1;
            obj.Sigma(:,:,end+1) = sigma_1;

        end
    end
end
