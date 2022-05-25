classdef CPUSim
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        StateDim
        MeasDim
        dt
        iter
        state_hist
        error_hist
        meas_hist
        Interface
        x
        w_mean
        Q
        v_mean
        R
        m
        I
        F_int
        F_ext
        M_int
        M_ext
    end

    methods
        function obj = CPUSim(dt,num_iter)
            %UNTITLED5 Construct an instance of this class
            %   Detailed explanation goes here
            obj.StateDim = 11;
            obj.MeasDim = 3;
            obj.dt = dt;
            obj.state_hist = zeros(obj.StateDim, num_iter);
            obj.error_hist = zeros(obj.StateDim, obj.StateDim, num_iter);
            obj.meas_hist = zero(obj.MeasDim, num_iter);
            obj.Interface = None;
            obj.initialize_properties(obj)
        end

        function obj = run_one_timestep(obj, t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            disp('Run State Estimation...')
            disp('Output Control...')

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