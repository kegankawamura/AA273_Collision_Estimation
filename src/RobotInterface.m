classdef RobotInterface
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties

        MeasDim
        v_mean
        R

        % add measurements
        IMUAccel
        IMUAngRate

        % add control outputs
        Force
        Moment

    end

    methods

        function obj = RobotInterface( )
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.MeasDim = 3;

            obj.v_mean = zeros(obj.MeasDim,1);
            obj.R = 0.1 * eye(obj.MeasDim);

            obj.Force = [0; 0];
            obj.Moment = 0;
        end

        function obj = setMeasAccel(obj,accel)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.IMUAccel = accel;
        end

        function accel = getMeasAccel(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            accel = obj.IMUAccel + mvnrnd(obj.v_mean(1:2), obj.R(1:2,1:2), 1)';
        end

        function obj = setMeasAngRate(obj,angRate)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.IMUAngRate = angRate;
        end

        function angRate = getMeasAngRate(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            angRate = obj.IMUAngRate + mvnrnd(obj.v_mean(3), obj.R(3,3), 1)';
        end

        function obj = setControlForce(obj,force)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.Force = force;
        end

        function force = getControlForce(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            force = obj.Force;
        end

        function obj = setControlMoment(obj,moment)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.Moment = moment;
        end

        function moment = getControlMoment(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            moment = obj.Moment;
        end
    end
end