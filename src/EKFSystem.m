% this is just a dumb idea
classdef EKFSystem < handle & dynamicprops
    % used to efficiently share parameters between the ekfs
    properties 
        system;
    end
    methods
        function obj = EKFSystem(ekf_system)
            obj.system = ekf_system;
            % so dumb
            fs = fields(ekf_system);
            for i=1:numel(fs)
                f = fs{i};
                obj.addprop(f);
                obj.(f) = ekf_system.(f);
            end
        end

    end
end
