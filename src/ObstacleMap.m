classdef ObstacleMap
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        corners
        coeff
        line_c
    end

    methods
        function obj = ObstacleMap(map_corners)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.corners = map_corners;
            line_coeff = zeros(length(obj.corners),3);
            % line_c = zeros(length(obj.corners),1);
            for i = 1:length(map_corners)
                ipp = mod(i, length(map_corners)) + 1;
                line_coeff(i,1) = -(obj.corners(ipp,2) - obj.corners(i,2));
                line_coeff(i,2) = obj.corners(ipp,1) - obj.corners(i,1);
                line_coeff(i,3) = obj.corners(ipp,2) * obj.corners(i,1) - obj.corners(ipp,1) * obj.corners(i,2);
                %line_c(i) = obj.corners(ipp,1) * obj.corners(i,2) - obj.corners(ipp,2) * obj.corners(i,1);
            end
            obj.coeff = line_coeff;
            disp(obj.coeff)
            %obj.line_c = line_c;
        end

        function [hitsWall, min_wall_normal] = hit_wall(obj,r)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            count_intersection = 0;
            ray = zeros(1,3);
            ray(1:2) = [-0.1924    0.8886]; % randn(1,2);
            ray(3) = -1 * ray(1:2) * reshape(r,[2,1]);
            %disp(ray)
            intersecting_walls = [];

            for i = 1:length(obj.coeff)
                l = obj.coeff(i,:);
                %disp('Loop')
                %disp(l)
%                C(1,:) = obj.line_c(i);
                P_h = cross(l, ray);
                %disp(P_h)
                p = [P_h(1)/P_h(3), P_h(2)/P_h(3)];

                ipp = mod(i, length(obj.coeff)) + 1;
                %disp(p)
                %disp(obj.corners(i,:))
                if p(1) <= r(1) && obj.num_between(p(1), obj.corners(i,1), obj.corners(ipp,1)) && obj.num_between(p(2), obj.corners(i,2), obj.corners(ipp,2))
                    count_intersection = count_intersection + 1;
                    disp(i)
                end
            end

            hitsWall = false;
            min_wall_normal = [0, 0];
            if mod(count_intersection, 2) == 0
                hitsWall = true;
                r_h = [r(1), r(2), 1];
                disp(r_h)
                min_dist = realmax;
                min_wall = 0;
                min_dir = 0;
                for i = 1:length(obj.coeff)
                    %disp('Loop')
                    a = obj.coeff(i,1);
                    b = obj.coeff(i,2);
                    T = obj.coeff(i,:) * r_h';
                    %disp(T)
                    d = abs(T)/sqrt(a^2 + b^2);
                    sgn = sign(T);
                    wall_normal = obj.coeff(i,1:2);
                    wall_normal = -1 * sgn * wall_normal/norm(wall_normal);
                    xp = r(1) - sgn*d*a/sqrt(a^2 + b^2);
                    yp = r(2) - sgn*d*b/sqrt(a^2 + b^2);
                    if ~obj.num_between(xp, obj.corners(i,1), obj.corners(ipp,1)) && ~obj.num_between(yp, obj.corners(i,2), obj.corners(ipp,2))
                        d = sqrt(min(sum((r - obj.corners(i,:)).^2), sum((r - obj.corners(ipp,:)).^2)));
                    end

                    %disp(d)

                    if d < min_dist
                        min_dist = d;
                        min_wall = i;
                        min_dir = -1*sgn;
                        min_wall_normal = wall_normal;
                    end
                end
            end
        end

        function is_between = num_between(obj, x, a, b)
            is_between = (x >= a && x <= b) || (x <= a && x >= b);
        end
    end
end