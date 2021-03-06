classdef ObstacleMap < handle
    properties
        corners
        coeff
        line_c
    end

    methods
        function obj = ObstacleMap(map_corners)
            if nargin ==0
                return
            end
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

        function [hitsWall, wall_normal, collision_loc] = hit_wall(obj,r_pre, r, R)
            count_intersection = 0;
            ray = zeros(1,3);
            rel_vec = r_pre - r;
            ray(1:2) = [rel_vec(2), -rel_vec(1)]; % randn(1,2);
            ray(3) = -1 * ray(1:2) * reshape(r,[2,1]);
            %disp(ray)
            closest_intersection = r;
            closest_wall = 0;

            hitsWall = false;
            for i = 1:length(obj.coeff)
                l = obj.coeff(i,:);
                denom = sqrt(l(1)^2 + l(2)^2);
                P_h = cross(l, ray);
                p = [P_h(1)/P_h(3); P_h(2)/P_h(3)];

                A = [l(1) l(2);
                    ray(1) ray(2)];

                b1 = [R*denom - l(3); -ray(3)];
                b2 = [-R*denom - l(3); -ray(3)];

                ipp = mod(i, length(obj.coeff)) + 1;
                if (obj.num_between(p(1), r_pre(1), r(1)) && ...
                        obj.num_between(p(2), r_pre(2), r(2)) && ...
                        obj.num_between(p(1), obj.corners(i,1), obj.corners(ipp,1)) && ...
                        obj.num_between(p(2), obj.corners(i,2), obj.corners(ipp,2))) || ...
                        R >= abs(l(1)*r(1) + l(2)*r(2) +l(3))/denom

                    hitsWall = true;
                    count_intersection = count_intersection + 1;
                    %disp(p)
                    %disp(r_pre)
                    %disp(r)
                    %disp(closest_intersection)
                    %disp(norm(p-r_pre))
                    %disp(norm(closest_intersection - r_pre))
                    if norm(p-r_pre) <= norm(closest_intersection - r_pre)
                        closest_intersection = p;
                        closest_wall = i;
                    elseif closest_intersection ==r
                        closest_intersection = p;
                        closest_wall = i;
                    end
                end
            end

            wall_normal = [0; 0];
            p_rollback = [NaN;NaN];
            if hitsWall
                %disp('Hit Wall')
                %disp(closest_wall)
                r_h = [r(1), r(2), 1];
                a = obj.coeff(closest_wall,1);
                b = obj.coeff(closest_wall,2);
                c = obj.coeff(closest_wall,3);
                denom = sqrt(a^2 + b^2);
                rel_vec = r_pre - closest_intersection;
                sgn = 1;
                if rel_vec(1)*a + rel_vec(2)*b < 0
                    sgn = -1;
                end
                wall_normal = obj.coeff(closest_wall,1:2);
                wall_normal = sgn * wall_normal/norm(wall_normal);

                A = [a b;
                    ray(1) ray(2)];

                b1 = [R*denom - c; -ray(3)];
                b2 = [-R*denom - c; -ray(3)];

                p_rollback = A\b1;
                if ~obj.num_between(p_rollback(1), r_pre(1), closest_intersection(1)) || ...
                        ~obj.num_between(p_rollback(2), r_pre(2), closest_intersection(2))
                    p_rollback = A\b2;
                end

            end

            collision_loc = p_rollback;
            wall_normal = reshape(wall_normal, 2,1);
        end





        function [hitsWall, min_wall_normal] = outside_map(obj, r)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            count_intersection = 0;
            ray = zeros(1,3);
            ray(1:2) = randn(1,2);
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
                    %disp(i)
                end
            end

            hitsWall = false;
            min_wall_normal = [0, 0];
            if mod(count_intersection, 2) == 0
                hitsWall = true;
                r_h = [r(1), r(2), 1];
                %disp(r_h)
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

        function valid = is_valid_loc(obj,r,R)
            valid = obj.is_inside(r)  &&  ~obj.inside_wall(r,R);
        end

        function inside = is_inside(obj, r)
            count_intersection = 0;
            ray = zeros(1,3);
            ray(1:2) = randn(1,2);
            %ray(1:2) = [0.49,0.58];
            ray(3) = -1 * ray(1:2) * reshape(r,[2,1]);

            for i = 1:length(obj.coeff)
                l = obj.coeff(i,:);
                P_h = cross(l, ray);
                p = [P_h(1)/P_h(3), P_h(2)/P_h(3)];
                ipp = mod(i, length(obj.coeff)) + 1;
                for j = 1:2
                    if obj.corners(i,j) == obj.corners(ipp,j)
                        p(j) = obj.corners(i,j);
                    end
                end
                if p(1) <= r(1) && ...
                        obj.num_between(p(1), obj.corners(i,1), obj.corners(ipp,1)) && ...
                        obj.num_between(p(2), obj.corners(i,2), obj.corners(ipp,2))
                    count_intersection = count_intersection + 1;
                end
            end

            inside = true;

            if mod(count_intersection, 2) == 0
                inside = false;
            end
        end

        function [P] = sample_pt(obj,mu,sigma,R)
            P = mvnrnd(mu,sigma,1)';
            while ~obj.is_valid_loc(P(1:2),R)
                P = mvnrnd(mu,sigma,1)';
            end
            return
        end

function in_wall = inside_wall(obj, r, R)
            r_h = [r(1); r(2); 1];
            in_wall = false;
            for i = 1:length(obj.coeff)
                l = obj.coeff(i,:);
                a = l(1);
                b = l(2);
                c = l(3);
                denom = sqrt(a^2 + b^2);
                d = abs(l*r_h)/denom;

                xp = (b^2*r(1) - a*b*r(2) - a*c)/denom^2;
                yp = (a^2*r(2) - a*b*r(1) - b*c)/denom^2;

                ipp = mod(i, length(obj.coeff)) + 1;
                if d < R && ...
                        obj.num_between(xp, obj.corners(i,1), obj.corners(ipp,1)) && ...
                        obj.num_between(yp, obj.corners(i,2), obj.corners(ipp,2)) 
                    in_wall = true;
                    break
                end
            end
        end
        function is_between = num_between(obj, x, a, b)
            is_between = (x >= a && x <= b) || (x <= a && x >= b);
        end

        function plt = plot_map(obj)
                plt = plot( [obj.corners(:,1);obj.corners(1,1)],...
                    [obj.corners(:,2);obj.corners(1,2)]);
        end
    end
end
