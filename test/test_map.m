addpath('../src/')

corners = [[0,0]; [0,2]; [-2,2]; [-2, 4]; [2, 6]; [1,4]; [3,1]];
show_corner = [corners; corners(1,:)];

close all
figure
hold on 
plot(show_corner(:,1), show_corner(:,2))

points = [[1,1]; [-1.5,3]; [-1,1.9]; [-2,3]];
scatter(points(:,1), points(:,2))

map = ObstacleMap(corners);

%[hits_wall, normal] = map.outside_map([1,1]);
%[hits_wall, normal] = map.outside_map([-1.5,3]);
%[hits_wall, normal] = map.outside_map([-1,1.9]);
%[hits_wall, normal] = map.outside_map([-2,3]);


[hitsWall, wall_normal, collision_loc] = map.hit_wall([2,1], [1,1]);
disp(hitsWall)
disp(collision_loc)
scatter(collision_loc(1), collision_loc(2))
[hitsWall, wall_normal, collision_loc] = map.hit_wall([-1.5,3], [-2.5,3.1]);
disp(hitsWall)
disp(collision_loc)
scatter(collision_loc(1), collision_loc(2))
[hitsWall, wall_normal, collision_loc] = map.hit_wall([-0.9, 2.1], [-1,1.9]);
disp(hitsWall)
disp(collision_loc)
scatter(collision_loc(1), collision_loc(2))
[hitsWall, wall_normal, collision_loc] = map.hit_wall([-1,3],[-2,3]);
disp(hitsWall)
disp(collision_loc)
scatter(collision_loc(1), collision_loc(2))
[hitsWall, wall_normal, collision_loc] = map.hit_wall([1,1],[-1.5,3]);
disp(hitsWall)
disp(collision_loc)
scatter(collision_loc(1), collision_loc(2))
