
function [pose_0,map_coord] = maps(name)
    if strcmp('square',name)
        map_coord = [0 0;
                     0 10;
                     10 10;
                     10 0];
        pose_0 = [5;5;pi/4; 2;0;0];

    elseif strcmp('diamond',name)
        map_coord = [   0,  0;
                       3,  2;
                       5, 5;
                        3, 4];
        
        pose_0 = [3.5;3.5;0;...
                  2;0;0.15];
          
    elseif strcmp('asym',name)
        map_coord = [0     0;
                     2      2;
                     5      0;
                     3      6;
                     2      4;
                     0      3;
                     ];
        pose_0 = [2.5;2.5;0;...
                  2;0;0.15];
    elseif  strcmp('v',name)
        map_coord = [0 0;
                     -3 2;
                     0 -3;
                     3 2;
                     ];
        pose_0 = [0;-1;-pi/4;...
                 -1;-1;0];
    else

    end
end
