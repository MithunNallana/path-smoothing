% Example 1

clear;
clc;

% 2d waypoints to build path
fx        = [0,10,10,20,20,30,30,40,40,-10,-10,  0,0];
fy        = [0, 0,10, 0,10,10, 0, 0,20, 20,-20,-20,0];
curvature = 1;

% calling function for 2d points
[X,Y,Z] = pathsmoothbezier(fx,fy,curvature);
%[X,Y,Z] = pathsmoothbezier(fx,fy);

% plotting points
figure(2);
hold on;
plot3(X,Y,Z)
axis equal;
plot3(fx,fy,fz,'rO');

% % % %% Example 2
% % % 
% % % clear;
% % % clc;
% % % 
% % % % 2d waypoints to build path
% % % fx = [0,10,10, 0,0];
% % % fy = [0, 0,10,10,0];
% % % fz = [0, 0,10,10,0];
% % % 
% % % % calling function for 3d points
% % % [X,Y,Z] = pathsmoothbezier(fx,fy,fz,1);
% % % 
% % % % plotting points
% % % hold on;
% % % plot3(X,Y,Z)
% % % axis equal;
% % % plot3(fx,fy,fz,'rO');