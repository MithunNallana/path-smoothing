function [ X,Y,Z ] = pathsmoothbezier(wayPointsX,wayPointsY,varargin)
% PATHSMOOTHBEZIER generates a smooth path from a set of waypoints
% usage: [ X,Y,Z ] = pathsmoothbezier(wayPointsX,wayPointsY)
% usage: [ X,Y,Z ] = pathsmoothbezier(wayPointsX,wayPointsY,maximumCurvature)
% usage: [ X,Y,Z ] = pathsmoothbezier(wayPointsX,wayPointsY,wayPointsZ,maximumCurvature)
%
% This function computes a smooth path 
% using a cubuc Bezier curve. This path
% satisfies both curvature continuity and
% maximum curvaure requirements. This method
% provides efficient analytical solution.
%
% arguments: (input)
%
%  wayPointsX,wayPointsY - column vectors of length n. These define
%        x and y coordinates of points to develop bezier curve on them.
%        Number of points must be atleast 3 for this method to work.
%
%  wayPointsZ - column vector of length n. This vector define z coordinates
%        of way points to develop 3d bezier curve.
%
%  maximumCurvature - number that defines maximum curvature allowed on 
%        path. This is especially useful in developing a Global path for 
%        cars which have limits on turning radius
%        DEFAULT VALUE: 3.
%
% arguments: (output)
%
%  X - x-coordinates of waypoints along the smooth path generated
%
%  Y - y-coordinates of waypoints along the smooth path generated
%
%  Y - y-coordinates of waypoints along the smooth path generated.
%      DEFAULT OUTPUT (FOR 2D CURVE): zero vector of length m
%
%  FOR EXAMPLE USAGE REFER "example.m"
%
%  See also: bezier, bernstein
%
%  References:
%   [1] Yang, Kwangjin, and Salah Sukkarieh. "An analytical continuous-
%   curvature path-smoothing algorithm." IEEE Transactions on Robotics 
%   26.3 (2010): 561-568.
%   [2] Walton, D. J., D. S. Meek, and J. M. Ali. "Planar G2 transition 
%   curves composed of cubic Bézier spiral segments." Journal of 
%   Computational and Applied Mathematics 157.2 (2003): 453-476.
%
% -------------------------------------------------------------
% Authors : Mithun Nallana, Balaji C.A.
% Email   : mithunbabu1141995@gmail.com
% Release : 1.0
%
% HISTORY:
%  15/06/2017 - Base version for 2d curves [Private].
%  10/08/2017 - Error for collinear points fixed [Private].
%  20/12/2017 - Error in angle corrected. [Private]
%  12/02/2018 - Implemented for 3d curves and input checks [Private].
%  18/04/2018 - Initial release of code [Public].



if nargin < 2
    error ('PATHSMOOTHBEZIER:insufficientarguments', ...
        'wayPointsX and wayPointsY must be provided')
end

num = size(wayPointsX,2);
if ~isvector(wayPointsX) || ~isvector(wayPointsY) || (length(wayPointsY) ~= num)
    error('PATHSMOOTHBEZIER:improperwaypoints',...
        'wayPointsX and wayPointsY must be vectors of same size')
elseif num < 3
    error('PATHSMOOTHBEZIER:improperwaypoints', ...
        'wayPointsX and wayPointsY must have atlease three points')
end

maximumCurvature = 3;
wayPointsZ       = zeros(1,num);

if (numel(varargin) > 0) && (numel(varargin) < 3) 
    if numel(varargin) == 1
        maximumCurvature = varargin{1};
    else
        arg = varargin{1};
        if ~isvector(arg) || (length(arg) ~= num)
            error('PATHSMOOTHBEZIER:improperwaypoints', ...
                'wayPointsZ is not a vector of same size of wayPointsX')
        end
        wayPointsZ = arg;
        maximumCurvature = varargin{2};
    end
elseif numel(varargin) > 2
    error('PATHSMOOTHBEZIER:toomanyarguments', ...
        'Number of arguments are more that allowed by this function')
end

% parameters required for constuction of bezier curve
% derivation can be provided in [2].
c1 = 7.2364;
c2 = (2*(sqrt(6)-1))/5;
c3 = (c2+4)/(c1+6);

% counter for size of outputvectors X,Y and Z

X(1) = wayPointsX(1);
Y(1) = wayPointsY(1);
Z(1) = wayPointsZ(1);

dt    = 0.01;
tBase = 0:dt:1;
tParm = [nchoosek(3,0)*((1-tBase).^3); nchoosek(3,1)*(tBase.^1).*((1-tBase).^2); ...
    nchoosek(3,2)*(tBase.^2).*((1-tBase).^1); nchoosek(3,3)*(tBase.^3)];
tParmRev = flipud(tParm);
tParmRev = tParmRev(:,2:end);

% selecting three waypoints at a time
for i=1:num-2
    % Waypoint one (W1)
    w1 = [wayPointsX(1,i); wayPointsY(1,i); wayPointsZ(1,i)];
    % Waypoint two (W2)
    w2 = [wayPointsX(1,i+1); wayPointsY(1,i+1); wayPointsZ(1,i+1)];
    % Waypoint three (W3)
    w3 = [wayPointsX(1,i+2); wayPointsY(1,i+2); wayPointsZ(1,i+2)];
    
    % ut,un,ub are three perpendicular unit vectors 
    % at w1 where ut(tangent) is pointing towards w2 and
    % un(normal) is in the direction of w3
    ut = (w2-w1)/norm(w2-w1);
    up = (w2-w3)/norm(w2-w3);
    ub = cross(up,ut);
    ub = ub/norm(ub);
    un = cross(ub,ut);
    
    % build transformation matrix to transform coordinates 
    % from global to local and local to global
    rotationMat     = [ut,un,ub];
    positionMat     = w1;
    transformMat    = [[rotationMat;zeros(1,3)],[positionMat;1]];
    transformMatInv = [[rotationMat';zeros(1,3)],[-rotationMat'*positionMat;1]];
    
    
    m1 = transformMatInv*[w1;1];
    m2 = transformMatInv*[w2;1];
    m3 = transformMatInv*[w3;1];
    
    m1 = m1(1:end-1,:);
    m2 = m2(1:end-1,:);
    m3 = m3(1:end-1,:);
    
    m1m2 = m2-m1;
    m2m3 = m3-m2;
    
    m1m2 = m1m2/norm(m1m2);
    m2m3 = m2m3/norm(m2m3);
    u1 = -m1m2;
    u2 = m2m3;
    
    gamma = acos(sum(m1m2 .* m2m3));
    beta = gamma/2;
    if(beta == 0)
        beta = 0.00001;
    end
    d = ((c2+4)^2/(54*c3))*(sin(beta)/(maximumCurvature*(cos(beta))^2));
    
    hb = c3*d;
    he = hb;
    gb = c2*c3*d;
    ge = gb;
    kb = ((6*c3*cos(beta))/(c2+4))*d;
    ke = kb;
    
    B0 = m2 + d*u1;
    B1 = B0 - gb*u1;
    B2 = B1 - hb*u1;
    
    E0 = m2 + d*u2;
    E1 = E0 - ge*u2;
    E2 = E1 - he*u2;
    
    B2E2 = E2-B2;
    ud = B2E2/norm(B2E2);
    
    B3 = B2 + kb*ud;
    E3 = E2 - ke*ud;
    
    Bmatrix = transformMat*[[B0;1],[B1;1],[B2;1],[B3;1]];
    Ematrix = transformMat*[[E0;1],[E1;1],[E2;1],[E3;1]];
    
    Bmatrix   = Bmatrix(1:3,:);
    
    bezierB = Bmatrix*tParm;
    X = [X,bezierB(1,:)];
    Y = [Y,bezierB(2,:)];
    Z = [Z,bezierB(3,:)];

    Ematrix   = Ematrix(1:3,:);
    
    bezierE = Ematrix*tParmRev;
    X = [X,bezierE(1,:)];
    Y = [Y,bezierE(2,:)];
    Z = [Z,bezierE(3,:)];
end
X = [X,wayPointsX(1,end)];
Y = [Y,wayPointsY(1,end)];
Z = [Z,wayPointsZ(1,end)];
end

