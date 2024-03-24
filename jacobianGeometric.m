%==============================================================================
% Author: Carl Larsson
% Description: Computes the Geometrical Jacobian
% Date: 17-02-2024
%==============================================================================
% Inputs:
% HTto0: cell array of (4x4) Homogeneous Transformation matrices. The
% i-th element of HTto0 converts points from i-th frame to the base
% frame, where i = 1 .. n (symbolic in joint variables).
% kindJoints: cell array of chars specifying if i-th joint is revolute
% ('r') or prismatic ('p').
%==============================================================================
% Outputs:
% Jg: array (matrix), the Geometric Jacobian (6 x n) for the n-th frame
%==============================================================================
function [Jg] = jacobianGeometric(HTto0, kindJoints)

% Number of joints
n = size(kindJoints, 2);
% z in its own frame (whichever that may be)
z = [0;0;1];

% Jacobian for position
J_p = sym(zeros(3,n));
for i=1:n
    % z_0^0 can't be calculated since we don't have A_0^0 in HTto0, however
    % it's just [0;0;1]
    % p_0^0 can't be calculated since we don't have A_0^0 in HTto0, however
    % it's just [0;0;0]
    if(i == 1)
        % Revolute
        if(kindJoints{i} == 'r')
            % J_p(i) = z_{i-1}^0 x (p_e - p_{i-1})
            J_p(:,i) = simplify(cross(z, (HTto0{end}(1:3, end) - [0;0;0])));
        % Prismatic
        else
            % J_p(i) = z_{i-1}^0
            J_p(:,i) = z;
        end
    else
        % Revolute
        if(kindJoints{i} == 'r')
            % J_p(i) = z_{i-1}^0 x (p_e - p_{i-1})
            J_p(:,i) = simplify(cross(HTto0{i-1}(1:3,1:3)*z, (HTto0{end}(1:3, end) - HTto0{i-1}(1:3, end))));
        % Prismatic
        else
            % J_p(i) = z_{i-1}^0
            J_p(:,i) = simplify(HTto0{i-1}(1:3,1:3)*z);
        end
    end
end

% Jacobian for orientation
J_o = sym(zeros(3,n));
for i=1:n
    % Revolute first
    if((i == 1) && (kindJoints{i} == 'r'))
        % z_0^0 can't be calculated since we don't have A_0^0 in HTto0, however
        % it's just [0;0;1]
        % J_o(i) = z_{i-1}^0 = R_{i-1}^0 * z^{i-1}
        J_o(:,i) = z;
    % Revolute
    elseif(kindJoints{i} == 'r')
        % J_o(i) = z_{i-1}^0 = R_{i-1}^0 * z^{i-1}
        J_o(:,i) = simplify(HTto0{i-1}(1:3,1:3)*z);
    % Could make another for prismatic just for clarity, 
    % but that is just 0, and since array
    % is initialized to all 0s, it wouldn't do anything
    end
end


Jg = [J_p;J_o];
end