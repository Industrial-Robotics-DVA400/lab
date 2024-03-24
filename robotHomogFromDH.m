%==============================================================================
% Author: Carl Larsson
% Description: Computes the DH homogeneous transformation matrices A_i from
% frame i to frame 0
% Date: 17-02-2024
%==============================================================================
% Inputs:
% dhparams: struct with field-names names 'd', 'theta', 'a', 'alpha'.
% Each field has field-data an array of size 1 x n containing the
% correspondent entry for a Denavit Hartenberg table, where the joint
% variables are symbolic.
%==============================================================================
% Outputs:
% HTto0: cell array of (4x4) Homogeneous Transformation matrices. The
% i-th element of HTto0 converts points from i-th frame to the base
% frame, where i = 1 .. n (symbolic in joint variables).
%==============================================================================
function [HTto0] = robotHomogFromDH(dhparams)

% Number of joints
n = size(dhparams.d, 2);
% HTto0 is a cell array
HTto0_temp = {1, n};

% Compute all the Homogeneous Transformation matrices from frame i to frame
% i-1, A_i^{i-1}
for i=1:n
    HTto0_temp{i} = simplify([cos(dhparams.theta(i)), -cos(dhparams.alpha(i))*sin(dhparams.theta(i)), sin(dhparams.alpha(i))*sin(dhparams.theta(i)) , dhparams.a(i)*cos(dhparams.theta(i));
                              sin(dhparams.theta(i)), cos(dhparams.alpha(i))*cos(dhparams.theta(i)) , -sin(dhparams.alpha(i))*cos(dhparams.theta(i)), dhparams.a(i)*sin(dhparams.theta(i));
                              0                     , sin(dhparams.alpha(i))                        , cos(dhparams.alpha(i))                        , dhparams.d(i)                       ;
                              0                     , 0                                             , 0                                             , 1                                   ]);
end

% Compute all the Homogeneous Transformation matrices from frame i to frame
% 0, A_i^0
for i = n:-1:1
    temp = i-1;
    while(temp > 0)
        HTto0_temp{i} = simplify(HTto0_temp{temp} * HTto0_temp{i});
        temp = temp - 1;
    end
end


HTto0 = HTto0_temp;
end