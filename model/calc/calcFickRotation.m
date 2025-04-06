function [yaw, pitch, torsion] = calcFickRotation(R, worldTarget)
% Calculates Fick angles required to direct a ray towards a target
%
% Syntax:
%  [yaw, pitch, torsion] = calcFickRotation(R, worldTarget)
%
% Description:
%   Calculates the 3D rotation required to direct a ray towards a target
%   point. First, the 3x3 rotation matrix is obtained. This is then
%   decomposed into the Fick angles.
%
%   GKA used the Gemini 2.5 LLM to generate this code.
%
% Inputs:
%   R                     - 3x2 matrix that specifies a vector of the form 
%                           [p; u], corresponding to
%                               R = p + t*u
%                           where p is vector origin, u is the direction
%                           expressed as a unit step, and t is unity for a
%                           unit vector.
%   worldTarget           - 3x1 vector that is the point the rotated ray
%                           should intersect.
%
% Outputs:
%   yaw, pitch, torsion   - Scalars. The Fick angles.
%

% --- Constants ---
TOLERANCE = 1e-10; % Tolerance for checking zero norms, alignment etc.
GIMBAL_TOL = 1e-7; % Tolerance for detecting gimbal lock singularity

ray_origin        = R(:,1);
initial_direction = R(:,2);

% --- Validate Initial Direction ---
norm_initial = norm(initial_direction);
if norm_initial < TOLERANCE
    error('Initial direction vector must have non-zero magnitude.');
end
u_initial = initial_direction / norm_initial; % Normalize initial direction

% --- Calculate Target Direction ---
target_direction_vector = worldTarget - ray_origin;
norm_target_dir = norm(target_direction_vector);

if norm_target_dir < TOLERANCE
    pitch = 0; yaw = 0; torsion = 0;
    return;
end
u_target = target_direction_vector / norm_target_dir; % Normalize target direction

% --- Calculate Rotation Axis and Angle ---
dot_prod = dot(u_initial, u_target);
dot_prod = max(-1.0, min(1.0, dot_prod)); % Clamp dot product for acos

if dot_prod > (1.0 - TOLERANCE)
    % Vectors are already aligned, no rotation needed.
    pitch = 0; yaw = 0; torsion = 0;
    return;
elseif dot_prod < (-1.0 + TOLERANCE)
    % Vectors are antiparallel (180 degree rotation).
    angle = pi;
    % Axis is arbitrary and perpendicular to u_initial. Find one robustly.
    axis = find_perpendicular(u_initial);
    if norm(axis) < TOLERANCE
         error('Could not find perpendicular axis for antiparallel vectors.');
    end
    axis = axis / norm(axis); % Ensure unit vector
else
    % General case: Calculate axis and angle using cross product.
    angle = acos(dot_prod);
    axis = cross(u_initial, u_target);
    axis_norm = norm(axis);
    if axis_norm < TOLERANCE
        % This case should ideally be caught by the dot_prod check,
        % but for robustness: if cross product is zero and vectors aren't
        % identical/antiparallel (which they shouldn't be here),
        % something is wrong. Assume identity or error.
        warning('Cross product near zero for non-aligned vectors.');
        pitch = 0; yaw = 0; torsion = 0;
        return
    end
    axis = axis / axis_norm; % Normalize rotation axis
end

% --- Convert Axis-Angle to Rotation Matrix R ---
% Using Rodrigues' formula (can replace with axang2rotm if toolbox available)
Rmat = rodrigues_rotation(axis, angle);


% --- Extract Angles ---
% From R = Rz(yaw) * Ry(pitch) * Rx(torsion)
%       | cy*cp    cy*sp*sr - sy*cr    cy*sp*cr + sy*sr |
% R =   | sy*cp    sy*sp*sr + cy*cr    sy*sp*cr - cy*sr |
%       | -sp           cp*sr              cp*cr       |
% cy=cos(yaw), sy=sin(yaw), cp=cos(pitch), sp=sin(pitch), cr=cos(torsion), sr=sin(torsion)

% Calculate Pitch (beta)
% R(3,1) = -sin(pitch)
sin_pitch = -Rmat(3,1);
% Clamp sin_pitch to [-1, 1] to avoid complex results from asin due to numerical error
sin_pitch = max(-1.0, min(1.0, sin_pitch));
pitch_rad = asin(sin_pitch); % Pitch is in [-pi/2, pi/2]

% Check for Gimbal Lock (pitch is +/- 90 degrees)
cos_pitch = cos(pitch_rad);

if abs(cos_pitch) < GIMBAL_TOL
    % Gimbal Lock: Pitch is +/- pi/2, cos(pitch) is near 0.
    % Yaw and torsion are coupled. Conventionally, set torsion to 0
    % and calculate the combined rotation as yaw (representing yaw +/- torsion).
    torsion_rad = 0;

    if sin_pitch > (1.0 - GIMBAL_TOL) % Pitch is close to +pi/2 (R(3,1) close to -1)
        pitch_rad = pi/2; % Force exact value
        % Calculate yaw = alpha - gamma = atan2(-R(1,2), R(1,3))
        % R(1,2) = cy*sp*sr - sy*cr = cy*sr - sy*cr = -sin(yaw-torsion)
        % R(1,3) = cy*sp*cr + sy*sr = cy*cr + sy*sr = cos(yaw-torsion)
        yaw_rad = atan2(-Rmat(1,2), Rmat(1,3));
    else % Pitch is close to -pi/2 (R(3,1) close to +1)
        pitch_rad = -pi/2; % Force exact value
        % Calculate yaw = alpha + gamma = atan2(-R(1,2), -R(1,3))
        % R(1,2) = -cy*sr - sy*cr = -sin(yaw+torsion)
        % R(1,3) = -cy*cr + sy*sr = -cos(yaw+torsion)
        yaw_rad = atan2(-Rmat(1,2), -Rmat(1,3));
    end
else
    % Normal Case: cos(pitch) is not near zero.
    % Calculate Yaw (alpha)
    % R(1,1) = cos(yaw)*cos(pitch)
    % R(2,1) = sin(yaw)*cos(pitch)
    yaw_rad = atan2(Rmat(2,1), Rmat(1,1)); % atan2(sy*cp, cy*cp)

    % Calculate Torsion (gamma)
    % R(3,2) = cos(pitch)*sin(torsion)
    % R(3,3) = cos(pitch)*cos(torsion)
    torsion_rad = atan2(Rmat(3,2), Rmat(3,3)); % atan2(cp*sr, cp*cr)
end

% --- Convert Radians to Degrees ---
yaw = rad2deg(yaw_rad);
pitch = rad2deg(pitch_rad);
torsion = rad2deg(torsion_rad);


end % End of main function


% --- Helper Functions ---

function R = rodrigues_rotation(axis_in, angle)
    % Rodrigues' rotation formula implementation to get rotation matrix.
    if norm(axis_in) < 1e-10
        % If axis is zero (e.g., angle is zero), return identity
        R = eye(3);
        return;
    end
    axis = axis_in(:) / norm(axis_in(:)); % Ensure unit vector axis
    K = [  0,      -axis(3),  axis(2);
           axis(3),  0,      -axis(1);
          -axis(2), axis(1),  0      ]; % Skew-symmetric matrix
    R = eye(3) + sin(angle) * K + (1 - cos(angle)) * (K*K); % Rodrigues' formula
end

function perp_axis = find_perpendicular(vec)
    % Finds an arbitrary vector perpendicular to input vector vec.
    % Handles cases where vec aligns with world axes.
    vec = vec(:);
    abs_vec = abs(vec);
    TOL = 1e-10; % Tolerance

    % Try crossing with the world axis corresponding to the smallest component of vec
    if abs_vec(1) <= abs_vec(2) + TOL && abs_vec(1) <= abs_vec(3) + TOL
        ref_axis = [1; 0; 0]; % Try X-axis first
    elseif abs_vec(2) <= abs_vec(1) + TOL && abs_vec(2) <= abs_vec(3) + TOL
        ref_axis = [0; 1; 0]; % Try Y-axis first
    else
        ref_axis = [0; 0; 1]; % Try Z-axis first
    end

    perp_axis = cross(ref_axis, vec);

    % If cross product is near zero (vec is parallel to ref_axis), try another axis
    if norm(perp_axis) < TOL
        if abs(ref_axis(1) - 1) < TOL % If we first tried X-axis (and failed)
             perp_axis = cross([0; 1; 0], vec); % Try Y-axis
        else % If we first tried Y or Z axis (and failed)
             perp_axis = cross([1; 0; 0], vec); % Try X-axis
        end
         % If it's still zero (e.g. vec is [0;0;0]), norm check below handles it
    end

     % Final check in case vec was zero vector initially
     if norm(perp_axis) < TOL
         % This should only happen if vec is the zero vector.
         % Return an arbitrary axis like X, although this case
         % should ideally be handled before calling find_perpendicular.
         perp_axis = [1; 0; 0];
     end

end