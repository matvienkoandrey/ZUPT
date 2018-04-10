function [ins_e] = ins(imu, att_mode, precision)
% ins_gps: loosely-coupled integrated navigation system.
%
% ins_gps integrates IMU and GPS measurements by using an Extended Kalman filter.
%
% INPUT:
%   imu, IMU data structure.
%         t: Ix1 time vector (seconds).
%        fb: Ix3 accelerations vector in body frame XYZ (m/s^2).
%        wb: Ix3 turn rates vector in body frame XYZ (radians/s).
%       arw: 1x3 angle random walks (rad/s/root-Hz).
%       vrw: 1x3 angle velocity walks (m/s^2/root-Hz).
%      gstd: 1x3 gyros standard deviations (radians/s).
%      astd: 1x3 accrs standard deviations (m/s^2).
%    gb_fix: 1x3 gyros static biases or turn-on biases (radians/s).
%    ab_fix: 1x3 accrs static biases or turn-on biases (m/s^2).
%  gb_drift: 1x3 gyros dynamic biases or bias instabilities (radians/s).
%  ab_drift: 1x3 accrs dynamic biases or bias instabilities (m/s^2).
%   gb_corr: 1x3 gyros correlation times (seconds).
%   ab_corr: 1x3 accrs correlation times (seconds).
%     gpsd : 1x3 gyros dynamic biases PSD (rad/s/root-Hz).
%     apsd : 1x3 accrs dynamic biases PSD (m/s^2/root-Hz);
%      freq: 1x1 sampling frequency (Hz).
% ini_align: 1x3 initial attitude at t(1).
% ini_align_err: 1x3 initial attitude errors at t(1).
%
%	gps, GPS data structure.
%         t: Mx1 time vector (seconds).
%       lat: Mx1 latitude (radians).
%       lon: Mx1 longitude (radians).
%         h: Mx1 altitude (m).
%       vel: Mx3 NED velocities (m/s).
%       std: 1x3 position standard deviations (rad, rad, m).
%      stdm: 1x3 position standard deviations (m, m, m).
%      stdv: 1x3 velocity standard deviations (m/s).
%      larm: 3x1 lever arm from IMU to GNSS antenna (x-fwd, y-right, z-down) (m).
%      freq: 1x1 sampling frequency (Hz).
%
%	att_mode, attitude mode string.
%      'quaternion': attitude updated in quaternion format. Default value.
%             'dcm': attitude updated in Direct Cosine Matrix format.
%
%   precision, finite number precision string.
%      double: double float precision (64 bits). Default value.
%      single: single float precision (32 bits).
%
% OUTPUT:
%   ins_e, INS/GPS estimates data structure.
%         t: Ix1 time vector (seconds).
%        tk: Mx1 time when Kalman filter is executed
%      roll: Ix1 roll (radians).
%     pitch: Ix1 pitch (radians).
%       yaw: Ix1 yaw (radians).
%       vel: Ix3 NED velocities (m/s).
%       lat: Ix1 latitude (radians).
%       lon: Ix1 longitude (radians).
%         h: Ix1 altitude (m).
%        Pp: Mx441 Kalman filter a posteriori covariance matrices.
%        Pi: Mx441 Kalman filter a priori covariance matrices.
%         A: Mx441 Kalman filter transition-state matrices.
%         B: Mx12 Kalman filter biases compensations.
%        In: Mx6  Kalman filter innovations.
%        Xi: Mx21 Kalman filter a priori states.
%        Xp: Mx21 Kalman filter a posteriori states.
%
%   Copyright (C) 2014, Rodrigo Gonzalez, all rights reserved.
%
%   This file is part of NaveGo, an open-source MATLAB toolbox for
%   simulation of integrated navigation systems.
%
%   NaveGo is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License (LGPL)
%   version 3 as published by the Free Software Foundation.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this program. If not, see
%   <http://www.gnu.org/licenses/>.
%
% Reference:
%   R. Gonzalez, J. Giribet, and H. Patiño. NaveGo: a
% simulation framework for low-cost integrated navigation systems,
% Journal of Control Engineering and Applied Informatics, vol. 17,
% issue 2, pp. 110-120, 2015. Alg. 2.
%
% Version: 003
% Date:    2017/05/10
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego

if nargin < 3, att_mode  = 'quaternion'; end
if nargin < 4, precision = 'double'; end

%% ZUPT detection algorithm
ZUPT = 'ON';
ZUPT_H = 'OFF';
% ZUPT_THRELHOLD = 1;%0.5;   % m/s
% ZUPT_DELAY = 0; % seconds
% ZUPT_WINDOW = 0.5;        % seconds
zupt = false;       
%%
Mi = (max(size(imu.t)));
global step;
step = ones(1, Mi);

if strcmp(precision, 'single')  % single precision
    
    ti = single(imu.t);
    tg = single(gps.t);
    
    % Preallocate memory for estimates
    roll_e  = single(zeros (Mi, 1));
    pitch_e = single(zeros (Mi, 1));
    yaw_e   = single(zeros (Mi, 1));
    vel_e   = single(zeros (Mi, 3));
    h_e     = single(zeros (Mi, 1));
    
    % Constant matrices
    I = single(eye(3));
    Z = single(zeros(3));
    
    % Kalman matrices for later analysis
    In = single(zeros(Mg, 6));       % Kalman filter innovations
    Pi = single(zeros(Mg, 441));     % Elements from a priori covariance matrix, Pi
    Pp = single(zeros(Mg, 441));     % Elements from a posteriori covariance matrix, Pp
    A  = single(zeros(Mg, 441));     % Elements from transition-state matrix, A
    Xi = single(zeros(Mg, 21));      % Evolution of Kalman filter a priori states, xi
    Xp = single(zeros(Mg, 21));      % Evolution of Kalman filter a posteriori states, xp
    B  = single(zeros(Mg, 12));      % Biases compensantions after Kalman filter correction
    x  = single([ zeros(1,9), imu.gb_fix, imu.ab_fix, imu.gb_drift, imu.ab_drift ]');  % Kalman filter error vector state
    
    % Initialize biases variables
    gb_drift = single(imu.gb_drift');
    ab_drift = single(imu.ab_drift');
    gb_fix   = single(imu.gb_fix');
    ab_fix   = single(imu.ab_fix');    
    
    % Initialize estimates at tti=1
    roll_e (1) = single(imu.ini_align(1));
    pitch_e(1) = single(imu.ini_align(2));
    yaw_e(1)   = single(imu.ini_align(3));
    vel_e(1,:) = single(gps.vel(1,:));
    h_e(1)     = single(gps.h(1));
        
else % double precision
    
    ti = (imu.t);
    
    % Preallocate memory for estimates
    roll_e  = zeros (Mi, 1);
    pitch_e = zeros (Mi, 1);
    yaw_e   = zeros (Mi, 1);
%     yawm_e  = zeros (Mi, 1);
    vel_e   = zeros (Mi, 3);
    h_e     = zeros (Mi, 1);
    
    % Constant matrices
    I = eye(3);
    Z = zeros(3);
    
    % Kalman matrices for later analysis
    In = zeros(Mi, 6);         % Kalman filter innovations
    Pi = zeros(Mi, 441);       % Elements from a priori covariance matrices, Pi
    Pp = zeros(Mi, 441);       % Elements from a posteriori covariance matrices, Pp
    A  = zeros(Mi, 441);       % Elements from transition-state matrices, A
    Xi = zeros(Mi, 21);        % Evolution of Kalman filter a priori states, xi
    Xp = zeros(Mi, 21);        % Evolution of Kalman filter a posteriori states, xp
    B  = zeros(Mi, 12);        % Biases compensantions after Kalman filter correction
    x  = [ zeros(1,9), imu.gb_fix, imu.ab_fix, imu.gb_drift, imu.ab_drift ]';  % Kalman filter error vector state
    
    % Initialize biases variables
    gb_drift = imu.gb_drift';
    ab_drift = imu.ab_drift';
    gb_fix   = imu.gb_fix';
    ab_fix   = imu.ab_fix';
    
    % Initialize estimates at tti = 1
    roll_e(1)  = imu.ini_align(1);
    pitch_e(1) = imu.ini_align(2);
    yaw_e(1)   = imu.ini_align(3);
%     yawm_e(1)  = imu.ini_align(3);
    vel_e(1,:) = imu.ini_vel;
    h_e(1)     = imu.ini_pos(3);
end


% Lat and lon cannot be set in single precision. They need full (double) precision.
lat_e    = zeros (Mi,1);
lon_e    = zeros (Mi,1);
lat_e(1) = double(imu.ini_pos(1));
lon_e(1) = double(imu.ini_pos(2));




DCMnb = euler2dcm([roll_e(1); pitch_e(1); yaw_e(1);]);
DCMbn = DCMnb';
qua   = euler2qua([roll_e(1) pitch_e(1) yaw_e(1)]);

% Initialize Kalman filter matrices

std = [7.8326e-08 1.5639e-07 1];
stdm = [0.5 0.5 1];
stdv = [0.0514 0.0514 0.0514]; 

S.R  = diag([stdv, stdm].^2);
S.Q  = diag([imu.arw, imu.vrw, imu.gb_psd, imu.ab_psd].^2);
S.Pp = diag([imu.ini_align_err, stdv, std, imu.gb_fix, imu.ab_fix, imu.gb_drift, imu.ab_drift].^2);

% UD filter matrices
% [Up, Dp] = myUD(S.P);
% dp = diag(Dp);

% DEC = 0.5 * 180/pi;             % Magnetic declination (deg)

% Initialize matrices for INS/GPS performance analysis
Pp(1,:) = reshape(S.Pp, 1, 441);
B(1,:)  = [gb_fix', ab_fix', gb_drift', ab_drift'];

i = 1;
global zupt_time;
global last_zupt;
zupt_time = zeros(1,Mi);
last_zupt = 0;
% GPS clock is the master clock
for i = 1:(Mi - 1)
    
    %% INERTIAL NAVIGATION SYSTEM (INS)

    % Print a dot on console every 10,000 INS executions
    if (mod(i,10000) == 0), fprintf('. ');  end
    % Print a return on console every 200,000 INS executions
    if (mod(i,200000) == 0), fprintf('\n'); end

    % Index for INS navigation update
    i = i + 1;

    % INS period
    dti = ti(i) - ti(i-1);

    % Correct inertial sensors
    wb_corrected = (imu.wb(i,:)' + gb_fix + gb_drift );
    fb_corrected = (imu.fb(i,:)' + ab_fix + ab_drift );

    % Turn-rates update
    omega_ie_n = earthrate(lat_e(i-1), precision);
    omega_en_n = transportrate(lat_e(i-1), vel_e(i-1,1), vel_e(i-1,2), h_e(i-1));

    % Attitude update
    [qua_n, DCMbn_n, euler] = att_update(wb_corrected, DCMbn, qua, ...
        omega_ie_n, omega_en_n, dti, att_mode);
    roll_e(i) = euler(1);
    pitch_e(i)= euler(2);
    yaw_e(i)  = euler(3);
    DCMbn = DCMbn_n;
    qua = qua_n;

    % Gravity update
    gn = gravity(lat_e(i-1), h_e(i-1));

    % Velocity update
    fn = (DCMbn_n * fb_corrected);
    vel_n = vel_update(fn, vel_e(i-1,:), omega_ie_n, omega_en_n, gn', dti);
    vel_e (i,:) = vel_n;
%         virtual_vel = DCMbn_n * [1 0.1 0.1]'.* (fb_corrected - DCMbn_n' * gn') * dti; % 
%         vel_e (i,:) = vel_e (i-1,:) + virtual_vel';

    % Position update
    pos = pos_update([lat_e(i-1) lon_e(i-1) double(h_e(i-1))], double(vel_e(i,:)), double(dti) );
    lat_e(i) = pos(1);
    lon_e(i) = pos(2);
    h_e(i)   = pos(3);

    % Magnetic heading update
%         yawm_e(i) = hd_update (imu.mb(i,:), roll_e(i),  pitch_e(i), D);

    % ZUPT detection algorithm
    zupt = is_stance_phase( i, dti, imu.fb, imu.wb, 'real_data' );
    
    %% INNOVATIONS
    
    [RM,RN] = radius(lat_e(i), precision);
    Tpr = diag([(RM + h_e(i)), (RN + h_e(i)) * cos(lat_e(i)), -1]);  % radians-to-meters
    
    % Innovations
%     zp = Tpr * ([lat_e(i); lon_e(i); h_e(i);] - [gps.lat(j); gps.lon(j); gps.h(j);]) ...
%         + (DCMbn_n * gps.larm);
    
     zv = (vel_e(i,:))';
    
    
    %% KALMAN FILTER
    
    % Vector to update matrix F
    upd = [vel_e(i,:) lat_e(i) h_e(i) fn'];
    
    % Update matrices F and G
    [S.F, S.G] = F_update(upd, DCMbn_n, imu, dti);
    
    if strcmp(ZUPT, 'ON')
        %if zupt == true
        if step(i) == 1
            if strcmp(ZUPT_H, 'ON')
                H_h = diag([0 0 -1]);
                S.H = [ Z I Z Z Z Z Z;
                Z Z H_h Z Z Z Z;];
                zupt_stdv = [0.01 0.01 0.01];
                zupt_stdh = [1 1 0.1];
                zv = vel_e(i,:)';
                zp = [0 0 h_e(i)-100]';
                S.R = diag([zupt_stdv zupt_stdh]).^2;        
                z = [ zv' zp' ]';
            else
                S.H = [ Z I Z   Z Z Z Z];
%                 zupt_stdv = [0.01 0.01 0.01];
%                 S.R = diag(zupt_stdv).^2;    
                z = zv;
            end
        else
            S.H = [ Z Z Z Z Z Z Z];
            S.R = diag(ones(3,1)).^2;        
            z = zeros(1,3)'; 
        end
    else
        S.H = [ Z Z Z Z Z Z Z];
        S.R = diag(ones(3,1)).^2;        
        z = zeros(1,3)';     
    end       
    
    % Execute the extended Kalman filter
    S = kalman(x, z, S, dti);
    x(10:21) = S.xp(10:21);
    
    %% INS CORRECTIONS
    
    % Quaternion corrections
    % Crassidis. Eq. 7.34 and A.174a.
    antm = [0 qua_n(3) -qua_n(2); -qua_n(3) 0 qua_n(1); qua_n(2) -qua_n(1) 0];
    qua = qua_n + 0.5 .* [qua_n(4)*eye(3) + antm; -1.*[qua_n(1) qua_n(2) qua_n(3)]] * S.xp(1:3);
    qua = qua / norm(qua);       % Brute-force normalization
    
    % DCM correction
    DCMbn = qua2dcm(qua);
    %     E = skewm(S.xp(1:3));
    %     DCMbn = (eye(3) + E) * DCMbn_n;
    
    % Attitude corrections
    %     euler = qua2euler(qua);
    %     roll_e(i) = euler(1);
    %     pitch_e(i)= euler(2);
    %     yaw_e(i)  = euler(3);
    
    % Attitude corrections
    roll_e(i)  = roll_e(i)  - S.xp(1);
    pitch_e(i) = pitch_e(i) - S.xp(2);
    yaw_e(i)   = yaw_e(i)   - S.xp(3);
    
    % Velocity corrections
    vel_e (i,1) = vel_e (i,1) - S.xp(4);
    vel_e (i,2) = vel_e (i,2) - S.xp(5);
    vel_e (i,3) = vel_e (i,3) - S.xp(6);
    
    % Position corrections
    lat_e(i) = lat_e(i) - double(S.xp(7));
    lon_e(i) = lon_e(i) - double(S.xp(8));
    h_e(i)   = h_e(i)   - S.xp(9);
    
    % Biases corrections
    gb_fix   = S.xp(10:12);
    ab_fix   = S.xp(13:15);
    gb_drift = S.xp(16:18);
    ab_drift = S.xp(19:21);
    
    % Matrices for later INS performance analysis
    Xi(i,:) = S.xi';
    Xp(i,:) = S.xp';
    Pi(i,:) = reshape(S.Pi, 1, 441);
    Pp(i,:) = reshape(S.Pp, 1, 441);
    A(i,:)  = reshape(S.A, 1, 441);
     
%     if(zupt == true)
%         
%         In(j,:) = [ zv; zeros(3,1);]';
%     else
%         In(j,:) = z';   
%     end
    B(i,:)  = [gb_fix', ab_fix', gb_drift', ab_drift'];
    
    zupt = false;
    
end


%% Estimates from INS/GPS integration

ins_e.t     = ti(1:i, :);       % IMU time
ins_e.tk    = ti;               % Time when Kalman filter is executed
ins_e.roll  = roll_e(1:i, :);   % Roll
ins_e.pitch = pitch_e(1:i, :);  % Pitch
ins_e.yaw   = yaw_e(1:i, :);    % Yaw
% ins_e.yawm  = yawm_e(1:i, :);    % Magnetic heading
ins_e.vel   = vel_e(1:i, :);    % NED velocities
ins_e.lat   = lat_e(1:i, :);    % Latitude
ins_e.lon   = lon_e(1:i, :);    % Longitude
ins_e.h     = h_e(1:i, :);      % Altitude
ins_e.Pi    = Pi;       % Pi matrices
ins_e.Pp    = Pp;       % Pp matrices
ins_e.A     = A;        % A matrices
ins_e.B     = B;        % Kalman filter biases compensations
ins_e.In    = In;       % Kalman filter innovations
ins_e.Xi    = Xi;       % Kalman filter a priori states
ins_e.Xp    = Xp;       % Kalman filter a posteriori states

fprintf('\n');

end