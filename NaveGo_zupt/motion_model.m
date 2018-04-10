dt = 0.01;
sample_rate = 1 / dt;

t = 0:dt:(100 - dt);
ref.t = t';

% step period = 1 s
step_period = 10;
step_samples_count = int16(step_period * sample_rate);

duty_cycle = 0.6;

swing_duration = step_period * duty_cycle;
swing_samples_count = int16(swing_duration * sample_rate);

stance_duration = step_period - swing_duration;
stance_samples_count = int16(stance_duration * sample_rate);
f_swing = 1 / swing_duration;

T = 0:dt:(swing_duration - dt);

% full time
steps_count = 10;

% swing phase acceleration values (Y-axis)
ay_swing = 10 * sin(2 * pi * f_swing * T);

% stance phase acceleration values (Y-axis)
ay_stance = zeros(1, stance_samples_count);

% all acceleration values (Y-axis)
% for i = 1:20
%     j = i * sample_rate - (sample_rate - 1);
%     if mod(i, 2) == 0
%         ay(1, j:(j + (sample_rate - 1))) = ay_stance;
%     else
%         ay(1, j:(j + (sample_rate - 1))) = ay_swing;
%     end
% end
j = 1;
swing = 1;
while j <= length(t)
    if swing
        ay(1, j:(j + (swing_samples_count - 1))) = ay_swing;
        j = j + swing_samples_count;
        swing = 0;
    else
        ay(1, j:(j + (stance_samples_count - 1))) = ay_stance;
        j = j + stance_samples_count;
        swing = 1;
    end
end
ax = zeros(1, length(t));
az = zeros(1, length(t));

ref.kn = length(t);
ref.freq = 100;

ref.roll(:, 1) = zeros(ref.kn, 1)';
ref.pitch(:, 1) = zeros(ref.kn, 1)';
% ref.yaw(:, 1) = pi / 2 * ones(ref.kn / 2, 1)';
% ref.yaw((ref.kn/2)+1:ref.kn, 1) = 2 * pi / 3 * ones(ref.kn / 2, 1)';
ref.yaw(:, 1) = 0 * pi / 2 * ones(ref.kn, 1)';

for i = 1:ref.kn
    DCM = euler2dcm([ref.roll(i); ref.pitch(i); ref.yaw(i);]);
    a_NED(i, 1:3) = DCM' * [ay(i); ax(i); az(i)];
    ref.DCMnb(i, 1:9) = reshape(DCM, 1, 9);
end

% V_NED(:, 1) = 100 * cos(ref.yaw) .* ones(ref.kn, 1);
% V_NED(:, 2) = 100 * sin(ref.yaw) .* ones(ref.kn, 1);
% V_NED(:, 3) = zeros(ref.kn, 1);

V_NED(1, 1:3) = zeros(1, 3);

ref.lat(1) = 60 * pi / 180;
ref.lon(1) = 40 * pi / 180;
ref.h(1) = 100;

for i = 2:ref.kn
    V_NED(i, 1:3) = V_NED(i-1, 1:3) + a_NED(i-1, 1:3) * dt;
    
    pos = pos_update([ref.lat(i-1) ref.lon(i-1) ref.h(i-1)], V_NED(i, :), dt);
    
    ref.lat(i, 1) = pos(1);
    ref.lon(i, 1) = pos(2);
    ref.h(i, 1) = pos(3);
end
[ref.vel, ~] = pllh2vned(ref);
%ref.vel = V_NED;

save model_data ref