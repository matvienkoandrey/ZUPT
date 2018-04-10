function ret_val = is_stance_phase( sample_idx, dt, acc_b, gyro_b, data_mode )
    global step;
    if strcmp(data_mode, 'model_data')
        ZUPT_THRESHOLD = 1;%0.5;  % m/s
        ZUPT_DELAY = 0;           % seconds
        ZUPT_WINDOW = 0.5;        % seconds

        global zupt_time;
        global last_zupt;

        ret_val = false;

        % ZUPT detection algorithm
        idz = floor( ZUPT_WINDOW / dt ); % Index to set ZUPT window time
        if ( sample_idx > idz )  
            acc_m = std( acc_b( sample_idx - idz : sample_idx, 1 ));
            if ( abs( acc_m ) <= ZUPT_THRESHOLD && ... 
                ( max( acc_b( sample_idx - idz : sample_idx, 1 )) - min( acc_b( sample_idx - idz : sample_idx , 1 ))) <= ( ZUPT_THRESHOLD - 0.8 ))       
                if abs( last_zupt - sample_idx ) * dt > ZUPT_DELAY
                    %disp('zupt')    % For debugging purposes
                    zupt_time( sample_idx ) = 1;
                    ret_val = true;
                    last_zupt = sample_idx;
                else
                    ret_val = false;
                end
            end
        end
    else
        d_size = length(acc_b);
        w_size = 15;
        acc_wndw = zeros(1, w_size);
        gyro_wndw = zeros(1, w_size);
        acc_var_threshold = 60;
        gyro_var_threshold = 10;
        acc_threshold = 2;
        gyro_threshold = 0.6;
        g_norm = 9;
        %time_delay = 0.05;
        %time(1) = 0;
        step_delay = 6;
        if(sample_idx < (d_size - w_size))
            %time(i + 1) = (TIME_StartTime(i + 1) - TIME_StartTime(i))/1000000 + time(i);
            acc_wndw = sqrt(acc_b(sample_idx:sample_idx+w_size,1).^2+acc_b(sample_idx:sample_idx+w_size,2).^2+acc_b(sample_idx:sample_idx+w_size,3).^2);
            gyro_wndw = sqrt(gyro_b(sample_idx:sample_idx+w_size,1).^2+gyro_b(sample_idx:sample_idx+w_size,2).^2+gyro_b(sample_idx:sample_idx+w_size,3).^2);
            acc_var(sample_idx) = var(acc_wndw);
            gyro_var(sample_idx) = var(gyro_wndw);
            %if acc_var(i) < acc_var_threshold && ang_var(i) < ang_var_threshold && ang_wndw(round(w_size/2)) < ang_threshold
            if (gyro_var(1) < gyro_threshold && ...
                var(acc_wndw) < acc_var_threshold && ...
                abs(mean(acc_wndw) - g_norm) < acc_threshold)
                    step(sample_idx:sample_idx+step_delay) = 0;
                    ret_val = true;
            else
                    ret_val = false;
            end
        else
            ret_val = false;
        end
    end
end

