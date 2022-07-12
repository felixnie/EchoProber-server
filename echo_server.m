% EchoServer
% 08 Jul 2022
% 
% This EchoServer implements:
% 1. Continuous data receiving and plotting.
% 2. Resolving both short message and long data from clients.
% 3. Multiple client connection management.
% 4. Multiple figure background refresh mechanism.
% 5. Server to client messaging.
% 6. Automatic disconnect upon restarting server.
% 
% Notes:
% Please upgrade to MATLAB R2021b or later to use tcpserver.
% 
% To-do:
% 1. Parallel processing for analysis tasks.
% 1. Parallel processing for each analysis tasks.
% 2. Fast and light data caching (smaller cache and avoid copying vars).

clear server
close all

% in light mode, data and data_filtered will not keep growing forever
% light_mode = true;

% set up server
localhost = "155.69.142.178";

global data data_filtered start_list f h server

% ask clients to disconnect
% confusingly this still works after server is cleared
for idx = 1:3
    if ~isempty(server) && ~isempty(server{idx}) && server{idx}.Connected
        try
            writeline(server{idx},"Disconnect.")
            prefix = strcat("[", num2str(idx), "] ");
            disp(strcat(prefix, "Message sent: Disconnect."))
        catch ME
            disp(ME)
        end
    end
end

% create global variables
% support multiple clients, please assign the ports manually
data            = cell(1,3);    % all data received
data_filtered   = cell(1,3);    % all data filtered
start_list      = cell(1,3);    % list of start point index
f               = cell(1,3);    % figure list
h               = cell(1,3);    % plot list
server          = cell(1,3);    % server list
port            = 8170:8172;    % device 1~3 on 8170, 8171 and 8172

% set callback functions
try
    for idx = 1:3
        data{idx}           = [];
        data_filtered{idx}  = [];
        start_list{idx}     = [];

        server{idx} = tcpserver(localhost, port(idx));
        server{idx}.ConnectionChangedFcn = @connectionInfo;
        configureTerminator(server{idx},"CR/LF")
        configureCallback(server{idx}, "terminator", @readData);
    end
catch ME
    switch ME.identifier
        case 'instrument:interface:tcpserver:cannotConnect'
            prefix = strcat("[", num2str(idx), "] ");
            disp(strcat(prefix, "Connection failed. Please try again."));
            return
        otherwise
            rethrow(ME)
    end
end

% callback function for receiving data
function readData(src,~)
    global data data_filtered start_list f h server

    %% get server index
    for idx = 1:3
        if src == server{idx}
            prefix = strcat("[", num2str(idx), "] ");
            break
        end
    end

    %% receive data
    bytes_received = src.NumBytesAvailable;
    
    if  bytes_received == 0 % do nothing
        
        % to-do: check where the empty packets are from
        % disp(strcat(prefix, "Bytes received: ", num2str(bytes_received), " - ", ...
        %                     "Empty packet"))
        flush(src)
        return

    elseif bytes_received < 20 % short message

        msg = native2unicode(readline(src), "UTF-8"); % msg without \r\n
        disp(strcat(prefix, "Bytes received: ", num2str(bytes_received), " - ", ...
                            "Message: ", msg))
        flush(src)
        % reset data when play button is pressed
        if strcmp(msg, "Start playing.")
            clearData(idx)
            disp(strcat(prefix, "Data cleared."))
        end
        return

    else % long data

        src.UserData = read(src, (bytes_received-2)/2, 'int16'); % drop \r\n
        disp(strcat(prefix, "Bytes received: ", num2str(bytes_received), " - ", ...
                            "Data length: ", num2str(length(src.UserData))))
        flush(src)

    end

    %% set up filter and chirp finder
    fc                  = 10000;
    fs                  = 44100;
    [b,a]               = butter(5, fc/(fs/2), 'high'); % same as Python ver?

    % data              : all raw data
    % data_filtered     : all filtered data
    data{idx}           = [data{idx}, src.UserData];
    data_filtered{idx}  = filter(b, a, data{idx});
    % chunk             : new raw data
    % chunk_filtered    : new filtered data
    chunk               = src.UserData;
    chunk_filtered      = data_filtered{idx}(end-length(chunk)+1 : end);  

    % maximum_dalay     : first start point should be within maximum_dalay
    maximum_delay       = 12000;
    % duration          : designed duration between two start points
    duration            = 4000;

    % start_processed   : last processed/displayed chunk index
    start_processed     = length(start_list{idx});

    % allow_silence     : display no matter if it's a chirp
    allow_silence       = true;
    % silence_threshold : the amplitude threshold when there's no chirp
    silence_threshold   = 500; % can be determined through the delay part

    % check if there's enough data for finding start point
    if length(data{idx}) < maximum_delay
        % disp(strcat(prefix, "Data length < ", num2str(maximum_delay)))
        disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - Lack of data, pass."))
        return % do nothing
    end
  
    %% find the first start point
    if isempty(start_list{idx})
        dynamic_threshold = max(data_filtered{idx}) / 5; % find maximum in filtered data
        start_point = 1;
        while data_filtered{idx}(start_point) < dynamic_threshold
            start_point = start_point + 1;
        end
        start_list{idx} = [start_list{idx} start_point]; % get first start point

        % set silence threshold
        silence_threshold = 2 * max(abs(data_filtered{idx}(1 : start_point-10))); % fetch background noise
        disp(strcat(prefix, "Silence threshold = ", num2str(silence_threshold)))
    end

    %% find the rest of the start points
    while length(data_filtered{idx}) >= start_list{idx}(end) + duration
        % set dynamic threshold
        % dynamic_threshold = max(chunk_filtered) / 5; % find maximum in filtered chunk
        % dynamic_threshold = max([dynamic_threshold silence_threshold]); % will not find start when there's no chirp
        % dynamic_threshold = silence_threshold;
        dynamic_threshold = max(abs(data_filtered{idx}(end-maximum_delay+1:end))) / 5; % find maximum in last 12000 points
        disp(strcat(prefix, "Dynamic threshold = ", num2str(dynamic_threshold)))

        % re-search for start point
        start_point = start_list{idx}(end) + duration - 200;
        no_start_point = false;
        while data_filtered{idx}(start_point) < dynamic_threshold
            start_point = start_point + 1;
            if start_point == length(data_filtered{idx})
                no_start_point = true;
                break % make sure start_point will not exceed data length
            end
        end
        if no_start_point % no data larger than threshold till the end
            if allow_silence
                start_point = start_list{idx}(end) + duration; % make up a start point
                start_list{idx} = [start_list{idx} start_point]; % add new start point
            end
        else % there exists data larger than threshold during re-search
            start_list{idx} = [start_list{idx} start_point]; % add new start point
        end
    end
    
    %% print the result for start point finding
    if length(start_list{idx}) > start_processed % new start point added
        if length(start_list{idx}) > 1 % have more than 2 start points, continue
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: ", num2str(start_list{idx}(end) - start_list{idx}(end-1)), " - ", ...
                                "New start point added."))
        else % have only 1 start point, do nothing
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: N/A", " - ", ...
                                "No enough start point."))
            return
        end
    else % no new start point added, do nothing
        if length(start_list{idx}) > 1
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: ", num2str(start_list{idx}(end) - start_list{idx}(end-1)), " - ", ...
                                "No new start point."))
            return
        else % have only 1 start point, do nothing
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: N/A", " - ", ...
                                "Wait for more data."))
            return
        end
    end
    
    %% initialize plot index and start analysis
    plt_idx = 0;

    % current full chunk. this will be used in all code blocks
    chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));

    % original chirp signal
    sig_len = 500;
    original_chirp = zeros(1, sig_len);
    f0 = 15000;
    f1 = 20000;
    w0 = f0 / fs * 2 * pi;
    w1 = f1 / fs * 2 * pi;
    for i = 1:sig_len
        K = sig_len * w0 / log(w1/w0);
        L = sig_len / log(w1/w0);
        phase = K * (exp(i/L) - 1.0);
        original_chirp(i) = sin(phase);
    end
    % 

    % Below are the scripts being called everytime when new data chunk is 
    % arrived. Change the index array to enable/disable each part. For
    % example:
    %     if ismember(idx, [1 2 3]) % enable script for device 1, 2, 3
    %     if ismember(idx, [])      % disable script for all devices
  
    %% 1.1 plot signal (chirp and echo)
    % plot time domain data after low-pass filtering
    % the length of chunk_plot is 4000
    if ismember(idx, [])
        plt_idx = plt_idx + 1; % current plot index
    
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx}) % figure not defined yet or figure being closed
            f{idx}{plt_idx} = figure(idx*100+plt_idx); % create or recreate figure handler
            h{idx}{plt_idx} = plot(chunk_plot); % create or recreate plot handler
            title("1.1 Filtered Signal (1~4000)")
            axis([0 4000 -10000 10000])
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', chunk_plot); % if figure defined and opened, reload data
    end

    %% 1.2 plot signal (echo part)
    % plot time domain data after low-pass filtering
    % the length of chunk_echo is 2420 (2920-550)
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_echo = chunk_plot(550:2920-1);
    
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            h{idx}{plt_idx} = plot(550:2920-1, chunk_echo);
            title("1.2 Filtered Signal (550~2919)")
            axis([550 2920 -2000 2000])
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', chunk_echo);
    end
    
    %% 1.3 cross-correlation
    % re-writing needed, please do not use this
    if ismember(idx, [])
        plt_idx = plt_idx + 1;

        % chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);
    
        % chunk_chirp_win = chunk_chirp .* hanning(length(chunk_chirp))';
        [c, lags] = xcorr(original_chirp, chunk_echo);

        % r = (1:2920-500) + length(chunk_echo); % ?
        % x = lags(r) / fs * 340 / 2;
        % [yupper, ~] = envelope(c(r), 10, 'rms');
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % h{idx}{plt_idx} = plot(x, yupper);
            h{idx}{plt_idx} = plot(lags, c);
            xlabel('lag in distance / m'); ylabel('cross-corr');
            title("1.3 Envelope of Cross-correlation")
            axis([lags(1) lags(end) -max(abs(c)) max(abs(c))])
            grid on
        end
        % set(h{idx}{plt_idx}, 'ydata', yupper);
        set(h{idx}{plt_idx}, 'ydata', c);
    end

    %% 2.1.1 STFT (echo part, 12x48)
    % short-time fourier transform, resolution: 12x48
    if ismember(idx, [1 2 3])
        plt_idx = plt_idx + 1;
        % chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        n_window = 96;
        n_overlap = 48;
        n_fft = 96;
        [ss,ff,tt] = stft(chunk_echo, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
        
        % to capture the frequency interval f0~f1
        f0 = 15000;
        f1 = 20000;
        segment = (fs/2) / n_overlap;
        lower_f = int8(f0/segment) - 1;
        upper_f = int8(f1/segment);
        ff = ff(lower_f:upper_f-1);
        ss = ss(lower_f:upper_f,:);
        
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            h{idx}{plt_idx} = imagesc(tt * 170, ff / 1000, I); % distance as x-axis
            xlabel('distance / m'); ylabel('frequency / kHz');
    
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff / 1000, I); % time as x-axis
            % xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("2.1.1 STFT of Filtered Signal (550~2919) 12x48")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end
    
    %% 2.1.2 STFT (echo part, higher resolution)
    % short-time fourier transform, resolution will be shown in the title
    if ismember(idx, [1 2 3])
        plt_idx = plt_idx + 1;
        % chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        n_window = 256;
        n_overlap = 63/64 * n_window;
        n_fft = n_window * 8;

        [ss,ff,tt] = stft(chunk_echo, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');

        ss = ss((ff>14000 & ff<21000),:); % frequency range you're interest in

        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            h{idx}{plt_idx} = imagesc((tt+550/44100) * 170, ff(ff>f0 & ff<f1) / 1000, I); % distance as x-axis
            xlabel('distance / m'); ylabel('frequency / kHz');
    
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I); % time as x-axis
            % xlabel('time / ms'); ylabel('frequency / kHz');
            
            resolution = [num2str(size(ss,1)), 'x', num2str(size(ss,2))];
            title(["2.1.2 STFT of Filtered Signal (550~2919)" resolution])
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 2.1.3 STFT (chirp part, higher resolution)
    % short-time fourier transform, resolution will be shown in the title
    if ismember(idx, [1 2 3])
        plt_idx = plt_idx + 1;
        chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        n_window = 256;
        n_overlap = 63/64 * n_window;
        n_fft = n_window * 8;

        [ss,ff,tt] = stft(chunk_chirp, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');

        ss = ss((ff>14000 & ff<21000),:); % frequency range you're interest in

        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            % h{idx}{plt_idx} = imagesc(tt * 170, ff(ff>f0 & ff<f1) / 1000, I); % distance as x-axis
            % xlabel('distance / m'); ylabel('frequency / kHz');
    
            h{idx}{plt_idx} = imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I); % time as x-axis
            xlabel('time / ms'); ylabel('frequency / kHz');
            
            resolution = [num2str(size(ss,1)), 'x', num2str(size(ss,2))];
            title(["2.1.3 STFT of Filtered Signal (1~500)" resolution])
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 2.2 spectrogram (echo part, 12x48)
    % spectrogram includes STFT, but with more options, resolution: 12x48
    % this is just to validate that the result is the same as 2.1
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        % chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        n_window = 96;
        n_overlap = 48;
        n_fft = 96;
        [ss,ff,tt] = spectrogram(chunk_echo, tukeywin(n_window,0.25), n_overlap, n_fft, fs, 'yaxis');
        f0 = 15000;
        f1 = 20000;
        segment = (fs/2) / n_overlap;
        lower_f = int8(f0/segment) - 1;
        upper_f = int8(f1/segment);
        ff = ff(lower_f:upper_f-1);
        ss = ss(lower_f:upper_f,:);
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            h{idx}{plt_idx} = imagesc(tt * 170, ff / 1000, I); % distance as x-axis
            xlabel('distance / m'); ylabel('frequency / kHz');
    
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff / 1000, I); % time as x-axis
            % xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("2.2 Spectrogram of Filtered Signal (550~2919) 12x48")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 2.3.1 pspectrum (chirp part)
    % power spectrum, but also supports spectrogram, power spectrum and persistence spectrum
    % plot the chirp part (1~500)
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        [ss,ff,tt] = pspectrum(chunk_chirp, fs, 'spectrogram', ...
                        'FrequencyLimits', [14000 20000], ...
                        'TimeResolution', 96/44100, ...
                        'OverlapPercent', 90);
        ff = ff(:);
        ss = ss(:,:);
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);

        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            % h{idx}{plt_idx} = imagesc(tt * 170, ff / 1000, I); % distance as x-axis
            % xlabel('distance / m'); ylabel('frequency / kHz');
    
            h{idx}{plt_idx} = imagesc(tt * 1000, ff / 1000, I); % time as x-axis
            xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("2.3.1 Power Spectrogram of Filtered Signal (1~500)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 2.3.2 pspectrum (chirp part, amplified)
    % plot the chirp part, amplify by M to light up nearby reflections
    if ismember(idx, [])
        plt_idx = plt_idx + 1;

        % light up nearby reflections, please adjust M
        M = 3;
        S = min(S * M, 1);
        I = mat2gray(S);

        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
    
            h{idx}{plt_idx} = imagesc(tt * 1000, ff / 1000, I); % time as x-axis
            xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("2.3.2 Power Spectrogram of Filtered Signal (1~500, amplified)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 2.3.3 pspectrum (echo part)
    % plot the echo part (550~2919)
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_echo = chunk_plot(550:2920-1);

        [ss,ff,tt] = pspectrum(chunk_echo, fs, 'spectrogram', ...
                        'FrequencyLimits', [14000 20000], ...
                        'TimeResolution', 4*96/44100, ...
                        'OverlapPercent', 90);
        ff = ff(:);
        ss = ss(:,:);
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);

        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            h{idx}{plt_idx} = imagesc((tt+550/44100) * 170, ff / 1000, I); % distance as x-axis
            xlabel('distance / m'); ylabel('frequency / kHz');
    
            title("2.3.3 Power Spectrogram of Filtered Signal (550~2919)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end
    
    %% 3.1 FMCW distance measurement
    % FMCW distance measurement with power spectrum (not spectrogram)
    % 'power' option in pspectrum is used
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));
        start_time = 700;
        interp_rate = 2;
        chunk_interp = spline(1:3000, chunk_plot(1:3000), 1 : 1/interp_rate : 3000+1/interp_rate);

        f0 = 15000;
        f1 = 20000 + 5000 * 5;
        fs = 44100 * 2;
        w0 = f0 / fs * 2 * pi;
        w1 = f1 / fs * 2 * pi;
        sig_len = 6000;
        my_chirp = zeros(1, sig_len);
        for i = 1:sig_len
            K = sig_len * w0 / log(w1/w0);
            L = sig_len / log(w1/w0);
            phase = K * (exp(i/L) - 1.0);
            my_chirp(i) = sin(phase);
        end
        disp([length(chunk_interp), length(my_chirp)])

        chunk_IF = chunk_interp .* my_chirp; % use plot_chirp.m to generate chirp
        chunk_chirp = chunk_IF(1:start_time * interp_rate);
        chunk_echo = chunk_IF(start_time * interp_rate : 2920 * interp_rate-1);

        % fs = 44100 * 2
        [ss,ff,tt] = pspectrum(chunk_echo, fs, 'power', ...
                        'Reassign', false, ...
                        'FrequencyLimits', [0 22050], ...
                        'FrequencyResolution',1000);
                        
%                         'OverlapPercent', 90);
%                         'TimeResolution', 2*96/fs, ...

%                         'Reassign', false, ...
%                         'FrequencyLimits', [14000 20000], ...
        ff = ff(:);
        ss = ss(:,:);
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);

        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            h{idx}{plt_idx} = imagesc((tt + start_time/44100), ff*170*(500+start_time)/(5000*44100), I); % distance as x-axis
            xlabel('time / s'); ylabel('distance / m');
    
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff / 1000, I); % time as x-axis
            % xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("PSpectrogram of Filtered Echoes (700~2919)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 0. add new analysis script here
    

end

% callback function on building connections
function connectionInfo(src,~)
    global server
    for idx = 1:3
        if src == server{idx}
            prefix = strcat("[", num2str(idx), "] ");
            break
        end
    end
    if src.Connected
        serverInfo = strcat(src.ServerAddress, ":", num2str(src.ServerPort));
        clientInfo = strcat(src.ClientAddress, ":", num2str(src.ClientPort));
        disp(strcat(prefix, "Connected to ", serverInfo, " from client ", clientInfo))
    else
        serverInfo = strcat(src.ServerAddress, ":", num2str(src.ServerPort));
        disp(strcat(prefix, "No client on ", serverInfo))
        flush(src)
    end
end

% function to reinitialize servers
function clearData(idx_list)
    global data data_filtered start_list
    for idx = idx_list
        data{idx}           = [];
        data_filtered{idx}  = [];
        start_list{idx}     = [];
    end
end