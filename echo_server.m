% EchoServer
% 06 Mar 2022
% 
% EchoServer 3 implements:
% 1. Continuous data receiving and plotting.
% 2. Resolving both short message and long data from clients.
% 3. Multiple client connection management.
% 4. Multiple figure background refresh mechanism.
% 5. Server to client messaging.
% 
% Notes:
% Please upgrade to MATLAB R2021b or later to use tcpserver.
% 
% To-do:
% 1. Parallel processing for analysis tasks.

clear server
close all

% set up server
localhost = "10.25.245.22";

global data data_filtered start_list f h server
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
            disp(strcat(prefix, "Please disconnect from Android client"));
            disp(strcat(prefix, "Please check host IP"));
            % command = "sudo kill 'sudo lsof -t -i:8170'";
            % [status, cmdout] = system(command);
            return
        otherwise
            rethrow(ME)
    end
end

% callback function on receiving data
function readData(src,~)
    global data data_filtered start_list f h server

    % get server index
    for idx = 1:3
        if src == server{idx}
            prefix = strcat("[", num2str(idx), "] ");
            break
        end
    end

    % receive data
    bytes_received = src.NumBytesAvailable;
    
    if  bytes_received == 0
        % do nothing
        % to-do: check where are empty packets from, with wireshark
        % disp(strcat(prefix, "Bytes received: ", num2str(bytes_received), " - ", ...
        %                     "Empty packet"))
        flush(src)
        return
    elseif bytes_received < 20
        % short message
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
    else
        % long data
        src.UserData = read(src, (bytes_received-2)/2, 'int16'); % drop \r\n
        disp(strcat(prefix, "Bytes received: ", num2str(bytes_received), " - ", ...
                            "Data length: ", num2str(length(src.UserData))))
        flush(src)
    end

    % set up filter
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
  
    % find the first start point
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

    % find the rest of the start points
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

    if length(start_list{idx}) > start_processed % new start point added
        if length(start_list{idx}) > 1 % continue
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: ", num2str(start_list{idx}(end) - start_list{idx}(end-1)), " - ", ...
                                "New start point added."))
        else % have only 1 start point
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: N/A", " - ", ...
                                "No enough start point."))
            % do nothing
            return
        end
    else % no new start point
        if length(start_list{idx}) > 1
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: ", num2str(start_list{idx}(end) - start_list{idx}(end-1)), " - ", ...
                                "No new start point."))
            % do nothing
            return
        else % have only 1 start point
            disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - ", ...
                                "Last ending: ", num2str(start_list{idx}(end)), " - ", ...
                                "Last length: N/A", " - ", ...
                                "Wait for more data."))
            % do nothing
            return
        end
    end
    
    % initialize plot index
    plt_idx = 0;

    % Below are the scripts being called everytime when new data chunk is 
    % arrived. Change the index array to enable/disable each part. For
    % example:
    %     if ismember(idx, [1 2 3]) % enable script for device 1, 2, 3
    %     if ismember(idx, [])      % disable script for all devices
  
    %% 1. plot signal (data_filtered, full chunk with length 4000)
    if ismember(idx, [1 2 3])
        plt_idx = plt_idx + 1; % plot filtered signal
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end)); % chunk to plot
    
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            % figure not defined yet or figure being closed
            f{idx}{plt_idx} = figure(idx*100+plt_idx); % create or recreate figure handler
            h{idx}{plt_idx} = plot(chunk_plot); % create or recreate plot handler
            title("Filtered Signal (1~4000)")
            axis([0 4000 -10000 10000])
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', chunk_plot); % if figure defined and opened, reload data
    end

    %% 2. plot signal (data_filtered, echo part with length 2370)
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));
        chunk_echo = chunk_plot(550:2920-1);
    
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            h{idx}{plt_idx} = plot(550:2920-1, chunk_echo);
            title("Filtered Signal (550~2919)")
            axis([550 2920 -2000 2000])
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', chunk_echo);
    end
    
    %% 3. cross-correlation (need to re-write)
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));
        
        % %%%%%%%%%% chunk_chirp %%%%%%%%%%
        sig_len = 500;
        chunk_chirp = zeros(1, sig_len);
        
        f0 = 15000;
        f1 = 20000;
        w0 = f0 / fs * 2 * pi;
        w1 = f1 / fs * 2 * pi;
        
        for i = 1:sig_len
            K = sig_len * w0 / log(w1/w0);
            L = sig_len / log(w1/w0);
            phase = K * (exp(i/L) - 1.0);
            chunk_chirp(i) = sin(phase);
        end
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);
    
        % chunk_chirp_win = chunk_chirp .* hanning(length(chunk_chirp))';
        [c, lags] = xcorr(chunk_chirp, chunk_echo);

        % r = (1:2920-500) + length(chunk_echo); % ?
        % x = lags(r) / fs * 340 / 2;
        % [yupper, ~] = envelope(c(r), 10, 'rms');
        
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % h{idx}{plt_idx} = plot(x, yupper);
            h{idx}{plt_idx} = plot(lags, c);
            xlabel('lag in distance / m'); ylabel('cross-corr');
            title("Envelope of Cross-correlation")
            axis([lags(1) lags(end) -max(abs(c)) max(abs(c))])
            grid on
        end
        % set(h{idx}{plt_idx}, 'ydata', yupper);
        set(h{idx}{plt_idx}, 'ydata', c);
    end

    %% 4. spectrogram 12x48
    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));
        chunk_chirp = chunk_plot(1:500);
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
    
            title("Spectrogram of Filtered Signal (550~2919)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 5. STFT 12x48
    if ismember(idx, [1 2 3])
        plt_idx = plt_idx + 1;
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));
        chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        n_window = 96;
        n_overlap = 48;
        n_fft = 96;
        [ss,ff,tt] = stft(chunk_echo, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
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
    
            title("STFT of Filtered Signal (550~2919)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 6. pspectrogram
    if ismember(idx, [1 2 3])
        plt_idx = plt_idx + 1;
        chunk_plot = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end));
        chunk_chirp = chunk_plot(1:500);
        chunk_echo = chunk_plot(550:2920-1);

        [ss,ff,tt] = pspectrum(chunk_echo, fs, 'spectrogram', ...
                        'FrequencyLimits', [14000 20000], ...
                        'TimeResolution', 16/44100, ...
                        'OverlapPercent', 50);
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
    
            h{idx}{plt_idx} = imagesc(tt * 170, ff / 1000, I); % distance as x-axis
            xlabel('distance / m'); ylabel('frequency / kHz');
    
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff / 1000, I); % time as x-axis
            % xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("PSpectrogram of Filtered Signal (550~2919)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 7. add new script here

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