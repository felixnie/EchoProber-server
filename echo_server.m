%% EchoServer
% 15 Aug 2022
%
% 13 Aug 2022 Updates:
% 1. Tested for EchoProber-iOS. Set duration to 4800 samples to comply with
%    iOS output pace.
% 2. Add frequency domain cross-correlation.
% 3. Add STFT binarization and filtering.
% 
% 25 Jul 2022 Updates:
% 1. New global variable time_list and freq to store timestamps of each frame.
%
% 15 Jul 2022 Updates:
% 1. GIF animation code added. See the example at 1.1.
% 2. Some old scripts are placed at echo_server_codebase.m.
% 3. New global variable fft_list to store fft_IF from each frame.
% 4. A linear chirp will upload to clients upon connecting.
% 
% This EchoServer implements:
% 1. Real-time data receiving and plotting.
% 2. Messaging between server and clients.
% 3. Parsing short socket as message and long data as recordings.
% 4. Multiple clients connection management.
% 5. Multiple figures refreshing mechanism.
% 6. Automatically disconnect upon restarting server.
% 
% Notes:
% Please upgrade to MATLAB R2021b or later to use tcpserver.
% 
% To-do:
% 1. Parallel processing for analysis tasks.
% 2. Parallel processing for each analysis tasks.
% 3. Light mode to save memory.
% 4. Label with device orientation. (Used MATLAB Drive instead)

clear server
close all


%% configuation

% % in light mode, data and data_filtered will not keep growing forever
% light_mode = false;
% buffer_size = 44100; % fix data and data_filtered to this size
% 
% in GIF_mode, a series of frames will be saved as a .gif file
gif_mode = false;

% set up server
% localhost = "192.168.1.103";
localhost = "155.69.142.8";


%% initialize global variables

% this will create empty variables
global f h server data data_filtered start_list fft_list time_list freq
% f:                (cell of arrays of figure) to hold the figures
% h:                (cell of arrays of plot) the handles for figures
% server:           (cell of objects) the tcpserver objects
% data:             (cell of arrays of float) the float data received 
%                                             between a 'play' and 'stop'
% data_filtered:    (cell of arrays of float) the low-pass filtered data
% start_list:       (cell of arrays of int) the index of each chirp arrival
% fft_list:         (cell of matrix of float) the fft results
% time_list:        (cell of arrays of timestamp) the timestamps
% freq:             (array of float) frequency points

% ask existing clients to disconnect
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

% support multiple clients (3 by default), please assign the ports manually
data            = cell(1,3);    % all data received
data_filtered   = cell(1,3);    % all data filtered
start_list      = cell(1,3);    % list of start point index
fft_list        = cell(1,3);    % list of fft results
time_list       = cell(1,3);    % list of timestamps
f               = cell(1,3);    % figure list
h               = cell(1,3);    % plot list
server          = cell(1,3);    % server list
port            = 8170:8172;    % device 1~3 on 8170, 8171 and 8172

% initialization
try
    for idx = 1:3
        data{idx}               = [];
        data_filtered{idx}      = [];
        start_list{idx}         = [];
        fft_list{idx}           = [];
        time_list{idx}          = [];
%         if light_mode
%             data{idx}           = zeros(1, buffer_size);
%             data_filtered{idx}  = zeros(1, buffer_size);
%         end
        % set callback functions for connection change and receiving data
        server{idx} = tcpserver(localhost, port(idx));
        server{idx}.ConnectionChangedFcn = @connectionInfo;
        configureTerminator(server{idx},"CR/LF")
        configureCallback(server{idx}, "terminator", @readData);
    end
catch ME
    switch ME.identifier
        case 'instrument:interface:tcpserver:cannotConnect'
            % in case of connection fail, ask user to retry
            prefix = strcat("[", num2str(idx), "] ");
            disp(strcat(prefix, "Connection failed. Please try again."));
            return
        otherwise
            rethrow(ME)
    end
end

% %% set up device
% 
% if ~exist('m', 'var')
%     if ~exist('m_list', 'var')
%         disp('[I] Loading device list')
%         m_list = mobiledevlist;
%         disp(m_list)
%     end
%     if size(m_list, 1) == 0
%         disp('[W] No device')
%         return
%     elseif size(m_list, 1) == 1
%         m = mobiledev;
%     elseif size(m_list, 1) > 1
%         i = input('[I] Select a device: ');
%         m = mobiledev(m_list.Device{i});
%     end
%     disp(['[I] Current device: ' m.Device])
% end
% 
% while ~m.Logging
%     disp('[I] Device has not started logging yet')
%     input('[I] START on smartphone then press ENTER: ')
% end



%% callback function for receiving data

function readData(src,~)
    % the reason of placing all data receiving and analysis code in one
    % place is to avoid data copying in MATLAB
    global f h server data data_filtered start_list fft_list time_list freq


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

    elseif bytes_received < 100 % short message

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
    fs                  = 44100;
    
    % duration          : designed duration between two start points
    %                     this will affect how we search the start points
    duration            = 4800;

    % chunk             : new raw data
    % chunk_filtered    : new filtered data
    % data              : all raw data
    % data_filtered     : all filtered data
    
    chunk               = src.UserData;
    data{idx}           = [data{idx}, chunk];
    
    % filter has to deal with a growing array. changed to only filter a part of data
    % data_filtered{idx}  = filter(b, a, data{idx}); 
    duration_to_filter  = 9600;
    start_of_filtered   = max(1, length(data{idx}) - duration_to_filter);
    [b,a]               = butter(5, 14000/(fs/2), 'high'); % high pass filter
    chunk_filtered      = filter(b, a, data{idx}(start_of_filtered : end));

    % strategy 1: add the filtered new data to data_filtered
    % chunk_filtered      = chunk_filtered(length(chunk_filtered)-length(chunk)+1 : end);
    % data_filtered{idx}  = [data_filtered{idx}, chunk_filtered];

    % strategy 2: also update the last chunk of data_filtered
    % this is to smooth the overall signal
    data_filtered{idx}  = [data_filtered{idx}, chunk]; % add some stuff to the end
    data_filtered{idx}(start_of_filtered : end) = chunk_filtered; % update

    % maximum_dalay     : first start point should be within maximum_dalay
    %                     usually set to 3x duration
    maximum_delay       = 12000;

    % start_processed   : last processed/displayed chunk index
    %                     sometimes more than 1 chirp chunk will arrive
    %                     within a data chunk when data chunk size
    %                     fluctuates. only the last chunk will be processed
    start_processed     = length(start_list{idx});

    % allow_silence     : display no matter if it's a chirp
    allow_silence       = true;
    % silence_threshold : the amplitude threshold when there's no chirp
    silence_threshold   = 500; % can be determined through the delay part

    % check if there's enough data for finding start point
    if length(data{idx}) < maximum_delay
        disp(strcat(prefix, "Data length: ", num2str(length(data_filtered{idx})), " - Lack of data, skip."))
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
        % 2x maximum background noise, 10 is the safe distance to chirp
        silence_threshold = 2 * max(abs(data_filtered{idx}(1 : start_point-10)));
        disp(strcat(prefix, "Silence threshold = ", num2str(silence_threshold)))
    end


    %% find the rest of the start points

    while length(data_filtered{idx}) >= start_list{idx}(end) + duration

        % set dynamic threshold strategy

        % strategy 1: find maximum in filtered chunk
        % dynamic_threshold = max(chunk_filtered) / 5;
        % dynamic_threshold = max([dynamic_threshold silence_threshold]); % will not find start when there's no chirp
        % dynamic_threshold = silence_threshold;

        % strategy 2: find maximum in last 12000 points
        dynamic_threshold = max(abs(data_filtered{idx}(end-maximum_delay+1 : end))) / 5;
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

    % ========== some data that will be used in all code blocks ===========

    % current full chunk
    % chunk_full = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end)); % may not be exactly 4000
    chunk_full = data_filtered{idx}(start_list{idx}(end-1) : start_list{idx}(end-1)+duration-1); % fixed length

    % % original chirp signal (original log chirp)
    % sig_len = 500;
    % original_chirp = zeros(1, sig_len);
    % f0 = 15000;
    % f1 = 20000;
    % w0 = f0 / fs * 2 * pi;
    % w1 = f1 / fs * 2 * pi;
    % for i = 1:sig_len
    %     K = sig_len * w0 / log(w1/w0);
    %     L = sig_len / log(w1/w0);
    %     phase = K * (exp(i/L) - 1.0);
    %     original_chirp(i) = sin(phase);
    % end
    
    % short linear chirp signal (original linear chirp)
    linear_chirp = chirp(1/fs:1/fs:500/fs, 15000, 500/fs, 15000+5000);

    % long linear chirp signal (2x sample rate, 6x duration)
    Fs = fs * 2;
    long_chirp = chirp(1/Fs:1/Fs:6000/Fs, 15000, 6000/Fs, 15000+5000*6);

    % % code to upload chirp signal (used in connectionInfo())
    % sig_str = sprintf('%f,', linear_chirp);
    % sig_str = sig_str(1:end-1);
    % writeline(server{device_idx}, sig_str)
    % pause(2)

    % =====================================================================


    % Below are the scripts being called everytime when new chirp chunk is 
    % arrived. Change the index array to enable/disable each part.
    % example:
    %     if ismember(idx, [1 2 3]) % enable script for device 1, 2, 3
    %     if ismember(idx, [])      % disable script for all devices
  

    %% 1.1 plot signal (chirp and echo)
    % plot time domain data after low-pass filtering
    % the length of chunk_full equals duration

    % here detailed comments are provided on:
    %   how to manage figures with f
    %   how to update plots with handlers h
    %   how to save figures as gif

    if ismember(idx, [])
        plt_idx = plt_idx + 1; % current plot index
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx}) % figure not defined yet or figure being closed
            f{idx}{plt_idx} = figure(idx*100+plt_idx); % create or recreate figure handler
            h{idx}{plt_idx} = plot(chunk_full); % create or recreate plot handler
            title("1.1 Filtered Signal (1~" + duration + ")")
            axis([0 duration -10000 10000])
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', chunk_full); % if figure defined and opened, reload data

        if gif_mode && ismember(idx, []) % saving current figure as a gif file
            drawnow
            frame = getframe(f{idx}{plt_idx});
            frame_idx = length(start_list{idx}) - 1;
            frame_img = frame2im(frame);
            [A,map] = rgb2ind(frame_img, 256);
            gif_file_name = [num2str(idx*100+plt_idx) '.gif'];
            if frame_idx == 1
                imwrite(A, map, gif_file_name, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
            else
                imwrite(A, map, gif_file_name, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
            end
        end
    end

    
    %% 1.2 plot signal (echo part)
    % plot time domain data after low-pass filtering
    % the length of chunk_echo is 2420 (2920-550)

    if ismember(idx, [])
        
        echo_start = 550;
        echo_end = 2920;
        chunk_echo = chunk_full(echo_start:echo_end);
    
        plt_idx = plt_idx + 1;
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            h{idx}{plt_idx} = plot(echo_start:echo_end, chunk_echo);
            title("1.2 Filtered Signal (" + echo_start + "~" + echo_end +")")
            axis([echo_start echo_end -2000 2000])
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', chunk_echo);
    end
    

    %% 1.3 cross-correlation
    % archived in echo_server_codebase.m
    
    
    %% 2.1 STFT (echo part, 12x48)
    % short-time fourier transform, resolution 12x48
    % this is the same as in EchoLoc source code and pcm2pickle.py

    % here detailed comments are provided on:
    %   how to update images with handlers h

    if ismember(idx, [])
        
        echo_start = 550;
        echo_end = 2920;
        chunk_echo = chunk_full(echo_start:echo_end-1);

        n_window = 96;
        n_overlap = 48;
        n_fft = 96;
        [ss,ff,tt] = stft(chunk_echo, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
        
        % capture the frequency interval f0~f1, as the Python version does
        f0 = 15000;
        f1 = 20000;
        segment = (fs/2) / n_overlap;
        lower_f = int8(f0/segment) - 1;
        upper_f = int8(f1/segment);
        ff = ff(lower_f:upper_f-1);
        ss = ss(lower_f:upper_f,:);
        
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        plt_idx = plt_idx + 1;
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            % imshow : image without axis
            % imagesc: image with axis, can rescale, no colormap
            % both support full speed real-time display on Lambda workstation
    
            % cmap = colormap('hot');
            % h{idx}{plt_idx} = imshow(I, 'Colormap', cmap, 'InitialMagnification', 1200);
    
            % echo part, distance as axis
            % h{idx}{plt_idx} = imagesc((tt+echo_start/fs) * 170, ff(ff>f0 & ff<f1) / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % echo part, time as axis
            % h{idx}{plt_idx} = imagesc((tt+echo_start/fs) * 1000, ff(ff>f0 & ff<f1) / 1000, I); xlabel('time / ms'); ylabel('frequency / kHz');
            % chirp (and echo), distance as axis
            h{idx}{plt_idx} = imagesc(tt * 170, ff(ff>f0 & ff<f1) / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % chirp (and echo), time as axis
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I); xlabel('time / ms'); ylabel('frequency / kHz');
            
            title("2.1 STFT of Filtered Signal (" + echo_start + "~" + echo_end-1 + ") 12x48")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 2.2 STFT (chirp and enhanced echo, higher resolution)
    % short-time fourier transform, resolution will be shown in the title
    % plot echo together with chirp, enhance echo part by rescaling

    if ismember(idx, [1 2 3])

        n_window = 128;
        n_overlap = 31/32 * n_window;
        n_fft = n_window * 4;

        [ss,ff,tt] = stft(chunk_full, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');

        % if you uncomment the frequency selection below, please also
        % change the y-axis settings below from 'ff(ff>f0 & ff<f1) / 1000' to 'ff / 1000'
        % f0 = 14000; f1 = 21000;
        % ss = ss((ff>f0 & ff<f1),:);

        S = abs(ss) / max(max(abs(ss)));

        % enhance echo part
        echo_start = 600; echo_end = 4000;
        ss = S(:, (tt>echo_start/fs & tt<echo_end/fs));
        S(:, (tt>echo_start/fs & tt<echo_end/fs)) = abs(ss) / max(max(abs(ss)));
        % % enhancement another part if necessary
        % echo_start = 700; echo_end = 2919;
        % ss = S(:, (tt>echo_start/fs & tt<echo_end/fs));
        % S(:, (tt>echo_start/fs & tt<echo_end/fs)) = abs(ss) / max(max(abs(ss)));

        I = mat2gray(S);
        
        plt_idx = plt_idx + 1;
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);

            % chirp (and echo), distance as axis
            h{idx}{plt_idx} = imagesc(tt * 170, ff / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % h{idx}{plt_idx} = imagesc(tt * 170, ff(ff>f0 & ff<f1) / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % chirp (and echo), time as axis
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I); xlabel('time / ms'); ylabel('frequency / kHz');
            
            title(strcat("2.2 STFT of Filtered Signal (1~4000) ", num2str(size(I,1)), 'x', num2str(size(I,2))))
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end





    %% 3.1.1 FMCW distance measurement - STFT of IF signal (No-Interp)


    if ismember(idx, [])
        
        interp_rate = 1;
        echo_start = 1; echo_end = 500;
    
        chunk_IF = chunk_full(1:echo_end) .* linear_chirp;
        % chunk_IF_chirp = chunk_IF(1 : echo_start * interp_rate);
        chunk_IF_echo = chunk_IF(echo_start * interp_rate : echo_end * interp_rate);
        
        % filter out frequency over 21kHz (ignore mirror part)
        [b,a] = butter(5, 6000/(fs*interp_rate/2), 'low');
        % chunk_IF = filter(b, a, chunk_IF); % use full signal (1~3000) afterwards
        chunk_IF = filter(b, a, chunk_IF_echo); % use echo part (600~3000) afterwards

        % filter out frequency below 1kHz (ignore chirp part)
        [b,a] = butter(5, 1000/(fs*interp_rate/2), 'high');
        chunk_IF = filter(b, a, chunk_IF);


        n_window = 256;
        n_overlap = 31/32 * n_window;
        n_fft = n_window * 2;
    
        [ss,ff,tt] = stft(chunk_IF, Fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
    
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);

        plt_idx = plt_idx + 1;
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);

            % 1. echo part, distance as x-axis, frequency as y-axis
            % h{idx}{plt_idx} = imagesc((tt + echo_start/fs) *170, ff/1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % 2. echo part, time as x-axis, distance as y-axis
            h{idx}{plt_idx} = imagesc((tt + echo_start/fs)*1000, ff*170*500/(fs*5000) + echo_start/fs * 170, I); xlabel('time / ms'); ylabel('distance / m');
            % 3. chirp (and echo), distance as x-axis, frequency as y-axis
            % h{idx}{plt_idx} = imagesc(tt*170, ff/1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % 4. chirp (and echo), time as x-axis, frequency as y-axis
            % h{idx}{plt_idx} = imagesc(tt*1000, ff/1000, I); xlabel('time / ms'); ylabel('frequency / kHz');

            title(strcat("3.1.1 STFT of IF Signal (1~600) ", num2str(size(I,1)), 'x', num2str(size(I,2))))
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end

    %% 3.1.2 FMCW distance measurement - STFT of IF signal (Interp)
    % the echo part is from echo_start = 600 to echo_end = 3000,
    % which means it's time interval 600/44100 to 3000/44100
    %
    % to make sure the linear chirp can cover time 600/44100 ~ 3000/44100,
    % the original 500/44100 long linear chirp is upsampled to 1000/88200,
    % and then extended to 6000/88200. it's stored as long_chirp
    %
    % the echo part is also upsampled (interpolated) to Fs = fs * 2
    % 
    % note that since echo part start at 600/44100, the reflection within
    % 2.31m is dropped. you can also use the whole chunk_IF rather than
    % chunk_IF_echo. but when using the full signal, the chirp part might
    % pollute the IF signal (the chirp recorded is not a perfect linear
    % chirp, it includes burst sound when the chirp start/end, etc.)

    if ismember(idx, [1 2 3])

        interp_rate = 2;
        echo_start = 600; echo_end = 3000; % echo_start = 1: no enhance

        % enhance the part of echo_start ~ echo_end in 1 ~ echo_end
        chunk_enhance = chunk_full(1:echo_end);
        chunk_enhance(echo_start:echo_end) = chunk_enhance(echo_start:echo_end) * ...
                              max(abs(chunk_full)) / ...
                              max(abs(chunk_enhance(echo_start:echo_end)));

        % upsampling
        chunk_interp = spline(1 : echo_end, ... % sample points 
                              chunk_enhance(1 : echo_end), ... % sample values
                              1 : 1/interp_rate : echo_end + 1/interp_rate); % query points
    
        chunk_IF = chunk_interp .* long_chirp;
        % chunk_IF_chirp = chunk_IF(1 : echo_start * interp_rate);
        chunk_IF_echo = chunk_IF(echo_start * interp_rate : echo_end * interp_rate);
        
        % filter out frequency over 21kHz (ignore mirror part)
        [b,a] = butter(5, 21000/(Fs/2), 'low');
        % chunk_IF = filter(b, a, chunk_IF); % use full signal (1~3000) since now
        chunk_IF = filter(b, a, chunk_IF_echo); % use echo part (600~3000) since now

        % filter out frequency below 1kHz (ignore chirp part)
        [b,a] = butter(5, 1000/(Fs/2), 'high');
        chunk_IF = filter(b, a, chunk_IF);

        % stft
        n_window = 256;
        n_overlap = 31/32 * n_window;
        n_fft = n_window * 2;
    
        [ss,ff,tt] = stft(chunk_IF, Fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
    
        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);

        plt_idx = plt_idx + 1;
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);

            % 1. echo part, distance as x-axis, frequency as y-axis
            % h{idx}{plt_idx} = imagesc((tt + echo_start/fs) *170, ff/1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % 2. echo part, time as x-axis, distance as y-axis
            h{idx}{plt_idx} = imagesc((tt + echo_start/fs)*1000, ff*170*500/(fs*5000) + echo_start/fs * 170, I); xlabel('time / ms'); ylabel('distance / m');
            % 3. chirp (and echo), distance as x-axis, frequency as y-axis
            % h{idx}{plt_idx} = imagesc(tt*170, ff/1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % 4. chirp (and echo), time as x-axis, frequency as y-axis
            % h{idx}{plt_idx} = imagesc(tt*1000, ff/1000, I); xlabel('time / ms'); ylabel('frequency / kHz');

            title(strcat("3.1 STFT of IF Signal (600~3000) ", num2str(size(I,1)), 'x', num2str(size(I,2))))
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 3.2 FMCW distance measurement - FFT of IF signal
    % this should run only when 3.1.2 is enabled because chunk_IF is needed
    % superfft is my helper function. please include superfft.m

    if ismember(idx, [1 2 3])

        [fft_IF,ff] = superfft(chunk_IF, interp_rate * (echo_end-echo_start), ...
                               Fs, 'single', 'none');

        fft_list{idx} = [fft_list{idx}; fft_IF];
        time_list{idx} = [time_list{idx}; datetime('now', 'Format', 'd-MMM-y HH:mm:ss:SSS')];
        freq = ff;
        
        plt_idx = plt_idx + 1;
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            
            % 1. distance as axis
            h{idx}{plt_idx} = plot(ff*170*500/(fs*5000) + echo_start/fs * 170, fft_IF); xlabel('distance / m'); ylabel('amplitude');
            % h{idx}{plt_idx} = plot(ff*170*500/(fs*5000), fft_IF); xlabel('distance / m'); ylabel('amplitude');
            % 2. frequency as axis
            % h{idx}{plt_idx} = plot(ff/1000, fft_IF); xlabel('frequency / kHz'); ylabel('amplitude');

            title("3.2 FFT of IF Signal (" + echo_start + "~" + echo_end + ")")
            grid on
        end
        set(h{idx}{plt_idx}, 'ydata', fft_IF);
    end


    %% 3.3 FMCW distance measurement - calculation only
    % this is a compact version of 3.1 and 3.2

    if ismember(idx, [])

        interp_rate = 2;
        echo_start = 600; echo_end = 3000;

        chunk_enhance = chunk_full(1:echo_end);
        % not to enhance
        % chunk_enhance(echo_start:echo_end) = chunk_enhance(echo_start:echo_end) * ...
        %                       max(abs(chunk_full)) / ...
        %                       max(abs(chunk_enhance(echo_start:echo_end)));

        % upsampling
        chunk_interp = spline(1 : echo_end, ...
                              chunk_enhance(1 : echo_end), ...
                              1 : 1/interp_rate : echo_end + 1/interp_rate);

        chunk_IF = chunk_interp .* long_chirp;
        chunk_IF_echo = chunk_IF(echo_start * interp_rate : echo_end * interp_rate);

        % filter out frequency over 21kHz (ignore mirror part)
        [b,a] = butter(5, 21000/(fs*2/2), 'low');

        % chunk_IF = filter(b, a, chunk_IF); % use full signal (1~3000) afterwards
        chunk_IF = filter(b, a, chunk_IF_echo); % use echo part (600~3000) afterwards

%         % filter out frequency below 1kHz (ignore chirp part)
%         [b,a] = butter(5, 1000/(fs*2/2), 'high');
%         chunk_IF = filter(b, a, chunk_IF);

        [fft_IF,ff] = superfft(chunk_IF, interp_rate * (echo_end-echo_start), ...
                               Fs, 'single', 'none');

        fft_list{idx} = [fft_list{idx}; fft_IF];
        time_list{idx} = [time_list{idx}; ...
            datetime('now', 'Format', 'd-MMM-y HH:mm:ss:SSS')];
        freq = ff;

        disp(strcat(prefix, 'Samples: ', num2str(size(time_list{idx}, 1))))

    end
    
    %% 0. orientation
    if ismember(idx, [])
        plt_idx = plt_idx + 1;


        ori = m.Orientation;
        ori(2) = -ori(2);
        quat = quaternion(ori, 'eulerd', 'ZYX', 'frame');
        ori_list{idx} = [ori_list{idx}; ori];
        

        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);
            h{idx}{plt_idx} = poseplot("NED", "MeshFileName", "phoneMesh.stl");
            xlabel("North-x (m)")
            ylabel("East-y (m)")
            zlabel("Down-z (m)");
        end
        set(h{idx}{plt_idx}, "Orientation", quat);
    end
    
    

end


%% callback function on building connections

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

        % upload chirp signal
        uploadChirp(idx)
    else
        serverInfo = strcat(src.ServerAddress, ":", num2str(src.ServerPort));
        disp(strcat(prefix, "No client on ", serverInfo))
        flush(src)
    end
end


%% function to reinitialize servers

function clearData(idx_list)
    global data data_filtered start_list fft_list
    for idx = idx_list
        data{idx}           = [];
        data_filtered{idx}  = [];
        start_list{idx}     = [];
        fft_list{idx}       = [];
    end
end

%% function to upload chirp signal

function uploadChirp(idx)
    global server
    prefix = strcat("[", num2str(idx), "] ");
    disp(strcat(prefix, "Uploading chirp."))

    fs = 44100;
    linear_chirp = chirp(1/fs:1/fs:500/fs, 15000, 500/fs, 15000+5000);

    sig_str = sprintf('%f,', linear_chirp);
    sig_str = sig_str(1:end-1);
    writeline(server{idx}, sig_str)
    disp(strcat(prefix, "Upload finished."))
end