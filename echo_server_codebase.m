% EchoServer Codebase
% a collection of unused code blocks

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
    

    %% 2.3 STFT (chirp part, higher resolution)
    % short-time fourier transform, resolution will be shown in the title
    % this is to monitor the chirp

    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        chunk_chirp = chunk_full(1:550);
        % chunk_echo = chunk_full(550:2920-1);

        n_window = 128;
        n_overlap = 31/32 * n_window;
        n_fft = n_window * 8;

        [ss,ff,tt] = stft(chunk_chirp, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
        
        % frequency range you're interest in
        f0 = 14000; f1 = 21000;
        ss = ss((ff>f0 & ff<f1),:);

        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);

            % chirp (and echo), distance as axis
            h{idx}{plt_idx} = imagesc(tt * 170, ff(ff>f0 & ff<f1) / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % chirp (and echo), time as axis
            % h{idx}{plt_idx} = imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I); xlabel('time / ms'); ylabel('frequency / kHz');

            title(strcat("2.3 STFT of Filtered Signal (1~550) ", num2str(size(I,1)), 'x', num2str(size(I,2))))
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 2.4 STFT (echo part, higher resolution)
    % short-time fourier transform, resolution will be shown in the title
    % this is to monitor the chirp

    if ismember(idx, [])
        plt_idx = plt_idx + 1;
        % chunk_chirp = chunk_full(1:550);
        chunk_echo = chunk_full(550:2920-1);

        n_window = 128;
        n_overlap = 31/32 * n_window;
        n_fft = n_window * 8;

        [ss,ff,tt] = stft(chunk_echo, fs, 'Window', tukeywin(n_window,0.25), ...
            'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
        
        % frequency range you're interest in
        f0 = 14000; f1 = 21000;
        ss = ss((ff>f0 & ff<f1),:);

        S = abs(ss) / max(max(abs(ss)));
        I = mat2gray(S);
        
        if length(f{idx}) < plt_idx || ~isvalid(f{idx}{plt_idx})
            f{idx}{plt_idx} = figure(idx*100+plt_idx);

            % echo part, distance as axis
            h{idx}{plt_idx} = imagesc((tt+550/44100) * 170, ff(ff>f0 & ff<f1) / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % echo part, time as axis
            % h{idx}{plt_idx} = imagesc((tt+550/44100) * 1000, ff(ff>f0 & ff<f1) / 1000, I); xlabel('time / ms'); ylabel('frequency / kHz');

            title(strcat("2.4 STFT of Filtered Signal (550~2919) ", num2str(size(I,1)), 'x', num2str(size(I,2))))
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 3.1 spectrogram (echo part, 12x48)
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
    
            % echo part, distance as axis
            h{idx}{plt_idx} = imagesc((tt+550/44100) * 170, ff(ff>f0 & ff<f1) / 1000, I); xlabel('distance / m'); ylabel('frequency / kHz');
            % echo part, time as axis
            % h{idx}{plt_idx} = imagesc((tt+550/44100) * 1000, ff(ff>f0 & ff<f1) / 1000, I); xlabel('time / ms'); ylabel('frequency / kHz');
    
            title("3.1 Spectrogram of Filtered Signal (550~2919) 12x48")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 4.1 pspectrum (chirp part)
    % power spectrum, but also supports spectrogram, power spectrum and persistence spectrum
    % plot the chirp part (1~500)
    % archived in echo_server_codebase.m


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
    
            title("4.1 Power Spectrogram of Filtered Signal (1~500)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 4.2 pspectrum (chirp part, amplified)
    % plot the chirp part, amplify by M to light up nearby reflections
    % archived in echo_server_codebase.m


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
    
            title("4.2 Power Spectrogram of Filtered Signal (1~500, amplified)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end


    %% 4.3 pspectrum (echo part)
    % plot the echo part (550~2919)
    % archived in echo_server_codebase.m


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
    
            title("4.3 Power Spectrogram of Filtered Signal (550~2919)")
            set(gca, 'YDir', 'normal');
            colormap('hot')
        end
        set(h{idx}{plt_idx}, 'CData', I);
    end
    

    %% 5.1 FMCW distance measurement
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

