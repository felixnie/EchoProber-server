% EchoLocal
% 12 Jul 2022
%
% EchoServer, but the offline version.

close all
clear

%% original chirp signal
sig_len = 500;
original_chirp = zeros(1, sig_len);
fs = 44100;
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
% w = hann(100)';
% w = [w(1:50) ones(1, sig_len - 100) w(51:100)];
% original_chirp = w .* original_chirp;
% figure
% plot(original_chirp)


%% long chirp signal under 2*fs
% sig_len = 6000;
% long_chirp = zeros(1, sig_len);
% fs_ = 44100 * 2;
% f0 = 15000;
% f1 = 20000 + 5000 * 5;
% w0 = f0 / fs_ * 2 * pi;
% w1 = f1 / fs_ * 2 * pi;
% for i = 1:sig_len
%     K = sig_len * w0 / log(w1/w0);
%     L = sig_len / log(w1/w0);
%     phase = K * (exp(i/L) - 1.0);
%     long_chirp(i) = sin(phase);
% end
Fs = fs * 2;
long_chirp = chirp(1/Fs:1/Fs:6000/Fs, 15000, 6000/Fs, 15000+5000*6);


%% full chirp signal
play_seconds = input("Seconds to play (10 by default): ");
if isempty(play_seconds)
    play_seconds = 10;
end

duration = 4000;
repeat_times = play_seconds * fs / duration;
repeat_times = ceil(repeat_times);

full_chirp = zeros(1, duration * repeat_times);
for i = 1:repeat_times
    r = (i-1)*duration+1 : (i-1)*duration+sig_len;
    full_chirp(r) = original_chirp;
end


% % save full chirp signal
% filename = ['full_chirp_', num2str(repeat_times), '.wav'];
% audiowrite(filename, full_chirp, fs);


%% show device list
info = audiodevinfo;

disp('================== Input ===================')
for i = 1:length(info.input)
    disp(info.input(i))
    if strcmp(info.input(i).Name(6:12), 'Realtek')
        iID_default = info.input(i).ID;
    end
end

disp('================== Output ==================')
for i = 1:length(info.output)
    disp(info.output(i))
    if strcmp(info.output(i).Name(5:13), 'Soundcore')
        oID_default = info.output(i).ID;
    end
end

iID = input("Choose an input device: ");
if isempty(iID)
    iID = iID_default;
end
oID = input("Choose an output device: ");
if isempty(oID)
    oID = oID_default;
end

%% setup recorder and player
% https://www.mathworks.com/help/matlab/ref/audiorecorder.html#d123e57619
fs = 44100;
nBits = 16;
nChannels = 1;
recorderObj = audiorecorder(fs, nBits, nChannels, iID);
playerObj = audioplayer(full_chirp, fs, nBits, oID);



%% start recorder and player
tic

sound(full_chirp, fs)
disp("start player")
toc

record(recorderObj)
disp("start recorder")
toc

tic
last_len = 0;
while isrecording(recorderObj)
    % toc

    if toc > play_seconds
        stop(recorderObj)
    end

    try
        data = getaudiodata(recorderObj);
        current_len = sum(abs(data) > 0.00001);
        disp(current_len)
    catch ME
        disp(ME)
        continue
    end

%     if current_len - last_len > duration
%         current_data = data(last_len+1:current_len);
%         last_len = current_len;
%     end
end


%% find start points

fc = 10000; fs = 44100;
[b,a] = butter(5, fc/(fs/2), 'high');
data = filter(b, a, data);

figure
plot(data)

threshold = 0.2*max(abs(data));
start_point = [];

i = 1;
while i < length(data)-duration
    if abs(data(i)) > threshold
        start_point = [start_point i];
        i = i + 3800;
    else
        i = i + 1;
    end
end

hold on
stem(start_point, threshold * ones(1, length(start_point)))

input('Are they properly clipped? ENTER to continue: ')
%% analysis
start_point_idx = 5:1:length(start_point);
GIF_MODE = true;
fft_list = zeros(3001, length(start_point_idx));

for i = 1:length(start_point_idx)
    % init fig_idx and chunks
    fig_idx = 0;

    chunk_full = data(start_point(start_point_idx(i)) : ...
                      start_point(start_point_idx(i))+4000-1);

    chunk_chirp = chunk_full(1:500);
    chunk_echo = chunk_full(550:2919);

    %% 1. STFT (enhanced, full)
    n_window = 128;
    n_overlap = 31/32 * n_window;
    n_fft = n_window * 8;

    [ss,ff,tt] = stft(chunk_full, fs, 'Window', tukeywin(n_window,0.25), ...
        'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');

    ss = ss((ff>14000 & ff<21000), :); % frequency range you're interest in
    S = abs(ss) / max(max(abs(ss)));
    
    % enhancement
    echo_start = 400; echo_end = 4000;
    ss = S(:, (tt>echo_start/44100 & tt<echo_end/44100));
    S(:, (tt>echo_start/44100 & tt<echo_end/44100)) = abs(ss) / max(max(abs(ss)));
    
    echo_start = 600; echo_end = 4000;
    ss = S(:, (tt>echo_start/44100 & tt<echo_end/44100));
    S(:, (tt>echo_start/44100 & tt<echo_end/44100)) = abs(ss) / max(max(abs(ss)));

    echo_start = 800; echo_end = 2919;
    ss = S(:, (tt>echo_start/44100 & tt<echo_end/44100));
    S(:, (tt>echo_start/44100 & tt<echo_end/44100)) = abs(ss) / max(max(abs(ss)));

    I = mat2gray(S);
    
    fig_idx = fig_idx + 1;
    fig = figure(start_point_idx(i) * 100 + fig_idx);
%     imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I);
%     xlabel('time / ms'); ylabel('frequency / kHz');
    imagesc(tt * 170, ff(ff>f0 & ff<f1) / 1000, I);
    xlabel('distance / m'); ylabel('frequency / kHz');
    title(sprintf("STFT of full signal (chirp: 1~500, echo: 550~2919) idx: %d", i))
    set(gca, 'YDir', 'normal');
    colormap('hot')
    
    if GIF_MODE
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        [A,map] = rgb2ind(im{i},256);
        filename = '01.gif';
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    
        close(figure(start_point_idx(i) * 100 + fig_idx))
    end

    %% 2. STFT (chirp or echo, full spectogram or parted)
    n_window = 128;
    n_overlap = 31/32 * n_window;
    n_fft = n_window * 8;

    [ss,ff,tt] = stft(chunk_chirp, fs, 'Window', tukeywin(n_window,0.25), ...
        'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');

    ss = ss((ff>14000 & ff<21000),:); % frequency range you're interest in

    S = abs(ss) / max(max(abs(ss)));
    I = mat2gray(S);
    
    fig_idx = fig_idx + 1;
    fig = figure(start_point_idx(i) * 100 + fig_idx);
%     imagesc((tt+550/44100) * 1000, ff(ff>f0 & ff<f1) / 1000, I);
%     xlabel('time / ms'); ylabel('frequency / kHz');
    imagesc(tt * 1000, ff(ff>f0 & ff<f1) / 1000, I);
    xlabel('time / ms'); ylabel('frequency / kHz');
    title(sprintf("STFT of chirp signal (1~500) idx: %d", i))
    set(gca, 'YDir', 'normal');
    colormap('hot')

    if GIF_MODE
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        [A,map] = rgb2ind(im{i},256);
        filename = '02.gif';
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    
        close(figure(start_point_idx(i) * 100 + fig_idx))
    end

    %% 3. FMCW


    interp_rate = 2;
    echo_start = 600; echo_end = 3000;
    chunk_enhance = chunk_full(1:echo_end);
    chunk_enhance(echo_start:echo_end) = chunk_enhance(echo_start:echo_end) * ...
                              max(abs(chunk_chirp)) / ...
                              max(abs(chunk_enhance(echo_start:echo_end)));
    chunk_interp = spline(1:echo_end, ...
                          chunk_enhance(1:echo_end), ...
                          1 : 1/interp_rate : echo_end + 1/interp_rate);

    chunk_IF = chunk_interp .* long_chirp;
%     chunk_chirp = chunk_IF(1 : echo_start * interp_rate);
%     chunk_echo = chunk_IF(echo_start * interp_rate : echo_end * interp_rate-1);

    n_window = 128;
    n_overlap = 31/32 * n_window;
    n_fft = n_window * 8;

    [ss,ff,tt] = stft(chunk_IF, Fs, 'Window', tukeywin(n_window,0.25), ...
        'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');

%     ss = ss((ff>14000 & ff<21000),:); % frequency range you're interest in

    S = abs(ss) / max(max(abs(ss)));
    I = mat2gray(S);

    fig_idx = fig_idx + 1;
    fig = figure(start_point_idx(i) * 100 + fig_idx);
    imagesc(tt * 1000, ff / 1000, I);
    xlabel('time / ms'); ylabel('frequency / kHz');
%     imagesc(tt * 1000, ff * 170 * 500 / (44100 * 5000), I);
%     xlabel('time / ms'); ylabel('distance / m');
    title(sprintf("STFT of FMCW IF signal (1~3000) idx: %d", i))
    set(gca, 'YDir', 'normal');
    colormap('hot')

    if GIF_MODE
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        [A,map] = rgb2ind(im{i},256);
        filename = '03.gif';
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    
        close(figure(start_point_idx(i) * 100 + fig_idx))
    end


    %% 4. FFT
    [b,a] = butter(5, 22100/(44100*2/2), 'low');
    y = filter(b, a, chunk_IF);
    L = length(y);
    Y = fft(y);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    f = f * 170 * 500 / (44100 * 5000);

    fig = figure(start_point_idx(i) * 100 + fig_idx);
    plot(f,P1) 
    title(sprintf('Amplitude spectrum of IF signal idx: %d', i))
%     xlabel('f (Hz)')
    xlabel('distance / m')
    ylabel('amplitude')

    if GIF_MODE
        drawnow
        frame = getframe(fig);
        im{i} = frame2im(frame);
        [A,map] = rgb2ind(im{i},256);
        filename = '04.gif';
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    
        close(figure(start_point_idx(i) * 100 + fig_idx))
    end

    fft_list(:, i) = P1';

%     % STFT of long chirp
%     new_chirp = chirp(1/(fs*2):1/(fs*2):6000/(fs*2), 15000, 6000/(fs*2), 15000+5000*6);
%     [ss,ff,tt] = stft(new_chirp, fs * interp_rate, 'Window', tukeywin(n_window,0.25), ...
%         'OverlapLength', n_overlap, 'FFTLength', n_fft, 'FrequencyRange', 'onesided');
% 
% %     ss = ss((ff>14000 & ff<21000),:); % frequency range you're interest in
% 
%     S = abs(ss) / max(max(abs(ss)));
%     I = mat2gray(S);
% 
%     fig_idx = fig_idx + 1;
%     fig = figure(start_point_idx(i) * 100 + fig_idx);
% %     imagesc(tt * 1000, ff / 1000, I);
% %     xlabel('time / ms'); ylabel('frequency / kHz');
%     imagesc(tt * 1000, ff / 1000, I);
%     xlabel('time / ms'); ylabel('frequency / kHz');
%     title(sprintf("STFT of FMCW IF signal (1~3000) idx: %d", i))
%     set(gca, 'YDir', 'normal');
%     colormap('hot')
end

figure
fft_list = fft_list ./ max(max(fft_list));
fft_list(200:end) = fft_list(200:end) ./ max(max(fft_list(200:end)));
I = mat2gray(fft_list);
imagesc(1:length(start_point_idx), f, I);
title(sprintf("FFT results by time"))
set(gca, 'YDir', 'normal');
xlabel('sample / n'); ylabel('distance / m')
colormap('hot')