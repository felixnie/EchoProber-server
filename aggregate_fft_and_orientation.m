% EchoServer - Aggregate fft_list & ori_list in a Timetable
idx = 3;
load sensorlog.mat
load fltl.mat
fl = fft_list{idx};
tl = time_list{idx};
N = length(tl);

% fft timetable
% ft = array2timetable(fl, 'RowTimes', tl);
% % time range
% tr = timerange(tl(1), tl(end));

ot = Orientation;
ot = retime(ot, tl, 'linear');

% check range?
disp([tl(1), tl(end)])
disp([ot(1,:).Timestamp, ot(end,:).Timestamp])

offset = 0; % degree, azimuth when arrow points to 0

azimuth = (ot.X - offset) * pi / 180;
elevation = -ot.Y * pi / 180;

echo_start = 600; fs = 44100;
rr = ff * 170 * 500 / (44100*5000) + echo_start / fs * 170;
R = length(rr);

R = length(rr);
X = [];
Y = [];
Z = [];

for i = 1:N
    [x,y,z] = sph2cart(azimuth(i) * ones(R, 1), elevation(i) * ones(R, 1), rr(:));
    X = [X; x];
    Y = [Y; y];
    Z = [Z; z];
end

C = fl';
C = C(:);

I = 1:N * R;
I = I(C>0.7);
% I = I(I<3000000);
% C = C(C>mean(C));

NR = N * R;

S = ones(NR, 1) * 10;

% % downsample
% I = downsample(1:NR, 100);
% I = I(I<3000000);

figure(1)
alpha = 0.05;
h = scatter3(X(I), Y(I), Z(I), S(I), C(I), 'o', 'filled', ...
    'MarkerEdgeAlpha', alpha, 'MarkerFaceAlpha', alpha);
colormap(jet);
colorbar;
set(gca,'XLim',[-20 20],'YLim',[-20 20],'ZLim',[-20 20])
view([0 90])