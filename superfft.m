function [Y,f] = superfft(varargin)
% SUPERFFT  An all-in-one solution of FFT.
%   [Y,f] = SUPERFFT(y,L,Fs,SIDES,GRAPH) executes FFT with length of L and
%   sampling frequency of Fs.
%
%   [Y,f] = SUPERFFT(y,L,Fs,SIDES,GRAPH,TITLE) generates a single-side or 
%   double-side spectrum figure with your specific title.
%
%   [Y,f] = SUPERFFT(y,L,Fs,SIDES,GRAPH,XLABEL,YLABEL,TITLE) generates a
%   single-side or double-side spectrum figure.
%
%   SIDES = 'single', 'double'
%   GRAPH = 'none', 'plot', 'stem'
%
%   WARNING: If you're getting errors like this:
%       Error using horzcat
%           Dimensions of arrays being concatenated are not consistent.
%       Error in superfft (line 72)
%           P1 = [P2(floor(L/2)+2:end) P2(1:floor(L/2)+1)];
%   Try to add transposition to array y.
%
%   See also FFT.

% disp(['nargin is ' num2str(nargin)]);
% for v=1:numel(varargin)
%     disp(['varargin{' num2str(v) '} class is ' class(varargin{v})]);
% end
if nargin ~= 5 && nargin ~= 6 && nargin ~= 8
    error('The number of parameters is incorrect!!!')
end
y = varargin{1};
L = varargin{2};
Fs = varargin{3};
SIDES = varargin{4};
GRAPH = varargin{5};
if nargin == 6
    TITLE = varargin{6};
elseif nargin == 8
    XLABEL = varargin{6};
    YLABEL = varargin{7};
    TITLE = varargin{8};
end

if strcmp(SIDES,'single')
    Y = fft(y);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    Y = P1;
    if strcmp(GRAPH,'none')
        return
    elseif strcmp(GRAPH,'plot')
        figure
        plot(f,Y)
        xlabel('f/Hz');ylabel('|Y(f)|')
        title('Single-Sided Amplitude Spectrum of Y(t)')
        if nargin == 6
           title(TITLE)
        elseif nargin == 8
           xlabel(XLABEL);ylabel(YLABEL);title(TITLE)
        end
    elseif strcmp(GRAPH,'stem')
        figure
        stem(f,Y)
        xlabel('f/Hz');ylabel('|Y(f)|')
        title('Single-Sided Amplitude Spectrum of Y(t)')
        if nargin == 6
           title(TITLE)
        elseif nargin == 8
           xlabel(XLABEL);ylabel(YLABEL);title(TITLE)
        end
    else
        error('GRAPH parameter must be "none", "plot" or "stem"!!!')
    end
elseif strcmp(SIDES,'double')
    Y = fft(y);
    P2 = abs(Y/L);
    P1 = [P2(floor(L/2)+2:end) P2(1:floor(L/2)+1)];
    f = Fs*(-(L/2-1):(L/2))/L;
    Y = P1;
    if strcmp(GRAPH,'none')
        return
    elseif strcmp(GRAPH,'plot')
        figure
        plot(f,Y)
        xlabel('f/Hz');ylabel('|Y(f)|')
        title('Double-Sided Amplitude Spectrum of Y(t)')
        if nargin == 6
           title(TITLE)
        elseif nargin == 8
           xlabel(XLABEL);ylabel(YLABEL);title(TITLE)
        end
    elseif strcmp(GRAPH,'stem')
        figure
        stem(f,Y)
        xlabel('f/Hz');ylabel('|Y(f)|')
        title('Double-Sided Amplitude Spectrum of Y(t)')
        if nargin == 6
           title(TITLE)
        elseif nargin == 8
           xlabel(XLABEL);ylabel(YLABEL);title(TITLE)
        end
    else
        error('GRAPH parameter must be "none", "plot" or "stem"!!!')
    end
else
    error('SIDES parameter must be "single" or "double"!!!')
end
