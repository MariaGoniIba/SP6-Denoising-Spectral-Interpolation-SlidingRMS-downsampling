clear; clc; close all

load('resample_codeChallenge.mat')

figure(1)
plot(time, signal, 'k-', 'linew',3)
title('Proposal exercise','FontSize',14)

%% Dealing with NaNs
% They are interspearsed (not whole chunk of them) so I interpolate (spline)

s = signal;

nans = isnan(signal);
t    = 1:numel(signal);
s(nans) = interp1(t(~nans), s(~nans), t(nans), 'linear');

% Original signal with interpolated
figure(2)
plot(signal,'ks-', 'linew',3)
hold on
plot(s+0.5,'r','linew',3) % offset for better visualization with the solution
legend({'With NaNs';'Without NaNs'},'FontSize',12)
title('Dealing with NaNs','FontSize',14)


%% detect bad segments using sliding RMS

% window size as percent of total signal length
pct_win = 1.1; % in percent, not proportion!

% convert to indices
k = round(size(s,1) * pct_win/2/100);


% initialize RMS time series vector
rms_ts = zeros(1,size(s,1));

for ti=1:size(s,1)
    
    % boundary points
    low_bnd = max(1,ti-k);
    upp_bnd = min(size(s,1),ti+k);
    
    % signal segment (and mean-center!)
    tmpsig = s(low_bnd:upp_bnd);
    tmpsig = tmpsig - mean(tmpsig);
    
    % compute RMS in this window
    rms_ts(ti) = sqrt(sum( tmpsig.^2 ));
end

% plot RMS
figure(3), clf, hold on
plot(rms_ts,'s-')
title('Detecting bad segments with sliding RMS','FontSize',14)

% pick threshold manually based on visual inspection
thresh = 120;
plot(get(gca,'xlim'),[1 1]*thresh,'r')
legend({'Local RMS';'Threshold'})

signalR = s;
signalR( rms_ts>thresh ) = NaN;

figure(4), clf, hold on
plot(s,'ks-', 'linew',3)
plot(signalR,'r','linew',3)
legend({'Original';'Thresholded'},'FontSize',12)

%% Spectral interpolation
% Detect beginning and end of window
indbeg = []; indend = [];
for i=2: size(signalR,1)-1
    if isnan(signalR(i)) && ~isnan(signalR(i-1))
        indbeg = [indbeg i];
    elseif isnan(signalR(i)) && ~isnan(signalR(i+1))
        indend = [indend i];
    end
end

filtsig = signalR;
% I assume each empty window has a beginning and an end
for i = 1: length(indbeg)
    % FFTs of pre- and post-window data
    fftPre = fft(signalR( indbeg(i)-diff([indbeg(i) indend(i)])-1:indbeg(i)-1) );
    fftPst = fft(signalR( indend(i)+1:indend(i)+diff([indbeg(i) indend(i)])+1) );

    % interpolated signal is a combination of mixed FFTs and straight line
    mixeddata = detrend(ifft( ( fftPre+fftPst )/2 ));

    % formula for a line: x*m +b 
    % m (slope: 1st valid point after the missing data chunk - first valid point before the missing point) 
    % b (offset: first point)
    linedata  = linspace(0,1,diff([indbeg(i) indend(i)])+1)'*(signalR(indend(i)+1)-signalR(indbeg(i)-1)) + signalR(indbeg(i)-1);

    % sum together for final result
    linterp = mixeddata + linedata;

    % put the interpolated piece into the signal
    %filtsig = signal;
    filtsig(indbeg(i):indend(i)) = linterp;
end

figure(5), clf
plot(1:length(signalR),[signal signalR filtsig+8],'linew',2)
legend({'Original';'With gap';'Filtered (+ offset)'},'FontSize',12)
title('Spectral interpolation','FontSize',14)
zoom on

%% Resample to regularly spaced

% With the derivative of the time vector we see that data is not regularly spaced data
% plot(diff(time)*1000) 
% To interpolate, we choose the highest frequency (minimum time)
new_fs = 1/min(diff(time));
npnts = new_fs*time(end); % 2 seconds
new_time  = (0:npnts-1)/new_fs;

F = griddedInterpolant(time,filtsig','linear');
interp_s = F(new_time);

figure(6)
plot(time, filtsig, 'k', 'linew',3)
hold on
plot(new_time, interp_s+0.3,'r','linew',3)
legend({'Irregular srate';'Regular srate+offset'},'FontSize',12)
title('Regularly spaced','FontSize',14)