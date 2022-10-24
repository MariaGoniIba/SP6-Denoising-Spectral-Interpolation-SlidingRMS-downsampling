# Signal Processing (SP) project 6: Denoising, Spectral Interpolation, Sliding RMS and downsampling
Spotting different issues of a signal and cleaning/reconstructing it.

# Credit
The dataset and proposal of the exercise is from the Udemy course [Signal Processing Problems, solved in Matlab and in Python](https://www.udemy.com/course/signal-processing/). I highly recommend this course for those interested in signal processing.

# Procedure
Let's plot the signal to spot the issues.

<p align="center">
    <img width="800" src="https://github.com/MariaGoniIba/SP6-Denoising-Spectral-Interpolation-SlidingRMS-downsampling/blob/main/Figure1-ProposalExercise.png">
</p>

I can see:
* A lot of missing data. A look on the data values show us a lot of NaN values.

* Two clear noise bursts.

* When plotting the derivative of the time vector, it appears that this data is not regularly spaced.

Let's do something about it!

## Dealing with NaNs values
Given that these NaNs values appear to be interpeased (not a whole chunk of them), I apply spline interpolation.
```
s = signal;

nans = isnan(signal);
t    = 1:numel(signal);
s(nans) = interp1(t(~nans), s(~nans), t(nans), 'linear');
```
<p align="center">
    <img width="800" src="https://github.com/MariaGoniIba/SP6-Denoising-Spectral-Interpolation-SlidingRMS-downsampling/blob/main/Figure2-NaNs.png">
</p>

## Bursts of noise

First I will identify the entire noisy window using sliding RMS (root mean square). 
This technique detects whole windows of outliers or noise, instead of single points. 

I specify a window based on the percent of the total signal. 
In a loop, I mean center that part of the signal. This way you focus on the variance of the signal.

```
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
```

When plotting it, it seems clear that this method has succesfully identified these 2 time windows where there are extreme data values that we want to get rid of.
I pick up the threshold (_thresh = 120_) manually based on visual inspection.

<p align="center">
    <img width="800" src="https://github.com/MariaGoniIba/SP6-Denoising-Spectral-Interpolation-SlidingRMS-downsampling/blob/main/Figure3-SlidingRMS.png">
</p>

I equal all values above that threshold to NaN.
```
signalR = s;
signalR( rms_ts>thresh ) = NaN;
```
<p align="center">
    <img width="800" src="https://github.com/MariaGoniIba/SP6-Denoising-Spectral-Interpolation-SlidingRMS-downsampling/blob/main/Figure4-DetectionBadSegments.png">
</p>

Finally, I apply spectral interpolation to fill the gaps. 
Since this consists of missed chunk of data unlike random missed data, I use spectral interpolation. 
The idea is that we make the assumption that the power spectrum of the missed window smoothly transitions from whatever the spectrum is before to whatever the spectrum is after.
Therefore, we take 2 Fourier transforms, one before the boundary of the missing window and another one after the boundary.
Then we average these 2 spectrum together. This gives us a third spectrum and we take the inverse FFT to get a time domain signal.
This time domain is probably not going to transition from boundaries. 
Then, we add a linear trend that will smoothly connect both points.
For better visualization of the difference between the non-filtered and filtered data, I add an offset.

```
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
```
<p align="center">
    <img width="800" src="https://github.com/MariaGoniIba/SP6-Denoising-Spectral-Interpolation-SlidingRMS-downsampling/blob/main/Figure5-SpectralInterpolation.png">
</p>

## Irregularly spaced
To resample to a regular space, I interpolate to the highest frequency (minimum time).
```
new_fs = 1/min(diff(time));
npnts = new_fs*time(end); % 2 seconds
new_time  = (0:npnts-1)/new_fs;

F = griddedInterpolant(time,filtsig','linear');
interp_s = F(new_time);
```

Finally, this is how the new signal looks like!
<p align="center">
    <img width="800" src="https://github.com/MariaGoniIba/SP6-Denoising-Spectral-Interpolation-SlidingRMS-downsampling/blob/main/Figure6-Solution.png">
</p>
