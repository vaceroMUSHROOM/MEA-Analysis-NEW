%% Spike Detection draft and notes 
MinSpikeTimeInterval = 2;%In # of samples, not time 
%Break recording into chunks, and calculate sigma for each chunk. The
%calculated value for sigma, for a very long recording, might vary, so its
%better to calculate "local" sigmas multiple times. After finding the
%spikes, these times will all be joined together. 
Chunk = length(DownsampledData(23,1:end-8))/10;
for i=1:10
    if i==1
Sigma(1,i) = median(abs(DownsampledData(23,1:(Chunk*i))-median(DownsampledData(23,1:(Chunk*i)))))/0.6745;
    else Sigma(1,i) = median(abs(DownsampledData(23,Chunk*(i-1):(Chunk*i))-median(DownsampledData(23,Chunk*(i-1):(Chunk*i)))))/0.6745;
    end
    end
DetectionThreshold = 5.*Sigma;
DetectionThresholdMatrix = reshape(repmat(DetectionThreshold,Chunk,1),1,[]);

for i=1:10
[locs pks]=peakseek(DownsampledData(23,:),MinSpikeTimeInterval,DetectionThreshold);
locations(i,:)=locs;
peaks(i,:)=pks;
end

figure (1)
plot(DownsampledData(23,1:3730739))
figure (2)
plot(DownsampledDataband(23,1:3730739))
figure(3)
plot(DownsampledDatalow(23,1:3730739))

%% Spike Detection via Thresholding (RMS)
RMS=rms(DownsampledData(23,32000:3072000));
DetectionThresholdMatrixRMS = reshape(repmat(RMS,(3072000-32000),1),1,[])*3;
[locs pks]=peakseek(DownsampledData(23,32000:3072000),MinSpikeTimeInterval,DetectionThreshold);


%% Spike Detection via Thresholding (MAD)

Sigma = median(abs(DownsampledData(23,32000:3072000)-median(DownsampledData(23,32000:3072000))))/0.6745;
MinSpikeTimeInterval = 1;%In # of samples, not time 
DetectionThreshold = 6.*Sigma;
DetectionThresholdMatrix = reshape(repmat(Sigma,(3072000-32000),1),1,[])*6;
[locs pks]=peakseek(DownsampledData(23,64000:3072000),MinSpikeTimeInterval,DetectionThreshold);
[nlocs npks] = peakseek(DownsampledData(23,64000:3072000)*-1,MinSpikeTimeInterval,DetectionThreshold);
%Comvine locs (positive) and nlocs (negative) to have the full spike time
%list

%% Enhance spike detection using slope analysis 

%% Binary Spike PLot and Firing Rate
binaryspike = zeros(1,3072000-64000); %Fill the binaryspike matrix with zeros for the length of recording
binaryspike(1,locs) = 1 ; %Place 1 where there is a spike (locs = positive threshold)
binaryspike(1,nlocs) = 1; %Place 1 where there is a spike (nlocs = negative threshold) 
%20 ms histogram bins 
%128 samples in a 20 ms time bin 
%Rangeoftime * dt = seconds total 
plot (binaryspike)

figure 
% 20 ms time bin firing rate 
N= hist(RelativeTimeNew(locs), 23500);
bar(N)
ylabel('Firing Rate (Spikes/20ms Bin)')
xlabel('Time (s)')
title('Firing Rate for Agg in NB (DIV 14)')
VarianceofFR=sqrt(mean((N-mean(N)).^2)); % Standard Deviation of FR 20 ms time bins 
%FR with Guassian Kernel 
r=-10*0.02:1/2000:10*.02;
y=normpdf(r,0,0.02)/2000;%Normalize gaussian so that integral = 1
%gaussmf(X, [SIGMA, C]) = EXP(-(X - C).^2/(2*SIGMA^2));
instantFR=conv(binaryspike,y,'same');



%% Graphs and Figures (8 minute, 1 minute, 1 second)
%8 minute trace of Aggregates in Neurobasal, DIV 14
figure (1)
plot(RelativeTimeNew(1,64000:3072000),DownsampledData(23,64000:3072000))

%Zoomed in 1 minute trace of above 
figure (2)
plot(RelativeTimeNew(1,64000:448000),DownsampledData(23,64000:448000))


%Zoomed in 1 second trace of above 
figure (3)
plot(RelativeTimeNew(1,204800:211200),DownsampledData(23,204800:211200))

%% Spike Train Analysis (FFT, Powerspectrum)
%Find the oscillations 
%Concern about xcorr** Should I do 'coeff'?
DataFFT = fftshift(abs(fft(DownsampledData(23,:))));
fshiftDATA = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershiftDATA = abs(DataFFT).^2/n;     % zero-centered power
nData=length(powershiftDATA)
nSpike=length(powershiftSPIKE)
figure (1)
semilogy(fshiftDATA,powershiftDATA)
figure (2)
semilogy(fshiftSPIKE,powershiftSPIKE)
figure (2)
%LogLog of FFT with every 30th data point taken (decimation) in order to
%smooth out the data and make it easier to determine what the dom frq. is
loglog(fshiftDATA(nData/2:30:nData),powershiftDATA(nData/2:30:nData))
%LogLog of FFT with every 10th data point taken, followed by moving average
%with a 10 sample window.
loglog(smooth(fshiftDATA(nData/2:10:nData),10),smooth(powershiftDATA(nData/2:10:nData),10))
%LogLog of FFT with every 10th data point taken, followed by moving average
%with a 20 sample window. ***THIS I THINK IS IDEAL FFT***
loglog(smooth(fshiftDATA(nData/2:10:nData),10),smooth(powershiftDATA(nData/2:10:nData),10))
%
%
%
% 
%Now I'm going to do the spectogram 


%% Find oscillations 
Autocorrelation=xcorr(binaryspike,binaryspike); %
plot(Autocorrelation);
MovingAve= smooth(Autocorrelation,'moving');
Fit= smooth(Autocorrelation, 'sgolay',4);

[S,F,T] = spectrogram(DownsampledData(23,:),25600,1250,25600,6400);

S = abs(S);
imagesc(T(1,1:100),F(1:100),10*log10(S(1:100,:)/max(S(1:100,:))))
[S,F,T] = spectrogram(DownsampledData(23,:),25600,25600/2,25600/4,6400);
S=abs(S);
plot(Fit(n/2:n))
plot(MovingAve(n/2:n))
%Since we don't care for any signal above Hz (we don't expect any real
%oscillations in this range), we will smooth the autocorrelation signal
%using the smooth() function. This function will essentially smooth the
%signal with a moving average (5) 
length(Autocorrelation)
figure
plot(Autocorrelation)
figure (1)
as=periodogram(binaryspike,'power')
plot(as)
figure (1)
plot(as(1:1000))
%Now I want to try to find the 
ppd=fftshift(abs(fft(Autocorrelation)));
plot(ppd(1:1000))
%n=number of samples/datapoints pwelsh pmtm
n=length(Autocorrelation)
fs=6400
fshift_ppd = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift_ppd = abs(ppd).^2/n; 
loglog(fshift_ppd,powershift_ppd)
this=pmtm(binaryspike);
length(this)
ppdW=pwelch(binaryspike,'power');
figure (2)
plot(ppdW(1:1000))
plot(ppdWthousand)% Very wide, frequencies go much higher (100-600)
figure (3)
plot(ppdWten) %very narrow, frequency doesn't go above 130 
%What is the ideal window size? 
%How does pmtm compare to pwelch???? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
plot(this(1:1000))

 

%% Burst Analysis FileName