FileLocation = uigetdir(path,'Select the folder where the data is ~.^')

ExtractDataFunction(FileLocation)

%SpikeDetection(NameOfExperiment,StandardDeviation)


%%

FileLocation = uigetdir(path,'Select the folder where the data is ~.^')
BaseName=[FileLocation];
BaseNameCSC=[BaseName,'\CSC'];
BaseNameEvent=[BaseName,'\Events.nev']
ExtensionnameCSC='.ncs'
FieldSelection=[1 0 0 0 1];
[eventTimeStamps, eventStrings, Header] = Nlx2MatEV( BaseNameEvent, FieldSelection, 1, 1, 1 );
FieldSelection=[1 0 0 0 1];
ExtractMode = 4;

%%
VoltageData=[]
[c,d] = butter(6, 2000/(32000/2),'low'); %0 - 2000 Hz filter (butterworth)
[clow,dlow]=butter(6,100/32000/2,'low');
[cband,dband]=butter(2,[150/32000/2 1500/32000/2],'stop');
SamplesPerAverage = 5;             % Number of elements to create the mean over



%% The other code, e.g. ExtractDataFunction and MasterCode, are trying to do this but for every single electrode. 

%This initial code will give you the filtered (2K Hz low bandpass), and
%downsample to a 1/5 the # o samples (~6,400)

%This data is for spike analysis 
for i = 1:23

FileName=[BaseNameCSC,num2str(i),ExtensionnameCSC];

%%% This is the range of records that will be extraced when using ExtractMode = 4
RangeOfTime = [eventTimeStamps(6) eventTimeStamps(7)];

%%% This is the total recording time (in seconds, thus 1e-6 term) 
SpanOfRecording=(eventTimeStamps(7)-eventTimeStamps(6))*1e-6;

%Nlx2MatCSC extracts time stamps and unconverted sample values
[TimeStampsForSamples, RawVoltageSamples, Header] = Nlx2MatCSC( FileName, FieldSelection, 1, ExtractMode, RangeOfTime);


ADConv = str2num(Header{16,1}(end-25:end)); %Pulls out the conversion to obtain voltage

VoltageData = RawVoltageSamples*ADConv*1e6;% convert matrix to row vector

ReshapedData(1,:)=reshape(VoltageData(:,:),1,[]);

filtered_data = filter(c,d,ReshapedData);

if i==1 %Here we will determine the dimension of and build the Average matrix. 

size = length(filtered_data);      % Find the next smaller multiple of n

m  = size - mod(size, SamplesPerAverage);

DownsampledData=zeros(64,m/SamplesPerAverage);  %Matrix with 64 electrodes.

end

DataReshapedForAveraging  = reshape(filtered_data(1:m), SamplesPerAverage, []);     % Reshape x to a [n, m/n] matrix
%1st 3rd 5th 7th .... | This is what the matrix would look like if I wanted
%2nd 4th 6th 8th ....| to downsample by 2 
DownsampledData(i,:) = sum(DataReshapedForAveraging, 1) / SamplesPerAverage;  % Calculate the mean over the 1st dim

end
%%

for i=1:23
    
%This will yield the filtered data for LFP detection

filtered_data_LFP = filter(clow,dlow,ReshapedData);

if i==1 %Here we will determine the dimension of and build the Average matrix. 

size = length(filtered_data_LFP);      % Find the next smaller multiple of n

m  = size - mod(size, SamplesPerAverage);

DownsampledData_LFP=zeros(64,m/SamplesPerAverage);  %Matrix with 64 electrodes.

end

DataReshapedForAveraging_LFP  = reshape(filtered_data_LFP(1:m), SamplesPerAverage, []);     % Reshape x to a [n, m/n] matrix
%1st 3rd 5th 7th .... | This is what the matrix would look like if I wanted
%2nd 4th 6th 8th ....| to downsample by 2 
DownsampledData_LFP(i,:) = sum(DataReshapedForAveraging_LFP, 1) / SamplesPerAverage;  % Calculate the mean over the 1st dim
end
%% Get correct time stamps 
tstampsData=TimeStampsForSamples;
data=RawVoltageSamples;
   [ts,~] = getCorrectTStamps(tstampsData,data);

%If you are averaging your data, then use the following code to get new
%time stamps.
TimeReshapedForAveraging  = reshape(ts(1:m), SamplesPerAverage, []);     % Reshape x to a [n, m/n] matrix
ExactTimeNew = 1e-6*sum(TimeReshapedForAveraging, 1) / SamplesPerAverage; 
RelativeTimeNew= ExactTimeNew - ExactTimeNew(1);
dt=RelativeTimeNew(1,2);
%% Check if your filter worked 
RawFFT = abs(fft(DownsampledData(23,:)));
FilteredFFT=abs(fft(filtered_data(23,:)));
n=length(DownsampledData(23,:))
Y = fftshift(RawFFT);
fshift = (-n/2:n/2-1)*(fs/n); % zero-centered frequency range
powershift = abs(Y).^2/n;     % zero-centered power
figure (2)
semilogy(fshift,powershift)
%plot(fshift,powershift) this is the plot of the fft centered around zero,
%without a log scale. its hard to interpret the figure without a log scale.
