%First we extract the data from the event file: Event Time Stamps & Event Strings.
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
%This loop allows you to extract data from between two time points. Good
%and meticulous commenting of your recordings will help make this part
%painless, quick, and dependable. 

% You need to search through the strings manually and find what are the
% time points you want to extract between. It is faster than an automated
% process, since you don't waste processing time on unneccessary data. 
%% 
% Add a user input for the ranges (have them click on what are the strings
% that they want to extract from, and the those strings will be loaded into
% a 2x(length(user input)/2) Matrix, and each column will be a range. 

%https://www.mathworks.com/help/matlab/ref/uitable.html
%%
%This version puts all electrode voltage values in one matrix with 64 rows
for i = 1:64

FileName=[BaseNameCSC,num2str(i),ExtensionnameCSC]

%%% This is the range of records that will be extraced when using ExtractMode = 4
RangeOfTime = [eventTimeStamps(6) eventTimeStamps(7)];

%%% This is the total recording time (in seconds, thus 1e-6 term) 
SpanOfRecording=(eventTimeStamps(7)-eventTimeStamps(6))*1e-6

%Nlx2MatCSC extracts time stamps and unconverted sample values
[TimeStampsForSamples, RawVoltageSamples, Header] = Nlx2MatCSC( FileName, FieldSelection, 1, ExtractMode, RangeOfTime);

%%% Number of samples 
NumberofSamples=length(RawVoltageSamples)

%%%%%%%%%%%%%%%%%%% Now we will generate the time and voltage vectors, and
%%%%%%%%%%%%%%%%%%% combine then into one matrix, for each unique
%%%%%%%%%%%%%%%%%%% experiment.






%Mean(Sample Value) * ADBitVolts (value in header) * V to mV conversion
rawvoltage(:,i) = mean(RawVoltageSamples,1)*   3.0519e-09 * 1e6; 
%NBprestim_time=NB_Prestim_TimeStamps/1e6; You can use Timestamp data as
%well, but it won't start at time = 0 seconds. 
time(:,i)=linspace(0,SpanOfRecording, NumberofSamples);

%Place Holder Time Voltage: 
% Placeho...= {'Name your time', time, 'name your stim',rawvoltage} 
PlaceholderTimeVoltage = {'NB_Prestim',time;'NB_Baseline',rawvoltage}
ArrayVoltageTime = PlaceholderTimeVoltage';

%Name this whatever you like
NB_Agg = struct(ArrayVoltageTime{:});

%Name this whatever you like
NB_Aggregate(1).DIV14(1,1) = NB_Agg

save('NameofFiletoSave','NB_Aggregate')

end
%This is a test of gitsssaa d   d
%% Averaging Data
FileLocation = uigetdir(path,'Select the folder where the data is ~.^')

ExtractDataFunction(FileLocation)
%% Plot Raw Data
%Ideally you want to plot a single graph, because if you do all 64 its
%going to look bad. 
electrode=[1:64] 
for i=1:64

figure (1) 
subplot (8,8,i)
plot(time,rawvoltage(i) )
 sgtitle('Raw Voltage Signal') 
 %title(['\Electrode Number = ', num2str(electrode(i))]);
end
%%
tms{1,:}=(NameofExperiment(23).PreStim_Time)
vlt{1,:}=(NameofExperiment(23).StimData)
spike=detectspike(vlt{1}',tms{1}')

