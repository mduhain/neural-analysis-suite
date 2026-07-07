%im2pconv.m
%
% im2pconv() converts .im2p files to .nwb files
% @mduhain 2023-07-24
%
% good test folder locally on the analysis machine is m748 / 2023-03-22 / Exp001 / 
% 
%% Importing Data

% hardcoded temp variables (EDIT AS NEEDED)
tFileNum = 1; %if multiple im2p files exist in this folder
tMouseID = 'm810';
tDate = '2023-04-20';
tExp = 'Exp001';
tSessID = strcat(tMouseID,'-',num2str(tFileNum));
tMouseSex = 'F';
tMouseGenotype = 'PV-Cre-tdTom';
tMouseStrain = 'C57BL/6';
tMouseVirus = 'AAV1-hSyn-jGCaMP7c-WPRE';
tFluorophore = 'GCaMP7c';
tRawFramerate = 30.05;

% find and load .im2p-analyzed file
cd(strcat('C:\Users\ramirezlab\Desktop\ExRawData\',tMouseID,'\',tDate,'\',tExp));
[depths, fileNames] = analyzeDir();
if any(contains(fileNames,'im2p')) && any(contains(fileNames,'analyzed'))
    %Load in im2p file for Exp Metadata
    targetList = fileNames(contains(fileNames,'im2p-analyzed'));
    disp(strcat("Loading ",targetList(tFileNum),"..."));
    load(targetList(tFileNum));
    a1 = a;
    clear a
end


% Find and load raw images contained in 'im2p' non-analyzed version, and
% the rfp traces.
% this is the quick way, but the data is down-sampled 75% space and 66% in time.
if any(contains(fileNames,'im2p')) && any(~contains(fileNames,'analyzed'))
    targetList = fileNames(contains(fileNames,'im2p') & ~contains(fileNames,'analyzed'));
    disp(strcat("Loading ",targetList(tFileNum),"..."));
    load(targetList(tFileNum)); %a contains raw images
end



%% Creating and Assign NWB
tNumTrials = length(a1.stimOnTimes);
tImgDepth = depths(tFileNum);

% Create a blank NWB instance
nwb = NwbFile( ...
    'session_description', 'Mouse forepaw vibrotactile stimulation',...
    'identifier', strcat(tMouseID,'-',tSessID), ...
    'session_start_time', datetime(2022, 03, 22, 12, 59, 59), ...
    'general_experimenter', 'mduhain', ... % optional
    'general_session_id', 'session_0002', ... % optional
    'general_experiment_description',['2P recording of vibrotactile ' ...
    'induced responses in forepaw somatosensory cortex'],...
    'general_institution', 'University of Rochester', ...
    'general_virus', tMouseVirus );

% TODO: Calculate mouse age based on birthday?
% TODO: Extract mouse genotype
% TODO: Figure out how I can add viral information to subject

% Add subject information
subject = types.core.Subject( ...
    'subject_id', tMouseID, ...
    'age', 'P100D', ...
    'date_of_birth', datetime(2021, 11, 29, 12, 59, 59), ...
    'species', 'Mus musculus', ...
    'sex', tMouseSex, ...
    'genotype', tMouseGenotype,...
    'strain',tMouseStrain);
nwb.general_subject = subject;

%---------------------------------------------------------------------------
%Check to see if the stim times are properly loaded
disp("Checking if stimulus times are correctly loaded...");
errorFlag = 0;
if any(isnan(a1.stimOnTimes))
    disp("FAIL: NaN Values found in a.stimOnTimes"); 
    errorFlag = 1;
end
if any(a1.stimOnFrameNum == 0)
    disp("FAIL: Values of 0 found in a.stimOnTimes");
    errorFlag = 1;
end
if any(isnan(a1.stimOnFrameNum))
    disp("FAIL: NaN Values found in a.stimOnFrameNum");
    errorFlag = 1;
end
if any(a1.stimOnFrameNum == 0)
    disp("FAIL: Values of 0 found in a.stimOnFrameNum");
    errorFlag = 1;
end


numStimTimes = length(find((a1.stimOnTimes ~= 0)));
numStimFrams = length(find((a1.stimOnFrameNum ~= 0)));

%find stim times (from rfpTrace)
rfpExpTrace = detrend(a.rfpExpTrace);
nFrames = length(rfpExpTrace); %should match the number of files in Z000
rfpTraceMean = mean(rfpExpTrace);
rfpTraceSTD = std(rfpExpTrace);
rfpTraceThreshold = rfpTraceMean + 2*rfpTraceSTD;
%plotting section to visualize RFP spikes
figure; plot(rfpExpTrace)
hold on;
plot([1 length(rfpExpTrace)],repmat(rfpTraceThreshold,2,1),'r-');
%calculate RFP event times
I = find(rfpExpTrace > rfpTraceThreshold);
dI = diff(I);
dI = [0; dI];
I(dI == 1) = [];
I(1:3) = []; %remove first three flashes (exp start id)
disp(strcat(num2str(length(I))," trials identified from RFP Trace."));
%more calculations of other relevent ata
nRealTrials = length(I);
stimOnTimes = I ./ tRawFramerate;
interStimInterval = diff(stimOnTimes);
trialStartTimes = zeros(nRealTrials,1);
trialStopTimes = zeros(nRealTrials,1);
for n = 1 : length(interStimInterval)
    trialStartTimes(n+1) = sum(interStimInterval(1:n));
    trialStopTimes(n) = sum(interStimInterval(1:n))-0.001;
end
frameTimes = 0:(1/tRawFramerate):((nFrames/tRawFramerate)-(1/tRawFramerate));
frameTimes = frameTimes';
trialStopTimes(end) = nFrames/tRawFramerate;

%---------------------------------------------------------------------------
% Add trial information
trials = types.core.TimeIntervals( ...
    'colnames', {'start_time', 'stop_time','stim_on_time','stim_duration', 'stim_frequency', 'stim_amplitude', 'audio-only_trial'}, ...
    'description', 'trial data and properties', ...
    'id', types.hdmf_common.ElementIdentifiers('data', (0:(nRealTrials-1))'), ...
    'start_time', types.hdmf_common.VectorData( ...
        'data', trialStartTimes, ...
   	    'description','start time of trial in seconds' ...
    ), ...
    'stop_time', types.hdmf_common.VectorData( ...
        'data', trialStopTimes, ...
   	    'description','end of each trial in seconds' ...
    ), ...
    'stim_on_time', types.hdmf_common.VectorData( ...
        'data', stimOnTimes, ...
   	    'description', 'stimulus ON time in seconds' ...
    ), ...
    'stim_duration', types.hdmf_common.VectorData( ...
        'data', a1.stimDurs(1:nRealTrials), ...
   	    'description', 'stimulus duration in milliseconds' ...
    ), ...
    'stim_frequency', types.hdmf_common.VectorData( ...
        'data', a1.stimFreqs(1:nRealTrials), ...
   	    'description', 'stimulus frequency in hertz' ...
    ), ...
    'stim_amplitude', types.hdmf_common.VectorData( ...
        'data', a1.stimAmps(1:nRealTrials), ...
   	    'description', 'stimulus amplitude in volts' ...
    ), ...
    'audio-only_trial', types.hdmf_common.VectorData( ...
        'data', a1.stimAmps(1:nRealTrials), ...
   	    'description', 'a trail containing no tactile stimulus, only the auditory component of the stimulus, in a boolean value') ...
    );

nwb.intervals_trials = trials;


%---------------------------------------------------------------------------
% Add static image(s) of RFP field, Vasculature, Overlay

% grayscale vasculature image 1x
if any(strcmp(fileNames,'vasc1.png'))
    image_data1 = double(imread('vasc1.png'))';
elseif any(strcmp(fileNames,'Vasc1.png'))
    image_data1 = double(imread('Vasc1.png'))';
end
grayscale_image1 = types.core.GrayscaleImage( ...
    'data', image_data1, ...  % required
    'resolution', 70.0, ... %resolution in pixels per centimeter
    'description', 'Grayscale image of session vasculature at 1x digital zoom' ...
);

% grayscale vasculature image 3x
if any(strcmp(fileNames,'vasc3.png'))
    image_data3 = double(imread('vasc3.png'))';
elseif any(strcmp(fileNames,'Vasc3.png'))
    image_data3 = double(imread('Vasc3.png'))';
end
grayscale_image2 = types.core.GrayscaleImage( ...
    'data', image_data3, ...  % required
    'resolution', 70.0, ... % resolution in pixels per centimeter
    'description', 'Grayscale image of session vasculature at 3x digital zoom' ...
);

%add into image container
image_collection = types.core.Images( ...
    'description', 'A collection of logo images presented to the subject.'...
);

image_collection.image.set('grayscale_image1', grayscale_image1);
image_collection.image.set('grayscale_image2', grayscale_image2);
 
nwb.acquisition.set('image_collection', image_collection);


%---------------------------------------------------------------------------
% Store the animals bar grasping and licking data
grasp_intervals = types.core.TimeIntervals( ...
    'description', 'Intervals when the animal was grasping the bar.', ...
    'colnames', {'start_time', 'stop_time'} ...
);
for n = 1 : length(a1.graspTimes)
    grasp_intervals.addRow('start_time', a1.graspTimes(n,1), 'stop_time', a1.graspTimes(n,2));
end

nwb.intervals.set('grasp_intervals', grasp_intervals);


%---------------------------------------------------------------------------
% Add imaging plane for 2P session
optical_channel = types.core.OpticalChannel( ...
    'description', 'description', ...
    'emission_lambda', 510);

device = types.core.Device();
nwb.general_devices.set('Device', device);

imaging_plane_name = strcat('imaging_plane_',tImgDepth);
imaging_plane = types.core.ImagingPlane( ...
    'optical_channel', optical_channel, ...
    'description', 'primary somatosensory cortex, forepaw', ...
    'device', types.untyped.SoftLink(device), ...
    'excitation_lambda', 870, ...
    'imaging_rate', 10., ...
    'indicator', tFluorophore, ...
    'location', 'primary somatosensory cortex, forepaw');

nwb.general_optophysiology.set(imaging_plane_name, imaging_plane);


%---------------------------------------------------------------------------
%The correct way is to use DataPipe an iterative writing
cd(strcat('C:\Users\ramirezlab\Desktop\ExRawData\',tMouseID,'\',tDate,'\',tExp,'\',depths(tFileNum)));
[~,frameNames]=analyzeDir();
numFrames = length(frameNames);
testImage = imread(frameNames(1),1); 
[pxSizeY, pxSizeX] = size(testImage);

% load in 1/divSiz of the image stack
divSize = 10; %the number of sections to divide image stack into
stepSize = floor(numFrames/divSize); %how many frames are contained in 1/divSize
pPool = parpool;
dataPart = zeros(pxSizeY,pxSizeX,stepSize);
disp("Loading first section of images..."); tic;
parfor n = 1 : stepSize
    dataPart(:,:,n) = imread(frameNames(n),1); % 1 because GFP is in channel 1
end
toc;

fullDataSize = [pxSizeY, pxSizeX, numFrames]; % this is the size of the TOTAL dataset

% compress the data
fData_use = types.untyped.DataPipe( ...
    'data', dataPart, ...
    'maxSize', fullDataSize, ...
    'axis', 3);

% set compressed data as a twoPhotonSeries
fdataNWB = types.core.TwoPhotonSeries( ...
    'imaging_plane', types.untyped.SoftLink(imaging_plane), ...
    'starting_time', 0.0, ...
    'starting_time_rate', 10, ...
    'data', fData_use, ...
    'data_unit', 'lumens');


nwb.acquisition.set('TwoPhotonSeries', fdataNWB);

cd(strcat('C:\Users\ramirezlab\Desktop\ExRawData\',tMouseID,'\',tDate,'\',tExp));
uniqueFileName = strcat('tactileStim_',tSessID,'.nwb');
nwbExport(nwb,uniqueFileName);

nwb = nwbRead(uniqueFileName); %load the nwb file with partial data

% "load" each of the remaining 1/divSizes of the image stack
for i = 2:divSize % iterating through parts of data
    cd(strcat('C:\Users\ramirezlab\Desktop\ExRawData\',tMouseID,'\',tDate,'\',tExp,'\',depths(tFileNum)));
    if i == divSize
        finalStep = ((i-1)*stepSize+1):numFrames;
        finalStepSize = length(finalStep);
        dataPart = zeros(pxSizeY,pxSizeX,finalStepSize);
        disp(strcat("Loading ",num2str(i),"/",num2str(divSize)," section of images...")); tic;
        parfor n = 1 : finalStepSize
            dataPart(:,:,n) = imread(frameNames(finalStep(n)),1); % 1 because GFP is in channel 1
        end
        toc;
    else 
        dataPart = zeros(pxSizeY,pxSizeX,stepSize);
        disp(strcat("Loading ",num2str(i),"/",num2str(divSize)," section of images...")); tic;
        parfor n = 1 : stepSize
            dataPart(:,:,n) = imread(frameNames(i*stepSize-(stepSize-n)),1); % 1 because GFP is in channel 1
        end
        toc;
    end
    cd(strcat('C:\Users\ramirezlab\Desktop\ExRawData\',tMouseID,'\',tDate,'\',tExp));
    disp("Appending next iteration of images..."); tic;
    nwb.acquisition.get('TwoPhotonSeries').data.append(dataPart); % append the loaded data
    toc;
end

%---------------------------------------------------------------------------





%% Saving



%% Exporting

%done above




