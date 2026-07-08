% im2p.m
% 2p image class
% 
% <2022-08-15> mduhain
%   - created fucntion and defined components

classdef im2p
    %======================================================================
    properties
    % RFP COMPONENTS    
    rfpPVImgPre = NaN; %rfp image of PV neurons in field pre-experiment.
    rfpPVImgPost = NaN; %rfp img of PV neurons in field post-experiment.
    rfpExpTrace = NaN; %raw RFP trace during exp to track stim times.
    overlayImg = NaN;
    
    % GFP COMPONENTS
    gfp3FrameAvg = NaN; %3 frame averaged gfp frames from experiment.
    gfpMotionCorrected = NaN; %A motion corrected version of gfp3FrameAvg
    templateImg4MotionCorrection = NaN; %template Image from Importer
    gfpXCorrImg = NaN; %Pixel correlated image from RoiMaker.XCorr();

    % FLUORESCENCE TRACES
    somaticF = NaN; %extracted somatic fluorescence traces;
    somaticFDT = NaN; %Somatic fluorescence traces de-trended (photobleach)
    redSomaticF = NaN; %extracted somatic fluorescence traces from RFP labeled neurons;
    redSomaticFDT = NaN; %Somatic fluorescence traces de-trended (photobleach) from RFP
    neuropilF = NaN; %extraced neuropil somatic surround (paired with somaticF);
    neuropilFDT = NaN; %extraced neuropil somatic surround (paired with somaticF) de-trended;
    npsubSomaticF = NaN; %neuropil subtracted somatic F;
    npsubSomaticFDT = NaN; %neuropil subtracted somatic FDT;
    
    % META DATA COMPONENTS
    digZoomLevel = NaN; %Value for digital zoom level of experiment
    expMetaData = NaN; %table to hold experiment meta data
    numFrames = NaN; %number of frames in experiment
    stimOnTimes = NaN; %tactile stimulus ON times
    stimOnFrameNum = NaN; %frames cooresponding to stim ON
    ledTimes = NaN; %clock time when rig LED was turned ON
    stimFreqs = NaN; %tactile stimulus frequencies (Hz)
    stimAmps = NaN; %tactile stimulus amplitude (voltage)
    stimDurs = NaN; %tactile stimulus durations (ms)
    audioTrials = NaN; %logical list of audio (stim) only trials
    imgDepth = NaN; % depth of focal plane (micrometers)
    graspTimes = NaN; %grasp times [grasp ON, grasp OFF; ...]
    graspDurs = NaN; %difference of pairwise graspON/OFF
    lickTimes = NaN; %array of clock times when a lick was detected
    vasculature1x = NaN; %photo of blood vessels at brain surface, 1x digital zoom
    vasculature2x = NaN; %photo of blood vessels at brain surface, 2x digital zoom (optional)
    vasculature3x = NaN; %photo of blood vessels at brain surface, 3x digital zoom
    tformRFPbasis = NaN; %tform object from
    zStackMapIndx = NaN; %new
    globalDepth = NaN; %new
    
    % GFP SOMA ROI DETAILS
    somaticROI_PixelLists = NaN; %from CDiesters extractor.m
    somaticROIBoundaries = NaN; %from CDiesters extractor.m
    somaticROICenters = NaN; %from CDiesters extractor.m
    somaticRoiCounter = NaN; %from CDiesters extractor.m
    somaticROIs = NaN; %from CDiesters extractor.m

    % RFP SOMA ROI DETAILS
    redSomaticROI_PixelLists = NaN; %from CDiesters extractor.m
    redSomaticROIBoundaries = NaN; %from CDiesters extractor.m
    redSomaticROICenters = NaN; %from CDiesters extractor.m
    redSomaticRoiCounter = NaN; %from CDiesters extractor.m
    redSomaticROIs = NaN; %from CDiesters extractor.m

    % NEUROPIL ROI DETAILS (PAIRED w/ GFP SOMAs)
    neuropilROI_PixelLists = NaN; %from CDiesters extractor.m
    neuropilROIBoundaries = NaN; %from CDiesters extractor.m
    neuropilROICenters = NaN; %from CDiesters extractor.m
    neuropilRoiCounter = NaN; %from CDiesters extractor.m
    neuropilROIs = NaN %from CDiesters extractor.m
    
      
    end
    %======================================================================
    methods %Ordinary Methods
%         function obj = ClassName(arg1)
%             obj.PropertyName = arg1;
%             ...
%         end
%     
%         function ordinaryMethod(obj,arg1)
%             ...
%         end
    end 
    %======================================================================
    methods (Static) %Static Methods
        
        function f = rawTrace(fSignal,frameTimes)
            if any(size(fSignal) == 1) && length(size(fSignal)) == 2
                %fluorescence signal is [1 x n] array
                f = figure;
                plot(fSignal); 
                if nargin > 1
                    hold on;
                    for n = 1 : length(frameTimes)
                        plot([frameTimes(n) frameTimes(n)],[min(fSignal) max(fSignal)],'k-');
                    end
                end
            else
                %fluorescence signal is an image stack [N1 x N1 x N2]
                if size(fSignal,1) == size(fSignal,2) && length(size(fSignal)) == 3
                    fTrace = squeeze(mean(fSignal,[1 2]));
                    f = figure;
                    plot(fTrace); 
                    if nargin > 1
                        hold on;
                        for n = 1 : length(frameTimes)
                            frame_time = floor(frameTimes(n) / 3);
                            plot([frame_time frame_time],[min(fTrace) max(fTrace)],'k-');
                        end
                    end
                else
                    error("fSignal provided to rawTrace() is unrecognized");
                end
            end
        end
        %------------------------------------------------------------------
        function f = psth(fSignal,stimFrameTimes)
            numReps = size(fSignal,1);
            xRange = [-0.9:0.1:6];
            frame_times = floor(stimFrameTimes/3);
            f = figure; hold on;
            for nr = 1 : numReps
                allResps = zeros(length(frame_times),length(xRange));
                for n=2:length(frame_times) %skip first so we can normalize
                   if frame_times(n) < 10 || frame_times(n) == 0
                       break
                   end
                   cur_bucket = fSignal(nr,frame_times(n)-9:frame_times(n)+60);
                   prior = mean(cur_bucket(1:9)); %0.9 sec prior to stim ON;
                   allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                   plot(xRange,allResps(n,:),'Color',[0.9 0.9 0.9]);
                end
                plot(xRange,mean(allResps,1),'k');
                title(strcat("PSTH of ROI#",num2str(nr)));
                hold off;
            end
        end
        %------------------------------------------------------------------
        function [tuningMat, tuningMatSEM, yRange] = psth2(fSignal,stimFrameTimes,stimFreqs,audioTrials,nROI,threshSTD)
            % PSTH FUNCTION FOR TACTILE STIM ONly Audio-stim trials are
            % handled by psth_audio()
            %generate key map for freq / trials

            preStimT = 1; %sec before sitm on;
            postStimT = 5; %sec after stim on;
            frameDur = 0.1;
            Hz = 10; %frameRate
            allFreqs = unique(stimFreqs);
            numTrials = length(stimFreqs);
            keyMap = zeros(numTrials,length(allFreqs));
            for n=1:length(allFreqs)
                keyMap(:,n) = stimFreqs == allFreqs(n);
            end
            numROIs = size(fSignal,1);
            xRange = [-1*(preStimT-frameDur):frameDur:postStimT];
            tuningMat = zeros(1,length(allFreqs));
            spX = ceil(length(allFreqs)/2);
            spY = 2;
            tuningMatSEM = zeros(1,length(allFreqs));
            tuningMatSTD = zeros(1,length(allFreqs));
            yRange = [0,0];
            for ns = 1 : length(allFreqs)
                subplot(spX,spY,ns); hold on;
                frame_times = floor(stimFrameTimes/3);
                frame_times(frame_times == 0) = [];
                thisKey = keyMap(:,ns);
                trialList = find(thisKey);
                frame_times = frame_times(logical(thisKey(1:length(frame_times))));
                allResps = zeros(length(frame_times),length(xRange));
                for n=2:length(frame_times) %skip first so we can normalize
                   numPreFrames = Hz*(preStimT-frameDur);
                   numPostFrames = postStimT * Hz;
                   if frame_times(n) < preStimT*Hz || frame_times(n) == 0
                       break
                   end
                   if frame_times(n)+numPostFrames > size(fSignal,2)
                       continue;
                   end
                   if audioTrials(trialList(n)) == 1
                       % SKIP ALL AUDIO TRIALS FOR THIS PSTH
                       continue;
                   end
                   cur_bucket = fSignal(frame_times(n)-numPreFrames:frame_times(n)+numPostFrames);
                   prior = mean(cur_bucket(1:numPreFrames)); %0.9 sec prior to stim ON;
                   allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                end
                audioTrialsFlag = (sum(allResps,2) == 0);
                allResps(audioTrialsFlag,:) = [];
                meanResp = mean(allResps,1);
                stdResp = std(allResps,1);
                semResp = std(allResps,1)/sqrt(size(allResps,1));
                aboveThresh = allResps>(meanResp+threshSTD.*stdResp);
                belowThresh = allResps<(meanResp-threshSTD.*stdResp);
                outlierTrials = sum(aboveThresh+belowThresh,2);
                outlierTrials = outlierTrials > 0;
                [tuningMat(ns), indxMax] = max(meanResp(10:30)); %highest dF/F in 2s response window post stim
                tuningMatSTD(ns) = std(allResps(:,indxMax));
                tuningMatSEM(ns) =std(allResps(:,indxMax))/sqrt(size(allResps,1)); %calculates SEM for max dF/F resp
                shadedErrorBar(xRange,meanResp,semResp); hold on;
                title(strcat("PSTH of ROI#",num2str(nROI),"@ ",num2str(allFreqs(ns)),"Hz"));
                newYRange = [min(meanResp)*1.2,max(meanResp)*1.2];
                if newYRange(1) < yRange(1)
                    yRange(1) = newYRange(1);
                end
                if newYRange(2) > yRange(2)
                    yRange(2) = newYRange(2);
                end 
                xlabel("Time (Seconds)");
                ylabel("dF/F")
                hold off;
            end
            for n = 1 : ns
                subplot(spX,spY,n);
                ax = gca;
                ax.YLim = yRange;
            end
        end
        %------------------------------------------------------------------
        function f = simpleTuningCurve(stimFreqs,tuningMat,tuningMatSTD,nROI)
            allFreqs = unique(stimFreqs);
            bar(allFreqs,tuningMat);
            hold on;
            errorbar(allFreqs,tuningMat,tuningMatSTD);
            title(strcat("Tuning Curve ROI ",num2str(nROI)));
            xlabel('Frequency (Hz)');
            xtickangle(45)
            ylabel('max dF/F');
            hold off;
        end
        %------------------------------------------------------------------
        function f = tuningCurve(fSignal,stimFrameTimes,stimFreqs,nROI)
            %generate key map for freq / trials
            allFreqs = unique(stimFreqs);
            numTrials = length(stimFreqs);
            keyMap = zeros(numTrials,length(allFreqs));
            for n=1:length(allFreqs)
                keyMap(:,n) = stimFreqs == allFreqs(n);
            end
            numROIs = size(fSignal,1);
            xRange = [-0.9:0.1:4];
            tuningMat = zeros(numROIs,length(allFreqs));
            for nr = 1 : numROIs
                spX = ceil(length(allFreqs)/2);
                spY = 2;
                tuningMatStd = zeros(1,length(allFreqs));
                for ns = 1 : length(allFreqs)
                    frame_times = floor(stimFrameTimes/3);
                    frame_times(frame_times == 0) = [];
                    thisKey = keyMap(:,ns);
                    frame_times = frame_times(logical(thisKey(1:length(frame_times))));
                    allResps = zeros(length(frame_times),length(xRange));
                    for n=2:length(frame_times) %skip first so we can normalize
                       if frame_times(n) < 10 || frame_times(n) == 0
                           break
                       end
                       cur_bucket = fSignal(nr,frame_times(n)-9:frame_times(n)+40);
                       prior = mean(cur_bucket(1:9)); %0.9 sec prior to stim ON;
                       allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                    end
                    response = mean(allResps,1);
                    [tuningMat(nr,ns), indxMax] = max(response(10:30)); %highest dF/F in 2s response window post stim
                    tuningMatStd(ns) = std(allResps(:,indxMax))/sqrt(size(allResps,1)); %calculates SEM for max dF/F resp
                end
            end
            for n=1:numROIs
                bar(unique(stimFreqs),tuningMat(n,:)); hold on;
                errorbar(unique(stimFreqs),tuningMat(n,:),tuningMatStd);
                title(strcat("Tuning Curve ROI ",num2str(nROI)));
                xlabel('Frequency (Hz)');
                xtickangle(45)
                ylabel('max dF/F');
                hold off;
            end
        end
        %------------------------------------------------------------------
        function [f1, f2] = psthAudioTrials(fSignal,stimFrameTimes,stimFreqs,audioTrials,nROI,yRange)
            %generate key map for freq / trials
            allFreqs = unique(stimFreqs);
            numTrials = length(stimFreqs);
            keyMap = zeros(numTrials,length(allFreqs));
            for n=1:length(allFreqs)
                keyMap(:,n) = stimFreqs == allFreqs(n);
            end
            numROIs = size(fSignal,1);
            xRange = [-0.9:0.1:7];
            tuningMat = zeros(numROIs,length(allFreqs));
            overall = zeros(length(allFreqs),size(xRange,2));
            for nr = 1 : numROIs
                spX = ceil(length(allFreqs)/2);
                spY = 2;
                tuningMatStd = zeros(1,length(allFreqs));
                for ns = 1 : length(allFreqs) 
                    frame_times = floor(stimFrameTimes/3);
                    frame_times(frame_times == 0) = [];
                    thisKey = keyMap(:,ns);
                    trialList = find(thisKey);
                    frame_times = frame_times(logical(thisKey(1:length(frame_times))));
                    allResps = zeros(length(frame_times),length(xRange));
                    for n=2:length(frame_times) %skip first so we can normalize
                       if frame_times(n) < 10 || frame_times(n) == 0
                           break
                       end
                       if frame_times(n)+70 > size(fSignal,2)
                           continue;
                       end
                       if audioTrials(trialList(n)) == 0
                           continue;
                       end
                       cur_bucket = fSignal(nr,frame_times(n)-9:frame_times(n)+70);
                       prior = mean(cur_bucket(1:9)); %0.9 sec prior to stim ON;
                       allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                       %plot(xRange,allResps(n,:),'Color',[0.7 0.7 0.7]); hold on;
                    end
                    overall(ns,:) = mean(allResps,1);
%                     plot(xRange,response,'k'); hold on;
%                     [tuningMat(nr,ns), indxMax] = max(response(10:30)); %highest dF/F in 2s response window post stim
%                     tuningMatStd(ns) = std(allResps(:,indxMax))/sqrt(size(allResps,1)); %calculates SEM for max dF/F resp
%                     title(strcat("PSTH of ROI#",num2str(nROI),"@ ",num2str(allFreqs(ns)),"Hz"));
%                     hold off;
                end
                %plot(xRange,mean(overall,1),'Color',[0,0,0]);
                signalMean = mean(overall,1);
                signalSTD = std(overall,1);
                signalSEM = signalSTD/sqrt(size(overall,1));
                shadedErrorBar(xRange,signalMean,signalSEM); hold on;
                title(strcat("Audio-Stim PSTH of ROI#",num2str(nROI)));
                xlabel('Time post-stim (sec)');
                ylabel('dF/F');
                ax = gca;
                ax.YLim = yRange;
                hold off;
            end
%             hold off;
%             f2 = figure; %tuning Curves
%             for n=1:numROIs
%                 subplot(1,numROIs,n);
%                 bar(unique(stimFreqs),tuningMat(n,:)); hold on;
%                 errorbar(unique(stimFreqs),tuningMat(n,:),tuningMatStd);
%                 title(strcat("Tuning Curve ROI ",num2str(n)));
%                 xlabel('Frequency (Hz)');
%                 xtickangle(45)
%                 ylabel('max dF/F');
%             end
        end
        %------------------------------------------------------------------
        function f = psthAll(fSignal,stimFrameTimes)
            numReps = size(fSignal,1);
            xRange = [-0.9:0.1:6];
            frame_times = floor(stimFrameTimes/3);
            f = figure; hold on;
            for nr = 1 : numReps
                allResps = zeros(length(frame_times),length(xRange));
                for n=2:length(frame_times) %skip first so we can normalize
                   if frame_times(n) < 10 || frame_times(n) == 0
                       break
                   elseif frame_times(n)+60 > size(fSignal,2)
                       break
                   end
                   cur_bucket = fSignal(nr,frame_times(n)-9:frame_times(n)+60);
                   prior = mean(cur_bucket(1:9)); %0.9 sec prior to stim ON;
                   allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                   %plot(xRange,allResps(n,:),'Color',[0.9 0.9 0.9]);
                end
                shadedErrorBar(xRange,mean(allResps,1),std(allResps,1)/sqrt(length(frame_times)))
                plot(xRange,mean(allResps,1),'k');
                title(strcat("PSTH of ROI#",num2str(nr)));
                hold off;
            end
        end
        %------------------------------------------------------------------
        function f = movie(images,fps,stimOnFrameNum)
            frameNums = floor(stimOnFrameNum/3);
            figure;
            sticky = 0;
            for n = 1 : size(images,3)
                sticky = sticky - 1;
                imagesc(images(:,:,n));
                if sum(frameNums == n) == 1
                    %this frame was a stim presentation
                    title(strcat("STIM Frame: ",num2str(n)));
                    sticky = 10;
                elseif sticky > 1
                    title(strcat("STIM Frame: ",num2str(n)));
                else
                    title(strcat("Frame: ",num2str(n)));
                end
                time = GetSecs;
                while GetSecs - time < (1/fps)
                     pause(0)
                end
            end 
        end
        %------------------------------------------------------------------
        function f = motionCorrection(imageStack)
%             figure;
            f = zeros(size(imageStack,1),size(imageStack,2),size(imageStack,3));
            for n = 2 : size(imageStack,3)
                if rem(n,10) == 0
                    fprintf(".")
                end
                if rem(n,100) == 0
                    fprintf(strcat(num2str(n/10),"% /n"));
                end
                fixed = imageStack(:,:,n-1);
                moving = imageStack(:,:,n);
                subplot(1,2,1);
                imshowpair(fixed,moving,'Scaling','joint');
                title("Pre-Correction");
                %imregister set-up for motion corrections
                [optimizer, metric] = imregconfig('monomodal');
%                 optimizer.InitialRadius = 0.009;
%                 optimizer.Epsilon = 1.5e-4;
%                 optimizer.GrowthFactor = 1.01;
%                 optimizer.MaximumIterations = 300;
                f(:,:,n) = imregister(moving, fixed, 'affine', optimizer, metric);
                %Output Plots
%                 subplot(1,2,2);
%                 imshowpair(fixed, f(:,:,n),'Scaling','joint');
%                 title("Corrected");
            end
        end
        %------------------------------------------------------------------
        function f = imStack(images,nStart,nEnd)
            f = figure;
            imAvg = mean(images(:,:,nStart:nEnd),3);
            imagesc(imAvg); hold on;
            title(strcat("Images ",num2str(nStart)," : ",num2str(nEnd)));
        end
        %------------------------------------------------------------------
        function f = imDual(imGFP, imRFP)
           if length(size(imGFP)) == 3
               %3 average through stack
               imGFP = mean(imGFP,3);
           end
           if size(imRFP,1) > size(imGFP,1)
               %imGFP was likely resized to 0.5 scale
               imRFP = imresize(imRFP,0.5);
           end
           f = figure; hold on;
           pxNum = size(imRFP,1);
           rfpAvg = mean(mat2gray(imRFP),'all');
           gfpAvg = mean(mat2gray(imGFP),'all');
           if rfpAvg < gfpAvg
              redBoost = (gfpAvg/rfpAvg)/2;
              imRGB = cat(3,mat2gray(imRFP).*redBoost,mat2gray(imGFP),zeros(pxNum, pxNum));
           else
              imRGB = cat(3,mat2gray(imRFP),mat2gray(imGFP),zeros(pxNum, pxNum));
           end
           imshow(imRGB);
           hold off;
        end
        %------------------------------------------------------------------
        function [outStruct] = getResponsivity(fSignal,stimFrameTimes,stimFreqs,audioTrials,respWindow)
            % getResponsivity
            % getResponsivity(fSignal,stimFrameTimes,stimFreqs,audioTrials,respWindow)
            %  - Seperates trials by frequency.
            %  - calculates ranksum of freq. avg. brightness 1(sec) pre-stim & n(sec) post-stim
            %  - Ranksum's z-stat is then used in freq. tuning curve
            %  - Also creates freq. tuning curve from dF/F
            
            if all(isnan(fSignal))
                disp("Found trial of all NaNs, skipping...");
                return
            end
            
            % generate key map for seperating freq / trials
            preStimT = 1; %sec before sitm on;
            postStimT = 5; %sec after stim on;
            frameDur = 0.1;
            Hz = 10; %frameRate
            allFreqs = unique(stimFreqs);
            numTrials = length(stimFreqs);
            keyMap = zeros(numTrials,length(allFreqs));
            for n=1:length(allFreqs)
                keyMap(:,n) = stimFreqs == allFreqs(n);
            end
            
            % Preallocate variables
            outStruct = struct();
            for i = 1:length(allFreqs)
                outStruct(1).(strcat('f',num2str(allFreqs(i)))) = struct();
            end
            respFreqZP = zeros(size(allFreqs,1),2);
            numROIs = size(fSignal,1);
            xRange = [-1*(preStimT-frameDur):frameDur:postStimT];
            tuningMat = zeros(1,length(allFreqs));
            spX = ceil(length(allFreqs)/2);
            spY = 2;
            tuningMatSEM = zeros(1,length(allFreqs));
            tuningMatSTD = zeros(1,length(allFreqs));
            yRange = [0,0];
            respCalcWin = [((preStimT+frameDur)*Hz) ((respWindow+preStimT)*Hz)];
            
            
            % Main loop through all frequencies
            for ns = 1 : length(allFreqs)
                %subplot(spX,spY,ns); hold on;
                structName = strcat('f',num2str(allFreqs(ns)));
                frame_times = floor(stimFrameTimes/3);
                frame_times(frame_times == 0) = [];
                thisKey = keyMap(:,ns);
                trialList = find(thisKey);
                frame_times = frame_times(logical(thisKey(1:length(frame_times))));
                allResps = zeros(length(frame_times),length(xRange));
                for n=2:length(frame_times) %skip first so we can normalize
                   numPreFrames = Hz*(preStimT-frameDur);
                   numPostFrames = postStimT * Hz;
                   if frame_times(n) < preStimT*Hz || frame_times(n) == 0
                       break
                   end
                   if frame_times(n)+numPostFrames > size(fSignal,2)
                       continue;
                   end
                   if audioTrials(trialList(n)) == 1
                       % SKIP ALL AUDIO TRIALS FOR THIS PSTH
                       continue;
                   end
                   cur_bucket = fSignal(frame_times(n)-numPreFrames:frame_times(n)+numPostFrames);
                   prior = mean(cur_bucket(1:numPreFrames)); %0.9 sec prior to stim ON;
                   allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                end
                audioTrialsFlag = (sum(allResps,2) == 0);
                allResps(audioTrialsFlag,:) = [];
                meanResp = mean(allResps,1);
                stdResp = std(allResps,1);
                semResp = std(allResps,1)/sqrt(size(allResps,1));
                outStruct.(structName).allTraces = allResps;
                outStruct.(structName).stdTraces = stdResp;
                outStruct.(structName).semTraces = semResp;
                outStruct.(structName).avgTrace = meanResp;
                outStruct.(structName).preStim = mean(allResps(:,1:10),2);
                outStruct.(structName).avgRespPostStim = mean(allResps(:,respCalcWin(1):respCalcWin(2)),2); % avg resp window for each trial
                outStruct.(structName).stdRespPostStim = std(mean(allResps(:,respCalcWin(1):respCalcWin(2)),2)); %std of avgResp2secWin
                outStruct.(structName).semRespPostStim = std(mean(allResps(:,respCalcWin(1):respCalcWin(2)),2))/sqrt(size(allResps,1));
                try 
                    [P,H,STATS] = ranksum(outStruct.(structName).avgRespPostStim,outStruct.(structName).preStim);
                catch
                    disp(strcat("Error for Freq. ",num2str(ns)));
                end
                    
                if max(contains(fieldnames(STATS),'zval')) == 1
                    outStruct.(structName).responsivityZ = STATS.zval;
                else
                    outStruct.(structName).responsivityZ = NaN;
                end
                outStruct.(structName).responsivityP = P;
            end
            outStruct.xRange = xRange;
            out = outStruct;
            
            
            % START MannyanalysisNEW
            fieldNames = fieldnames(out);
            fieldNames(end) = [];
            out.pVals = zeros(1,length(fieldNames));
            out.zVals = zeros(1,length(fieldNames));
            matAvgFreqTrials = NaN(6,50);
            for nn = 1:length(fieldNames)
                respFreqZP(nn,1) = out.(fieldNames{nn}).responsivityZ;
                respFreqZP(nn,2) = out.(fieldNames{nn}).responsivityP;
                matAvgFreqTrials(nn,1:size(out.(fieldNames{nn}).avgRespPostStim,1)) = out.(fieldNames{nn}).avgRespPostStim;
            end
            selectivity = (max(abs(respFreqZP(:,1)))-min(abs(respFreqZP(:,1))))/...
                (max(abs(respFreqZP(:,1)))+min(abs(respFreqZP(:,1))));
            
            
%             % Plotting
%             f = figure;
%             f.Position = [100 100 1000 400];
            
            %Tuning curve from Z Statistic
%             subplot(1,4,1);
%             bar([100,300,500,700,900,1100],respFreqZP(:,1)); hold on;
%             plot([0, 1200],[1.96,1.96],'g--');
%             plot([0, 1200],[-1.96,-1.96],'g--');
%             ax = gca;
%             ax.YLim(2) = ax.YLim(2)*1.1;
%             title(strcat("Tuning curve Z-Stat."));
%             subtitle(strcat("\Delta","F post-stim per freq."),'FontSize',6);
%             xlabel("Frequency (Hz)");
%             ylabel("Z Statistic");
%             hold off;

            %Tuning curve from Avg. dF/F & Median dF/F post-stim
            avgPrePostSem = zeros(size(respFreqZP,1),2); %avg responses / freq.
            medPrePostMad = zeros(size(respFreqZP,1),2); %median responses / freq.
            allResps = zeros(1000,2);
            counter = 0;
            for ne = 1 : size(respFreqZP,1)
                avgPrePostSem(ne,1) = mean(out.(fieldNames{ne}).avgRespPostStim) - mean(out.(fieldNames{ne}).preStim);
                avgPrePostSem(ne,2) = out.(fieldNames{ne}).semRespPostStim;
                medPrePostMad(ne,1) = median(out.(fieldNames{ne}).avgRespPostStim) - median(out.(fieldNames{ne}).preStim);
                medPrePostMad(ne,2) = mad(out.(fieldNames{ne}).avgRespPostStim,1);
                finalIndx = size(out.(fieldNames{ne}).avgRespPostStim,1);
                allResps(counter+1:counter+finalIndx,1) = out.(fieldNames{ne}).avgRespPostStim;
                allResps(counter+1:counter+finalIndx,2) = out.(fieldNames{ne}).preStim;
                counter = counter + finalIndx;
            end
            
            subplot(1,2,1);
            bar([100,300,500,700,900,1100],avgPrePostSem(:,1),'FaceColor',[0.15,0.85,0.72]); hold on;
            [tuningMeanSEM] = WithinSubj_StdError(matAvgFreqTrials); %SEM calc across Subj. (Manny)
            er = errorbar([100,300,500,700,900,1100],avgPrePostSem(:,1),tuningMeanSEM);
            er.Color = [0 0 0];
            er.LineStyle = 'none';
            title("Tuning curve (mean resp.)");
            subtitle(strcat("\Delta","F/F Post-stim w/ SEM"),'FontSize',6);
            xlabel("Frequency (Hz)");
            ylabel(strcat("\Delta","F/F"));
            hold off;
            
%             subplot(1,4,3);
%             bar([100,300,500,700,900,1100],medPrePostMad(:,1)); hold on;
%             er = errorbar([100,300,500,700,900,1100],medPrePostMad(:,1),medPrePostMad(:,2));
%             er.Color = [0 0 0];
%             er.LineStyle = 'none';
%             title("Tuning curve (med. resp.)");
%             subtitle(strcat("\Delta","F/F Post-stim w/ MAD"),'FontSize',6);
%             xlabel("Frequency (Hz)");
%             ylabel(strcat("\Delta","F/F"));
%             hold off;
                
            [allP,allH,allStats] = ranksum(allResps(:,1),allResps(:,2));
            subplot(1,2,2);
            if allP < 0.05
                boxplot(cat(2,allResps(:,2),allResps(:,1)),'Labels',{'F pre-stim','F post-stim'},'PlotStyle','compact','Colors',[0.15,0.85,0.72]);hold on;
            else
                boxplot(cat(2,allResps(:,2),allResps(:,1)),'Labels',{'F pre-stim','F post-stim'});hold on;
            end
            title("Responsivity (all frequencies)");
            subtitle(strcat("P=",num2str(allP)," Zstat=",num2str(allStats.zval)),'FontSize',6);
            ylabel(strcat("\Delta","F/F"));
            
            %Find significance threshold >1.96 <-1.96 zVal
            ax = gca;
            if ax.YLim(1)<1.96 && ax.YLim(2)>1.96
                plot([ax.XLim(1) ax.XLim(2)],[1.96 1.96],'g--');
            end
            if ax.YLim(1)<-1.96 && ax.YLim(2)>-1.96
                plot([ax.XLim(1) ax.XLim(2)],[-1.96 -1.96],'g--');
            end

            %find the max and min response trace values for commmon axes
            axMinMax = [0,0];
            for nn = 1 : 6
                if max(out.(fieldNames{nn}).avgTrace) > axMinMax(2)
                    axMinMax(2) = max(out.(fieldNames{nn}).avgTrace);
                end
                if min(out.(fieldNames{nn}).avgTrace) < axMinMax(1)
                    axMinMax(1) = min(out.(fieldNames{nn}).avgTrace);
                end
            end
            axMinMax(1) = axMinMax(1)-(diff(axMinMax)/10); %pad 10%
            axMinMax(2) = axMinMax(2)+(diff(axMinMax)/10); %pad 10%
            
            %add additional outputs
            outStruct.respFreqZP = respFreqZP;
            outStruct.selectivity_bestMinWorst = selectivity;
            outStruct.responsivityZ = allStats.zval;
            outStruct.responsivityP = allP;
            outStruct.yAxMinMax = axMinMax;
            outStruct.tuningMeanSem = tuningMeanSEM;
            outStruct.avgRespPost = avgPrePostSem(:,1);
            outStruct.medRespPost = medPrePostMad(:,1);
            
            %END MannyanalysisNEW
        end
        %------------------------------------------------------------------
        function [outStruct] = getTraces(fSignal,stimFrameTimes,stimFreqs,audioTrials,respWindow)
            % getTraces.m
            % getTraces(fSignal,stimFrameTimes,stimFreqs,audioTrials,respWindow)
            %  - Seperates trials by frequency.
            %  - calculates ranksum of freq. avg. brightness 1(sec) pre-stim & n(sec) post-stim
            %  - Ranksum's z-stat is then used in freq. tuning curve
            %  - Also creates freq. tuning curve from dF/F
            
            
            % generate key map for seperating freq / trials
            preStimT = 1; %sec before sitm on;
            postStimT = 5; %sec after stim on;
            frameDur = 0.1;
            Hz = 10; %frameRate
            allFreqs = unique(stimFreqs);
            numTrials = size(stimFrameTimes,1); %calculate numTrials from number of stimOnTimes

            %check for discrepancy between nTrials and other variables
            if numTrials < size(stimFreqs,1)
                stimFreqs = stimFreqs(1:numTrials); %trim off trials with no stimOnTime
            elseif numTrials < size(audioTrials,1)
                audioTrials = audioTrials(1:numTrials); %trim off trials with no stimOnTime
            end

            %create logical map of trial numbers for each frequency
            keyMap = zeros(numTrials,length(allFreqs));
            for n=1:length(allFreqs)
                keyMap(:,n) = stimFreqs == allFreqs(n);
            end
            
            % Preallocate variables
            outStruct = struct();
            for i = 1:length(allFreqs)
                outStruct(1).(strcat('f',num2str(allFreqs(i)))) = struct();
            end
            respFreqZP = zeros(size(allFreqs,1),2);
            numROIs = size(fSignal,1);
            xRange = [-1*(preStimT-frameDur):frameDur:postStimT];
            tuningMat = zeros(1,length(allFreqs));
            spX = ceil(length(allFreqs)/2);
            spY = 2;
            tuningMatSEM = zeros(1,length(allFreqs));
            tuningMatSTD = zeros(1,length(allFreqs));
            yRange = [0,0];
            respCalcWin = [((preStimT+frameDur)*Hz) ((respWindow+preStimT)*Hz)];
            
            
            % Main loop through all frequencies
            allData = {};      
            repNumber = 35; %max number of repetitions per stimulus;
            for ns = 1 : length(allFreqs)
                %subplot(spX,spY,ns); hold on;
                myDataAvg = [];
                structName = strcat('f',num2str(allFreqs(ns)));
                frame_times = floor(stimFrameTimes/3);
                frame_times(frame_times == 0) = [];
                thisKey = keyMap(:,ns);
                trialList = find(thisKey);
                frame_times = frame_times(logical(thisKey(1:length(frame_times))));
                allResps = zeros(length(frame_times),length(xRange));
                for n=2:length(frame_times) %skip first so we can normalize
                   numPreFrames = Hz*(preStimT-frameDur);
                   numPostFrames = postStimT * Hz;
                   if frame_times(n) < preStimT*Hz || frame_times(n) == 0
                       break
                   end
                   if frame_times(n)+numPostFrames > size(fSignal,2)
                       continue;
                   end
                   if audioTrials(trialList(n)) == 1
                       % SKIP ALL AUDIO TRIALS FOR THIS PSTH
                       continue;
                   end
                   cur_bucket = fSignal(frame_times(n)-numPreFrames:frame_times(n)+numPostFrames);
                   prior = mean(cur_bucket(1:numPreFrames)); %0.9 sec prior to stim ON;
                   allResps(n,:) = (cur_bucket - prior)./prior; %normalization function
                end
                audioTrialsFlag = (sum(allResps,2) == 0);
                allResps(audioTrialsFlag,:) = [];
                meanResp = mean(allResps,1);
                stdResp = std(allResps,1);
                semResp = std(allResps,1)/sqrt(size(allResps,1));
                outStruct.(structName).allTraces = allResps;
                outStruct.(structName).stdTraces = stdResp;
                outStruct.(structName).semTraces = semResp;
                outStruct.(structName).avgTrace = meanResp;
                outStruct.(structName).preStim = mean(allResps(:,1:10),2);
                outStruct.(structName).avgRespPostStim = mean(allResps(:,respCalcWin(1):respCalcWin(2)),2); % avg resp window for each trial
                outStruct.(structName).stdRespPostStim = std(mean(allResps(:,respCalcWin(1):respCalcWin(2)),2)); %std of avgResp2secWin
                outStruct.(structName).semRespPostStim = std(mean(allResps(:,respCalcWin(1):respCalcWin(2)),2))/sqrt(size(allResps,1));

                myData = outStruct.(structName).avgRespPostStim;
%                 [P,H,STATS] = ranksum(outStruct.(structName).avgRespPostStim,outStruct.(structName).preStim);
                [P,H,STATS] = signrank(outStruct.(structName).avgRespPostStim);
                if max(contains(fieldnames(STATS),'zval')) == 1
                    outStruct.(structName).responsivityZ = STATS.zval;
                else
                    outStruct.(structName).responsivityZ = NaN;
                end
                outStruct.(structName).responsivityP = P;
                if length(myData) > repNumber
                    repNumber = length(myData);
                end
                if length(myData) < repNumber
                    myData = cat(1,myData, nan(repNumber-length(myData),1));
                else
                end
                allData = cat(2,allData,myData);
            end
            outStruct.xRange = xRange;

            outStruct.isSelectiveP = kruskalwallis(cell2mat(allData),[],'off');
            if outStruct.isSelectiveP < 0.05
                [~,indx] = max(cellfun(@mean,allData));
                outStruct.prefFreq = allFreqs(indx);
            else
                outStruct.prefFreq = NaN;
            end
            out = outStruct;

            
            
            % START MannyanalysisNEW
            fieldNames = fieldnames(out);
            fieldNames(end-1:end) = [];
            out.pVals = zeros(1,length(fieldNames));
            out.zVals = zeros(1,length(fieldNames));
            matAvgFreqTrials = NaN(6,50);
            for nn = 1: length(allFreqs)
                respFreqZP(nn,1) = out.(fieldNames{nn}).responsivityZ;
                respFreqZP(nn,2) = out.(fieldNames{nn}).responsivityP;
                matAvgFreqTrials(nn,1:size(out.(fieldNames{nn}).avgRespPostStim,1)) = out.(fieldNames{nn}).avgRespPostStim;
            end
            selectivity = (max(abs(respFreqZP(:,1)))-min(abs(respFreqZP(:,1))))/...
                (max(abs(respFreqZP(:,1)))+min(abs(respFreqZP(:,1))));
            

            %Tuning curve from Avg. dF/F & Median dF/F post-stim
            avgPrePostSem = zeros(size(respFreqZP,1),2); %avg responses / freq.
            medPrePostMad = zeros(size(respFreqZP,1),2); %median responses / freq.
            allResps = zeros(1000,2);
            counter = 0;
            for ne = 1 : size(respFreqZP,1)
                avgPrePostSem(ne,1) = mean(out.(fieldNames{ne}).avgRespPostStim) - mean(out.(fieldNames{ne}).preStim);
                avgPrePostSem(ne,2) = out.(fieldNames{ne}).semRespPostStim;
                medPrePostMad(ne,1) = median(out.(fieldNames{ne}).avgRespPostStim) - median(out.(fieldNames{ne}).preStim);
                medPrePostMad(ne,2) = mad(out.(fieldNames{ne}).avgRespPostStim,1);
                finalIndx = size(out.(fieldNames{ne}).avgRespPostStim,1);
                allResps(counter+1:counter+finalIndx,1) = out.(fieldNames{ne}).avgRespPostStim;
                allResps(counter+1:counter+finalIndx,2) = out.(fieldNames{ne}).preStim;
                counter = counter + finalIndx;
            end
            

            %find the max and min response trace values for commmon axes
            axMinMax = [0,0];
            for nn = 1 : 6
                if max(out.(fieldNames{nn}).avgTrace) > axMinMax(2)
                    axMinMax(2) = max(out.(fieldNames{nn}).avgTrace);
                end
                if min(out.(fieldNames{nn}).avgTrace) < axMinMax(1)
                    axMinMax(1) = min(out.(fieldNames{nn}).avgTrace);
                end
            end
            axMinMax(1) = axMinMax(1)-(diff(axMinMax)/10); %pad 10%
            axMinMax(2) = axMinMax(2)+(diff(axMinMax)/10); %pad 10%

            %add subbplots for traces per frequency.
%             for nn = 1 : 6
%                 subplot(2,3,nn); hold on;
%                 if respFreqZP(nn,2) < 0.05
%                     shadedErrorBar(out.xRange,out.(fieldNames{nn}).avgTrace,out.(fieldNames{nn}).semTraces,...
%                         'lineProps',{'g-','markerfacecolor',[0.9290 0.6940 0.1250]});
%                 else 
%                     shadedErrorBar(out.xRange,out.(fieldNames{nn}).avgTrace,out.(fieldNames{nn}).semTraces,...
%                         'lineProps',{'r-','markerfacecolor',[0.9290 0.6940 0.1250]});
%                 end
%                 titleText = strcat("Response at ",fieldNames{nn}," Hz");
%                 title(titleText);
%                 xlabel("Time Post-Stim (sec.)");
%                 ylabel("dF/F");
%                 ax = gca;
%                 ax.YLim = axMinMax;
%                 yLoc = ax.YLim(2)-((ax.YLim(2)-ax.YLim(1))/10);
%                 plot([0 0],[ax.YLim(1) ax.YLim(2)],'b--');
%                 plot([respWindow respWindow],[ax.YLim(1) ax.YLim(2)],'b--');
%                 figText = strcat("Z-Val = ",num2str(respFreqZP(nn,1)));
%                 figText = [figText;strcat("P-Val = ",num2str(respFreqZP(nn,2)))];
%                 text(-0.8,yLoc,figText);
%                 hold off;
%             end
            
            %add additional outputs
            outStruct.respFreqZP = respFreqZP;
            outStruct.selectivity_bestMinWorst = selectivity;
            outStruct.yAxMinMax = axMinMax;
            %outStruct.tuningMeanSem = tuningMeanSEM;
            outStruct.avgRespPost = avgPrePostSem(:,1);
            outStruct.medRespPost = medPrePostMad(:,1);
            
            %END MannyanalysisNEW
        end
        %------------------------------------------------------------------
        function r = multiplyBy(obj,n)
            %example internal function 
            r = [obj.Value]*n;
        end
    end
end