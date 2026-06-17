% FAS_simultaneousDecode.m
%
% one more script to rule them all and in the darkness bind them
% decode all 25 stimuli simultaneously, with each as a unique index

%% perform population decoding

nIter = 100; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
targetAmp = 5; % 0.33,0.76,1.6,3.5,7.7
usePCA = true; % use PCA step prior to model training
kStart = 1;
kComp = 25; %number of PCs to use in model training / testing
useFreqEffNeu = true; % use neurons with significant anova freq effects
useAllNeu = false; % use all available neurons for modeling
includeErrBar = true; %include error bars in plots
% cellType = ["PV","SOM","EXC"]; % order for plotting
cellType = "EXC"; % order for plotting
selType = ["FreqMod","AmpMod","DualMod","InterMod"]; 

% is.allMod = (is.FreqMod | is.DualMod | is.InterMod);
% selType = "allMod";
% is.all = true(nCellTotal,1);
% selType = "all";
neuStepSize = 25; % num neurons to add to model each step
nCellDraws = 40; % over-estimation to poulate arrays
maxDrawEXC = 1000;
nStim = nAmp * nFreq;

% CREATE OUTPUT ARRAYS
PCTall = zeros(nIter,length(selType),nCellDraws,length(cellType)); % (nj,nt,nd,nc)
MDLall = cell(nIter,length(selType),nCellDraws,length(cellType)); % (nj,nt,nd,nc)

PC = struct();
for ns = 1 : length(selType)
    PC.(selType(ns)) = struct();
end

for nc = 1:length(cellType) %nt = 1 : length(selType)
disp(cellType(nc));

    for nt = 1 : length(selType) %nd = 1 : length(nCellDraws)
        fprintf(strcat(selType(nt)," "))
        %disp(strcat("Draw ",num2str(nCellDraws(nd))))
    
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDrawEXC;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
        end
        if useFreqEffNeu == false
            maxNeuAvail = 120;
        end
        nCellDraws = kComp:neuStepSize:maxNeuAvail;
        trainTrial = round(minTrialNum*trainPortion);
        testTrial = round(minTrialNum*(1-trainPortion));
        trainFreq = repmat(freqVal,trainTrial,1);
        testFreq = repmat(freqVal,testTrial,1);
        trainAmp = repmat(ampValAll,trainTrial,1);
        testAmp = repmat(ampValAll,testTrial,1);
        trainStim = repmat((1:25)',trainTrial,1);
        testStim = repmat((1:25)',testTrial,1);
        
        for nd = 1 : length(nCellDraws) %nc = 1:length(cellType)
            %fprintf(strcat(" ",cellType(nc)))
            fprintf(strcat(num2str(nCellDraws(nd))," "));
            nCellDraw = nCellDraws(nd);
    
            % select cell type
            if useFreqEffNeu == true
                if multiAmp == true
                    idx = find(Tc.identity == cellType(nc) & is.(selType(nt)));
                elseif multiAmp == false
                    idx = find(Tc.identity == cellType(nc) & Tc.selectivityPVal < 0.05);
                end
            elseif useFreqEffNeu == false
                idx = find(Tc.identity == cellType(nc));
            end
            
            % reiterate random draws
            for nj = 1:nIter
                
                % randomly draw cells without replacement
                if useAllNeu == true
                    selectedCells = idx;
                    nCellDraw = length(idx);
                elseif useAllNeu == false
                    selectedCells = randsample(idx, nCellDraw, false);
                end
               
                % pre-allocate training and testing data array
                trainResp = zeros(nStim,trainTrial,nCellDraw);
                testResp = zeros(nStim,testTrial,nCellDraw);
        
                % loop over selected cells
                for nk = 1:length(selectedCells)
                    cellResp = Tc.responses{selectedCells(nk)};
                    cellFreq = Tc.frequencies{selectedCells(nk)};
                    cellAmp = Tc.amplitudes{selectedCells(nk)};
                              
                    for nf = 1 : nFreq
                        for na = 1 : nAmp
                        % randomly draw trials at each frequency, no replacement
                        trialNum = find(cellFreq == freqVal(nf) & cellAmp == ampValAll(na));
                        selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
                        stimNum = ((nf*5)-5)+(na);
                        % populate training and testing arrays
                        trainResp(stimNum,:,nk) = cellResp(selectedTrials(1:trainTrial));
                        testResp(stimNum,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
                        end
                    end
                end
        
                % Reshape matricies
                trainResp = reshape(trainResp, nStim*trainTrial, nCellDraw);
                testResp = reshape(testResp, nStim*testTrial, nCellDraw);
        
                % Z-Scoring 
                for nz = 1 : nCellDraw
                    mu = mean(trainResp(:,nz));
                    sigma = std(trainResp(:,nz));
                    trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
                    testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
                end
        
                % Training
                % WITH PCA 
                [coeff,score,latent] = pca(trainResp);
                PC.(selType(nt)).coeff{nj} = coeff;
                PC.(selType(nt)).score{nj} = score;
                PC.(selType(nt)).latent{nj} = latent;
                totalVariance = sum(latent);
                CVA = cumsum((latent./totalVariance).*100); %cumulative variance explained


                % kComp = floor(trainTrial*nFreq*0.5);
                trainResp_red = score(:,kStart:kComp);
                mdl = fitcdiscr(trainResp_red,trainStim,'DiscrimType','linear');
                predStim = predict(mdl,testResp*coeff(:,kStart:kComp));
        
                % collect output
                PCTall(nj,nt,nd,nc) = mean(predStim == testStim); %accuracy
                MDLall{nj,nt,nd,nc} = mdl; % model params
            end 
        end
        fprintf(newline)
    end

end



%% PLOT ACCURACY (ALL EXC 25:25:1000)

%neuStepSize = 5;
%maxDrawEXC = 1000;
figure('Color',[1 1 1],'Theme','Light'); hold on;

% fake data for correct legend order
plot([-1 -2],[-1 -1],'k--');  % Inter & Dual
plot([-1 -2],[-1 -1],'k:'); % Freq & Amp
plot([-1 -2],[-1 -1],'k-.');  % Inter & Dual
plot([-1 -2],[-1 -1],'k-'); % Freq & Amp
% fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.PV,'FaceAlpha',0.3,'EdgeColor','none');
% fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.SOM,'FaceAlpha',0.3,'EdgeColor','none');
% fill([-1 -2 -2 -1],[-1.1 -1.1 -.9 -.9],colors.EXC,'FaceAlpha',0.3,'EdgeColor','none');

% PLOT DATA FROM PCTall2 (FREQ. SEL NEURONS BINNED)
for nt = 1 : length(selType)
    for nc = 1 : length(cellType)
        if strcmp(cellType(nc),'EXC')
            maxNeuAvail = maxDrawEXC;
        else % PV | SOM
            maxNeuAvail = sum(is.(cellType(nc)) & is.(selType(nt)));
        end
        % nCellDraws = 15:neuStepSize:maxNeuAvail;
        target = cellType(nc);
        nBins = 1:length(nCellDraws);
        yMean = mean(squeeze(PCTall(:,nt,nBins,nc)),1).*100;
        % ySems = sem(squeeze(PCTall2(targetAmp,:,nt,nBins,nc)),1).*100;
        yCIs = zeros(length(nBins),2); % confidence intervals (95%)
        for nd = 1 : length(nBins)
            yCIs(nd,:) = mkCI(squeeze(PCTall(:,nt,nd,nc))).*100;
        end
        xLooped = [nCellDraws flip(nCellDraws)]; 
        yLooped = [yCIs(:,1)' flip(yCIs(:,2))'];
        fill(xLooped,yLooped,colors.(target),'FaceAlpha',0.3,'EdgeColor','none');
        switch nt
            case 1
                plot(nCellDraws,yMean,'--','LineWidth',1,'Color',colors.(target));
            case 2
                plot(nCellDraws,yMean,':','LineWidth',1,'Color',colors.(target));    
            case 3
                plot(nCellDraws,yMean,'-.','LineWidth',1,'Color',colors.(target));    
            case 4
                plot(nCellDraws,yMean,'-','LineWidth',1,'Color',colors.(target));    
        end
    end
end


ylabel('model accuracy (%)','FontName','Arial');
xlabel('Neurons used','FontName','Arial');
title('Stimulus decoding performance across al 25 stimulus combinations','FontName','Arial');
subtitle('EXC neurons only, 25 PCs in each model','FontName','Arial');
set(gca,'XLim',[20 1020],'YLim',[0 90],'YTick',0:10:100,'XTick',25:75:1000);
legend({"Freq. Mod","Amp. Mod","Dual Mod.","Interaction"},'Location','Northwest',...
    'NumColumns',1,'box','off');




