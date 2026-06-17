
% FAS_decodeFreq_bySession_withTitr.m
%
% with titration from all-EXC (no Inh) to max-INH (replacing EXC)
%

% perform population decoding vary 

nIter = 100; % Number of model iterations to run
trainPortion = 0.8; % percentage of trials to be used in training (remaining for testing)
usePCA = false; % use PCA step prior to model training
kStart = 1; % first PC to include
kComp = 15; % total number of PCs to use in model training / testing
neuStepSize = 15; % num neurons to add to model each step
neuTiSteps = 50; % nueron titration steps, nINH +1 max
deltaSteps = 0; %[0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1];
targetAmp = 5; % [0.33 0.76 1.6 3.5 7.7]

% CREATE OUTPUT ARRAYS
PCTall = zeros(nIter,length(sessIDs),neuTiSteps); % (nj,ns,ng,ne,sess vs pseudopop)
MDLall = cell(nIter,length(sessIDs),neuTiSteps); % (nj,ns,ng,nd)
CVEall = zeros(nIter,length(sessIDs),neuTiSteps);
PCFall = cell(nIter,length(sessIDs),neuTiSteps); % (nj,ns,ng,nd)


for ns = 1:length(sessIDs) %nt = 1 : length(selType)
    locSess = char(sessIDs(ns)); % local session name 'm000_00'
    locIdx = is.mouseNum(:,strcmp(mouseNums,locSess(1:4))); % all neurons/rows from this mouse
    disp(strcat("Session ",num2str(ns),": ",locSess));
    % disp(selType(nt));
    if usePCA == true && sum(is.sessID(:,ns)) < kComp
        % fewer neurons than min PC number
        continue;
    end

    % nCellDraws = kComp : neuStepSize : maxNeuAvail;
    trainTrial = round(minTrialNum*trainPortion);
    testTrial = round(minTrialNum*(1-trainPortion));
    trainFreq = repmat(freqVal,trainTrial,1);
    testFreq = repmat(freqVal,testTrial,1);

    % loc neurons in session
    locPVs = find(is.sessID(:,ns) & is.PV); 
    locSOMs = find(is.sessID(:,ns) & is.SOM);
    locEXCs = find(is.sessID(:,ns) & is.EXC);
    % % loc neurons from this mouse 
    % otherPVs = find(locIdx & is.PV);
    % otherSOMs = find(locIdx & is.SOM);
    % otherEXCs = find(locIdx & is.EXC);
    % % all other neurons not in session
    % otherPVs = find(~is.sessID(:,ns) & is.PV);
    % otherSOMs = find(~is.sessID(:,ns) & is.SOM);
    % otherEXCs = find(~is.sessID(:,ns) & is.EXC);

    nEXC = length(locEXCs);
    nPV = length(locPVs);
    nSOM = length(locSOMs);
    nINH = max([nPV nSOM]); 
    binsEXC = nEXC:-1:nEXC-nINH;
    binsPV = 0:1:nPV;
    binsSOM = 0:1:nSOM;

    if nEXC == 0
        continue;
    end

    for nt = 1 : length(binsEXC)

        locnEXC = binsEXC(nt); % number of EXC neurons to use
        if nPV > 0
            locnPV = binsPV(nt); % number of PV neurons to use
        else
            locnPV = 0;
        end
        if nSOM > 0
            locnSOM = binsSOM(nt); % number of PV neurons to use
        else
            locnSOM = 0;
        end

        % reiterate random draws
        for nj = 1:nIter
            
            % build arrays: pseudo populations to match sess cells
            selectedPVs = randsample(locPVs, locnPV, false); % random draw from population
            selectedSOMs = randsample(locSOMs, locnSOM, false); % random draw from population
            selectedEXCs = randsample(locEXCs, locnEXC, false); % random draw from population
            selectedNeus = cat(1,selectedPVs,selectedSOMs,selectedEXCs); % random sampled pop
            numNeuSess = size(selectedNeus,1);
           
            % pre-allocate training and testing data array
            trainResp = zeros(nFreq,trainTrial,numNeuSess); % session
            testResp = zeros(nFreq,testTrial,numNeuSess);
    
            % loop over selected cells
            for nk = 1:length(selectedNeus)
    
                % FROM THIS SESSION
                cellResp = Tc.responses{selectedNeus(nk)};
                cellFreq = Tc.frequencies{selectedNeus(nk)};
                if multiAmp
                    cellAmp = Tc.amplitudes{selectedNeus(nk)};
                end
                          
                for nl = 1:nFreq
                    % randomly draw trials at each frequency, no replacement
                    if multiAmp
                        trialNum = find(cellFreq == freqVal(nl) & cellAmp == ampValAll(targetAmp));
                    else
                        trialNum = find(cellFreq == freqVal(nl));
                    end
                    % populate training and testing arrays
                    selectedTrials = randsample(trialNum, trainTrial+testTrial, false);
                    trainResp(nl,:,nk) = cellResp(selectedTrials(1:trainTrial));
                    testResp(nl,:,nk) = cellResp(selectedTrials(trainTrial+1:end));
                end
    
            end
    
            % Reshape matricies
            trainResp = reshape(trainResp, nFreq*trainTrial, numNeuSess);
            testResp = reshape(testResp, nFreq*testTrial, numNeuSess);
    
            % Z-Scoring 
            for nz = 1 : numNeuSess
                % session
                mu = mean(trainResp(:,nz));
                sigma = std(trainResp(:,nz));
                trainResp(:,nz) = (trainResp(:,nz) - mu) ./ sigma;
                testResp(:,nz) = (testResp(:,nz) - mu) ./ sigma;
            end
    
            % Percent variance explained
            [coeffS,scoreS,latentS] = pca(trainResp);
            % PC.(cellType(nc)).coeff{ni} = coeff;
            % PC.(cellType(nc)).score{ni} = score;
            % PC.(cellType(nc)).latent{ni} = latent;
            explainedS = 100 * latentS / sum(latentS);
            cumFracS = cumsum(sort(explainedS,'descend'))/sum(explainedS);
            % CVEall(nj,ns,ng,ne) = sum(explained(1:15));
            CVEall(nj,ns,nt,1) = find(cumFracS >= 0.80, 1); %num of PCs to explain 80% of variance (sess)
            % --------------------
            % [coeffP,scoreP,latentP] = pca(trainRespP);
            % explainedP = 100 * latentP / sum(latentP);
            % cumFracP = cumsum(sort(explainedP,'descend'))/sum(explainedP);
            % CVEall(nj,ns,nt,1,2) = find(cumFracP >= 0.80, 1); %num of PCs to explain 80% of variance (pseudopop)
            % ---------------------
            % % Bar plot + cumulative line
            % figure;
            % bar(explained,'FaceColor',[0.6 0.6 0.9]); hold on;
            % plot(cumsum(explained),'r-o','LineWidth',1.5);
            % xlabel('Principal Component');
            % ylabel('Variance Explained (%)');
            % legend('Individual','Cumulative','Location','best');
            % xlim([0.5 numel(explained)+0.5]);
            % grid on;
            % hold off;
    
            % TRAIN WITH PCA
            if usePCA == true
                % Session data
                trainResp_red = scoreS(:,kStart:kComp);
                PCtunings = zeros(length(freqVal),kComp);
                PCvariation = zeros(length(freqVal),kComp);
                for nf = 1 : length(freqVal)
                    PCtunings(nf,kStart:kComp) = mean(trainResp_red(trainFreq == freqVal(nf),:),1);
                    PCvariation(nf,kStart:kComp) = std(trainResp_red(trainFreq == freqVal(nf),:),1);
                end
                for ng = 1 : length(gammaSteps)
                    for ne = 1 : length(deltaSteps)
                        mdl = fitcdiscr(trainResp_red,trainFreq,'DiscrimType','linear','Gamma',gammaSteps(ng),...
                            'Delta',deltaSteps(ne));
                        predFreq = predict(mdl,testResp*coeffS(:,kStart:kComp));
                        % collect output
                        PCTall(nj,ns,nt,ne) = mean(predFreq == testFreq); % accuracy
                        MDLall{nj,ns,nt,ne} = mdl; % model parameters
                        % PCFall{nj,ns,ng,ne} = cat(3, PCtunings, PCvariation); % accuracy
                        PCFall{nj,ns,nt,ne} = coeffS; % coeffs from PCA
                    end
                end
    
            % TRAIN WITHOUT PCA
            elseif usePCA == false
                for ne = 1 : length(deltaSteps)
                    % Session
                    mdl = fitcdiscr(trainResp,trainFreq,'DiscrimType','linear',...
                        'Delta',deltaSteps(ne));
                    predFreq = predict(mdl,testResp);
                    PCTall(nj,ns,nt) = mean(predFreq == testFreq); % accuracy
                    MDLall{nj,ns,nt} = mdl; % model parameters
                end
            end
        end 

    end % nt number titration bins
    fprintf(newline)
end


% allPVals = zeros(ns,1);
% allOS = zeros(ns,1); % observed stats
% for n = 1 : size(PCTall,2)
%     array1 = PCTall(:,n,1,1,1);
%     array2 = PCTall(:,n,1,1,2);
%     [pVal, obsStat, permStats] = permTest2sample(array1, array2, 5000);
%     allPVals(n) = pVal;
%     allOS(n) = obsStat;
% end
% 
% PCT = squeeze(mean(PCTall,1));
% sig1 = allPVals < pValThresh & allOS < 0;
% sig2 = allPVals < pValThresh & allOS > 0;
% 
% if usePCA
%     wPCA = struct();
%     wPCA.PCTall = PCTall;
%     wPCA.MDLall = MDLall;
%     wPCA.allPVals = allPVals;
%     wPCA.allOS = allOS;
% else
%     woPCA = struct();
%     woPCA.PCTall = PCTall;
%     woPCA.MDLall = MDLall;
%     woPCA.allPVals = allPVals;
%     woPCA.allOS = allOS;
% end



%% Plotting

%How do models change when titrating in INH neurons: VIEW ALL
figure('Color',[1 1 1]); hold on;
for ns = 1 : size(PCTall,2)
    locAr = squeeze(mean(PCTall(:,ns,:),1));
    locAr = locAr(locAr ~= 0) .* 100;
    plot(0:length(locAr)-1,locAr);
end
xlabel('Number of INH neurons included','FontName','Arial');
ylabel('Frequency decoding performance (per session)','FontName','Arial');


% Find sessions with sig diff between noINH and maxINH
allPval = zeros(size(PCTall,2),1);
allObs = zeros(size(PCTall,2),1);
sessType = repmat("",size(PCTall,2),1);
figure('Color',[1 1 1]); 
tiledlayout(6,6,"TileSpacing","compact","Padding","compact");
for ns = 1 : size(PCTall,2)
    locMat0 = squeeze(PCTall(:,ns,:));
    locMat = locMat0(:,~all(locMat0 == 0,1)); % remove zero columns
    % % Permutation test
    % [allPval(ns), allObs(ns), permStats] = permTest2sample(locMat(:,1),locMat(:,end),1000);
    % Ranksum test
    allPval(ns) = ranksum(locMat(:,1),locMat(:,end));
    allObs(ns) = mean(locMat(:,end)) - mean(locMat(:,1));

    % sessType
    if any(is.sessID(:,ns) & is.EXC) && any(is.sessID(:,ns) & is.PV)
        sessType(ns) = "PV";
    elseif any(is.sessID(:,ns) & is.EXC) && any(is.sessID(:,ns) & is.SOM)
        sessType(ns) = "SOM";
    end
    % ploting swarm
    if allPval(ns) < 0.05
        nexttile(); hold on;
        % swarmchart(ones(size(PCTall,1),1),locMat(:,1).*100,'filled','o','MarkerFaceColor',colors.EXC);
        locCI1 = mkCI(locMat(:,1).*100);
        locMean1 = mean(locMat(:,1).*100);
        errorbar(1,locMean1,locCI1(2)-locMean1,'Color',colors.EXC);
        plot(1,locMean1,'.','MarkerSize',10,'Color',colors.EXC);
        if sessType(ns) == "PV"
            % swarmchart(ones(size(PCTall,1),1)+2,locMat(:,end).*100,'filled','o','MarkerFaceColor',colors.PV);
            locCI2 = mkCI(locMat(:,end).*100);
            locMean2 = mean(locMat(:,end).*100);
            errorbar(2,locMean2,locCI2(2)-locMean2,'Color',colors.PV);
            plot(2,locMean2,'.','MarkerSize',10,'Color',colors.PV);
        elseif sessType(ns) == "SOM"
            % swarmchart(ones(size(PCTall,1),1)+2,locMat(:,end).*100,'filled','o','MarkerFaceColor',colors.SOM);
            locCI2 = mkCI(locMat(:,end).*100);
            locMean2 = mean(locMat(:,end).*100);
            errorbar(2,locMean2,locCI2(2)-locMean2,'Color',colors.SOM);
            plot(2,locMean2,'.','MarkerSize',10,'Color',colors.SOM);
        end
        hold off;
        set(gca,'XTick',[],'XLim',[0.5 2.5]);
    end
end

% allObs(allPval<0.05)
% sessType(allPval<0.05)


%% Relationship: Inh Nueron Addition & Frequency Decoding Performance
% Only from sessions with high performance (top 50%) and many INH neurons (top50%)

accByNeus = squeeze(mean(PCTall,1)); % sessions by inh neuron additions
numInhNeuPerSess = sum(accByNeus ~= 0,2); % number of inhibitory neuron bins per session
inhNeuThresh = median(numInhNeuPerSess); % threshold for number of inhibitory neurons to include (top half)
accThresh = median(accByNeus(:,1)); % threshold of model accuracy to include 

% grab top half performing models & top half of sessions with most inh neurons
topPerfSess = accByNeus(:,1) > accThresh; % logical, 164x1 session(s) performance > threshjold
topInhSess = numInhNeuPerSess > inhNeuThresh; % logical 164x1 sessions(s) with many inh neurons
sessList = find(topPerfSess & topInhSess);

nMdlItr = size(PCTall,1); % number of model iterations calculated above
allMdls = cell(length(sessList),1);
slopes = [];
for ns = 1 : length(sessList)
    locNum = numInhNeuPerSess(sessList(ns));
    locX = repmat(1:locNum,nMdlItr,1); % inh neuron bins by model train repetitions
    locY = squeeze(PCTall(:,sessList(ns),1:locNum));
    mdl = fitlm(locX(:),locY(:)); % linear fit to the data, inh neuron bins (x-axis) and frequency decoding performance (y-axis)
    allMdls{ns} = mdl;
    % disp(strcat("sess ",num2str(sessList(ns))," p-val = ",num2str(mdl.ModelFitVsNullModel.Pvalue)));
    if mdl.ModelFitVsNullModel.Pvalue < pValThresh
        slopes = [slopes mdl.Coefficients.Estimate(2)];
    end
end

ns = 11;
predict(mdl,[1 locNum]')