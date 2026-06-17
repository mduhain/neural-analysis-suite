


cellTypes = ["EXC","PV","SOM"];

for nt = 1 : length(cellTypes)

    figure('Color',[1 1 1],'Position',[25 25 867 925]);
    hold on;
    tiledlayout(1,nAmp);
    target = cellTypes(nt);
    targetNeurons = find(is.(target));
    figureTitle = strcat("Frequency Tuning by Amplitude Heatmap (",target,")");
    
    for na0 = 1 : nAmp
    
        targetAmp = na0; %starting amplitude [7.7, 3.5, 1.6, 0.77, 0.33] um
    
        % BUILD MATRIX OF RESPONSE VALUES FROM TARGET NEURONS
        respMat = zeros(num.(target),numFreqs,nAmp); %pre-allocate
        for n = 1 : length(targetNeurons)
            for na = 1 : nAmp
                %check if this neuron and this amplitude have freq selectivity
                if M(targetNeurons(n)).ampGroupedAnova(na) == 1
                    respMat(n,:,na) = M(targetNeurons(n)).avgRespPost(na,:);
                else
                    respMat(n,:,na) = nan(1,5);
                end
            end
        end
        
        % FIND STARTING VALUES OF TARGET AMP
        goodRows = all(~isnan(respMat(:,:,targetAmp)),2); %find rows in RespMat that are not all nans
        locNeuronList = targetNeurons(goodRows); %list of local neuron numbers [M|R] 
        [~,fIndx] = max(respMat(goodRows,:,targetAmp),[],2); %find top frequency, per neuron
        [bestFreq,fIndxSorted] = sort(fIndx,'ascend'); % Sort by preferred frequency
        locMat = zeros(length(fIndxSorted),numFreqs); %matrix of normalized (1-0) responses
        selVal = zeros(length(fIndxSorted),1); % array of selectivity values
        neuNums = zeros(length(fIndxSorted),1); % array of neuron numbers [M|R]
        for nn = 1 : length(fIndxSorted) %Loop all neurons in amp group
            neuNums(nn) = locNeuronList(fIndxSorted(nn)); % neuron index
            selVal(nn) = quickSelectivity(M(neuNums(nn)).avgRespPost(targetAmp,:)); %neuron freq selectivity metric
            locMat(nn,:) = rescale(M(neuNums(nn)).avgRespPost(targetAmp,:)); %neurons z-scored response
        end
        
        % ADDITIONAL SORT BY SELECTIVITY METRIC
        neuNumsSelSort = neuNums;
        %locMat2 = zeros(size(locMat));
        for nf = 1 : numFreqs
            locTargets = bestFreq == nf; %Nx1 binary
            [~,locIndx] = sort(selVal(locTargets),'ascend'); %sort frequency group by selectivity
            tempMat = locMat(locTargets,:); %extract freq group response matrix
            locMat(locTargets,:) = tempMat(locIndx,:); %apply sort
            locNeurons = neuNums(locTargets); %copy local neuron numbers within freq tuning group
            neuNumsSelSort(locTargets) = locNeurons(locIndx); %apply selectivity sort to neurons numbers
        end
    
        %PLOTTING SECTION
        nexttile;
        imagesc(locMat); 
        subtitle(strcat(num2str(amps(targetAmp))," \mum"),'FontName','Arial');
        set(gca,'XTickLabel',freqs,'XTickLabelRotation',90);
        if na0 == 1
            ylabel('Neuron Number','FontName','Arial'); 
        end
        if na0 == 3
            xlabel('Frequency (Hz)','FontName','Arial'); 
            title(figureTitle,'FontName','Arial');
        end
        
    end

end
