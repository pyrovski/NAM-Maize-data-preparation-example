function [skip] = projection(imputedMarkerFilename, mapFilename, fastphaseFilename, ...
    phenoFilename, phenoCol, outputFilename)
% SNP prjection
%
% Liya Wang, 03/15/10
% Peter Bailey, 4/2012

%for data preparations, see prepare.sh

%! @todo the integer constants in this file are suspect; they work for
%    chromosome 10, but no others.
%see https://pods.iplantcollaborative.org/wiki/display/ipg2p/GLM+Report
%}

%% Do projection
% To project a SNP, use its physical position, snp.pos, to find the flanking NAM markers, based on agp_pos from the NAM map.
%  For each NAM line in the residuals file, find its projected SNP value as follows:
%
%    1. Determine the parent number for the line as the sample Z-number or 17, if the line name starts with M.
%    2. Find the parent SNP value for that line from the fastphase file.
%    3. for each position (snp.pos in fastphas), Determine the proportional distance, pd, of SNP from its left and right flanking markers as
%            pd = (snp.pos - left marker agp_pos) / (right marker agp_pos � left marker agp_pos).
%    4. If the parent value = 0, set the snp.value = 0 for that line.
%    5. If the parent value = 1, set snp.value = left marker value * (1 - pd) + right marker value * pd
%    6. Note that snp.pos comes from the fastphase file, the flanking markers and agp_pos are determined from the NAM map,
%            and the marker values come from the imputed marker file.

tStart = tic;
marker = load(imputedMarkerFilename);
newmap = load(mapFilename);
newmap(:, 3) = [];

%for whatever reason, uint64s are many times slower than doubles in matlab
%for various operations.
%newmap = uint64(newmap);
newfast = load(fastphaseFilename); % from fastphase_chr10.txt
phen = load(phenoFilename);

tReadStop = toc(tStart)
tCompStart = tic;

[m n] = size(phen);
markerLower = newmap(1, 2);
markerUpper = newmap(end, 2);
numMarkers = size(newmap, 1);
%    inputHigh = size(newfast,1);
inputHigh = find(newfast(:, 2) < newmap(end, 3), 1, 'last');
%    inputLow = 1;
inputLow = find(newfast(:, 2) > newmap(1, 3), 1, 'first');
numInputs = inputHigh - inputLow + 1; % number of SNPs
projectedSNP = fopen(outputFilename, 'w');
posLower = 1;
posUpper = newfast(end, 2);
fasts = newfast(inputLow:inputHigh, 2);
rightpos = zeros(numInputs, 1);
leftpos = zeros(numInputs, 1);
rightmark = zeros(numInputs, 1);
leftmark = zeros(numInputs, 1);
%skip = false(numInputs, 1);

snpIndex = 1;
snpPos = fasts(snpIndex);

for interval = 1:(numMarkers - 1)
%      if(interval >= 35 && interval <= 37)
%          disp('stop here');
%      end
    li = interval;
    ri = interval + 1;

%note that the interval is close on the left and open on the right
    while(snpPos >= newmap(li, 3) && snpPos < newmap(ri, 3))
        leftmark(snpIndex) = newmap(li, 2);
        rightmark(snpIndex) = newmap(ri, 2);
        leftpos(snpIndex) = newmap(li, 3);
        rightpos(snpIndex) = newmap(ri, 3);

        snpIndex = snpIndex + 1;
        if(snpIndex <= numInputs)
            snpPos = fasts(snpIndex);
        else
            snpPos = 0;
            break;
        end
    end
    
end

%li = 1;
%ri = 2;

%for sj = 1:numInputs
    %j = sj + inputLow - 1;
    
    % @todo the map should be sorted, so we can take this find() out of
    % the loop, and shouldn't need to use the skip array
    
    % find marker to the left of this SNP
    %li = find(newmap(:, 3) < fasts(sj), 1, 'last');
    
    % find marker to the right of this SNP
    %ri = find(fasts(sj) < newmap(:, 3), 1, 'first');
    
%     if fasts(sj) > ri
%         ri = ri + 1;
%         li = li + 1;
%     end
    % @todo discard positions with only one flanking marker
    %if(length(li) > 0)
        % at least one marker has a lower position than this SNP
        %if(length(ri) > 0)
            % this SNP is between markers
%             rightpos(sj) = newmap(ri, 3); leftpos(sj) = newmap(li, 3);
%             rightmark(sj) = newmap(ri, 2); leftmark(sj) = newmap(li, 2);
        %else
            % all markers have a lower position than this SNP
            %rightpos(sj) = posUpper; leftpos(sj) = newmap(end, 3);
            %rightmark(sj)= markerUpper; leftmark(sj) = newmap(1, 2);
            %skip(sj) = true;
            %disp('should never get here...');
        %end
    %else
        % all markers have a higher position than this SNP
        %rightpos(sj) = newmap(1, 3); leftpos(sj) = posLower;
        %rightmark(sj)= newmap(1, 2); leftmark(sj)=markerLower;  % use t1 marker
        %skip(sj) = true;
            %disp('should never get here...');
    %end
% end
pd = double(fasts - leftpos) ./ double(rightpos - leftpos);


if ~isempty(find(pd < 0 | pd > 1, 1))
    disp('pd error!')
    return
end
tPDStop = toc(tCompStart)

tPopCacheStart = tic;
newRow = zeros(1, numInputs);
popCount = length(unique(phen(:,1)));
popCache = cell(1, popCount);
for pop = 1:popCount
    popCache{pop} = find(newfast(inputLow:inputHigh, pop + 3) > 0);
end

tPopCacheStop = toc(tPopCacheStart)
tProjectStart = tic;

lastPop = 0;
for i = 1 : m
    if mod(i,1000) == 0
        i
    end
    pop = phen(i, 1); %population
    sam = phen(i, 2); %sample
    %        newRow(end - 1) = phen(i, phenoCol); % Chromosome
    %        newRow(end) = pop;
    
    % fmark is the only part of this loop that depends on i
    fmark = marker(find(marker(:, 1) == pop & marker(:, 2) == ...
        sam), 3:end);
    %t1 m1030-m1106 t2 -- marker value
    % question: m1030 is not in map file, use t1 and m1031 instead? m1106 is
    % not in map file either, use m1105 and m1106 instead?
    
    if(lastPop ~= pop)
        selj = popCache{pop};
    end
    
    %    projectedSNP(i, selj) = fmark(leftmark(selj) - 1026)' .* (1 - pd(selj)) + ...
    %        fmark(rightmark(selj) - 1026)' .* pd(selj);;
    
    % loop is faster for selj
    for sj = 1:length(selj)
        %if(skip(sj))
        %continue;
        %end
        j = selj(sj);
        newRow(j) = fmark(leftmark(j) - markerLower + 1) * ...
            (1 - pd(j)) + ...
            fmark(rightmark(j) - markerLower + 1) * pd(j);
    end
    fwrite(projectedSNP, newRow, 'double');
    lastPop = pop;
end
elapsedComp = toc(tProjectStart)
elapsed = toc(tStart);
output = sprintf('projected %d SNPs for %d individuals in %f s: %f ind/s, %f B/s out', m, ...
    numInputs, elapsed, numInputs / elapsed, numInputs ...
    * m * 8 / elapsed);
disp(output);



end

