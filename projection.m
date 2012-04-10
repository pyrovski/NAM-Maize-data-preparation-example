% SNP prjection
%
% Liya Wang, 03/15/10

%% rewrite map file
% replace this with tail -n +2 | cut -f1,3-|sed -re 's/m(0)*//'
%{
 fid = fopen('NAM_Map_20090730.txt');
 tline = fgetl(fid); %remove top line
 fid1 = fopen('newmap.txt', 'w');
 
 while 1
     tline = fgetl(fid);
     if ~ischar(tline), break, end
     [chr r] = strtok(tline);
     [marker r] = strtok(r);
     [marker r] = strtok(r);
     [pos r] = strtok(r);
     [strpos r] = strtok(r);
     fprintf(fid1, '%3d%7d%7.1f%11d\n', str2num(chr), str2num(marker(2:end)), str2num(pos), str2num(strpos));
 end
 fclose(fid);
 fclose(fid1);
%}

%% rewrite phen.txt and imputedMarkerchr10.txt

%% rewrite fastphase file
% this seems to be saving only fields 3,4,12-
% why not replace with with cut -f3,4,12- ?
%rs	alleles	chrom	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCcode	t0	t1	t2	t3	t4	t5	t6	t7	t8	t9	t10	t11	t12	t13	t14	t15	t16	t17	t18	t19	t20	t21	t22	t23	t24	t25	t26
%{
 fid = fopen('fastphase_chr10.txt');
 tline = fgetl(fid);
 fid1 = fopen('newfast.txt', 'w');
 while 1
     tline = fgetl(fid);
     if ~ischar(tline), break, end
     [rs r] = strtok(tline);
     [alleles r] = strtok(r);
     [chrom r] = strtok(r);
     [pos r] = strtok(r);
     [strand r] = strtok(r);
     [assembly r] = strtok(r);
     [center r] = strtok(r);
     [protLSID r] = strtok(r);
     [assayLSID r] = strtok(r);
     [panelLSID r] = strtok(r);
     [QCcode r] = strtok(r);
     fprintf(fid1, '%s %s %s\n', chrom, pos, r);
 end
 fclose(fid),
 fclose(fid1);
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

%marker = importdata('imputedMarkersGWAS.chr10.082809.txt');
%marker = double(marker.data); %
%{
split first columnt into population number and entry number (sample)
see https://pods.iplantcollaborative.org/wiki/display/ipg2p/GLM+Report
cat imputedMarkersGWAS.chr10.082809.txt|tail -n +2|sed -re 's/Z([[:digit:]]+)E([[:digit:]]+)/\1\t\2/' > imputedMarkers.chr10.082809.txt
%}

marker = load('imputedMarkerchr10.txt');
newmap = load('newmap.txt');
newfast = load('newfast.txt'); % from fastphase_chr10.txt
[p q] = size(newfast);
phen = load('phen.txt');
[m n] = size(phen);
inputHigh = size(newfast,1);
%inputHigh = 7500;
inputLow = 1;
%inputLow = 2501;
numInputs = inputHigh - inputLow + 1;
projectedSNP = fopen('projectedSNP.dat', 'w');
%projectedSNP = zeros(m, numInputs + 2);
tic;
fasts = [newfast(inputLow:inputHigh, 2)];
li = zeros(numInputs, 1);
rightpos = zeros(numInputs, 1);
leftpos = zeros(numInputs, 1);
rightmark = zeros(numInputs, 1);
leftmark = zeros(numInputs, 1);
pd = zeros(numInputs, 1);
for sj = 1:numInputs
    j = sj + inputLow - 1;
    ri = find(newmap(:, 4)<fasts(sj));
    li(sj) = length(ri);
    if(li(sj) > 0)
        if(li(sj) == size(newmap,1))
            rightpos(sj) = 148000162; leftpos(sj) = newmap(end, 4);
            rightmark(sj)= 1106; leftmark(sj) = newmap(1, 2);
        else
            rightpos(sj) = newmap(ri(end)+1, 4); leftpos(sj) = newmap(ri(end), 4);
            rightmark(sj) = newmap(ri(end)+1, 2); leftmark(sj) = newmap(ri(end), 2);
        end
    else
        rightpos(sj) = newmap(1, 4); leftpos(sj) = 0;
        rightmark(sj)= newmap(1, 2); leftmark(sj)=1029;  % use t1 marker
    end
end
pd = (fasts - leftpos) ./ (rightpos - leftpos);
newRow = zeros(1, numInputs + 2);
for i = 1 : m
    if mod(i,1000) == 0
        i
    end
    pop = phen(i, 1); %population
    sam = phen(i, 2); %sample
    pheno = phen(i, end); %Chromosome 10
    %newRow(end - 1) = pheno;
    %newRow(end) = pop;

    % fmark is the only part of this loop that depends on i
    fmark = marker(find(marker(:, 1) == pop & marker(:, 2) == sam), 1:end); %t1 m1030-m1106 t2 -- marker value
    %! todo: if fmark is empty, the loop below could probably be avoided
    % question: m1030 is not in map file, use t1 and m1031 instead? m1106 is
    % not in map file either, use m1105 and m1106 instead?
    %fasts = [newfast(:, 2) newfast(:, pop + 3)]; %snp.pos + parent snp value
    %projectedSNP(i,1:(end-2)) = newfast(inputLow:inputHigh, pop + 3)';
    %{
    This could be cached, as the number of populations is limited
    %}
    selj = find(newfast(inputLow:inputHigh, pop + 3) > 0);
    
    %    projectedSNP(i, selj) = fmark(leftmark(selj) - 1026)' .* (1 - pd(selj)) + ...
    %        fmark(rightmark(selj) - 1026)' .* pd(selj);;

    % loop is faster for selj
    for sj = 1:length(selj)
        j = selj(sj);


        newRow(j) = fmark(leftmark(j) - 1026) * (1 - pd(j)) + fmark(rightmark(j) - 1026) * pd(j);
    end
    fwrite(projectedSNP, newRow, 'double');
end
toc
%save projectedSNP projectedSNP
%{
We can calculate the rank of projectedSNP by calculating the rank of 
projectedSNP * projectedSNP'.  As nrows << ncols, this should result in a 
shorter computation time.
tic;
rank(single(projectedSNP * projectedSNP'))
toc
%}
