for i = 1:10
    iStr = int2str(i);
    markerFilename = ['imputedMarkers.chr' iStr '.merged']
    mapFilename = ['map.' iStr '.txt']
    phenoFilename = ['fastphase_chr' iStr '.txt.filtered']
    residFilename = 'residuals_distancebased_for_asi.txt.filtered'
    outFilename = ['projectedSNP.' sprintf('%04d', i) '.dat']
    projection(markerFilename, mapFilename, ...
    phenoFilename, ...
    residFilename, 2 + i, ...
    outFilename);
end
