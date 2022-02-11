% makeIndividualRGCDensityFiles
% Produces .mat files that contain RGC density data
%
% Description:
%   This routine operated upon an Excel spreadsheet of RGC density data
%   and saves a set of labeled and organized .mat files that contain these
%   data. It should not need to be run in the future, as its only purpose
%   was to create the files that are now stored within this directory.
%

tableFileNames = {'DENSITY5_gc_resorted_computedAverage.xls','curcioReportedGCAverage.xls'};
outputFileNameStem = 'curcioRawRGCDensity_';

for tt = 1:length(tableFileNames)
    T=readtable(tableFileNames{tt});
    
    support=table2array(T(1,4:end));
    
    uniqueSubjectNames = unique(table2array(T(2:end,1)));
    
    meridianNames = {'NASAL','SUPERIOR','TEMPORAL','INFERIOR'};
    
    for nn=1:length(uniqueSubjectNames)
        curcioRawRGCDensity = struct;
        curcioRawRGCDensity.support = support;
        subjectName= uniqueSubjectNames{nn};
        for mm = 1:length(meridianNames)
            thisMeridianName = meridianNames{mm};
            rowIdx=find( strcmp(table2array(T(:,1)),subjectName) .* strcmp(table2array(T(:,2)),thisMeridianName) );
            curcioRawRGCDensity.(lower(thisMeridianName)) = table2array(T(rowIdx,4:end));
        end
        curcioRawRGCDensity.meta.subjectName = subjectName;
        curcioRawRGCDensity.meta.dataTableName = tableFileNames{tt};
        curcioRawRGCDensity.meta.densityUnits = T.Units{rowIdx};
        curcioRawRGCDensity.meta.supportUnits = T.Units{1};
        save([outputFileNameStem subjectName], 'curcioRawRGCDensity');
    end % loop over subjects
    
end % loop over tables