% makeIndividualConeDensityFiles
% Produces .mat files that contain cone density data
%
% Description:
%   This routine operated upon an Excel spreadsheet of cone density data
%   and saves a set of labeled and organized .mat files that contain these
%   data. It should not need to be run in the future, as its only purpose
%   was to create the files that are now stored within this directory.
%

tableFileNames = {'DENSITY8_cones_resorted_computedAverage.xlsx','curcioReportedConeAverage.xlsx'};
outputFileNameStem = 'curcioRawConeDensity_';

for tt = 1:length(tableFileNames)
    T=readtable(tableFileNames{tt});
    
    support=table2array(T(1,4:end));
    
    uniqueSubjectNames = unique(table2array(T(2:end,1)));
    
    meridianNames = {'NASAL','SUPERIOR','TEMPORAL','INFERIOR'};
    
    for nn=1:length(uniqueSubjectNames)
        curcioRawConeDensity = struct;
        curcioRawConeDensity.support = support;
        subjectName=uniqueSubjectNames{nn};
        for mm = 1:length(meridianNames)
            thisMeridianName = meridianNames{mm};
            rowIdx=find( strcmp(table2array(T(:,1)),subjectName) .* strcmp(table2array(T(:,2)),thisMeridianName) );
            curcioRawConeDensity.(lower(thisMeridianName)) = table2array(T(rowIdx,4:end));
        end
        curcioRawConeDensity.meta.subjectName = subjectName;
        curcioRawConeDensity.meta.dataTableName = tableFileNames{tt};
        curcioRawConeDensity.meta.densityUnits = T.Units{rowIdx};
        curcioRawConeDensity.meta.supportUnits = T.Units{1};
        save([outputFileNameStem subjectName], 'curcioRawConeDensity');
    end % loop over subjects
    
end % loop over tables