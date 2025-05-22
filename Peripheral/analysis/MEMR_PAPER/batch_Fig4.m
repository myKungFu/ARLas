function [] = batch_Fig4()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% batchMEMR_groupData
% 
% Batch analysis of group MEMR data 
%
% Author: Shawn Goodman & Ehsan Khalili
% Date: January 23, 2024
% Last Updated: January 23, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    parentDrive = 'D'; % use 'C' or 'D'
    dataPathName = [parentDrive,':\MEMR_Analysis1\']; % location of raw data
    savePath = [parentDrive,':\check\']; % where to save analyzed data

    d = dir([dataPathName,'*.mat']); % read in file directory information
    nFiles = size(d,1); % number of subject folders
    counter = 1;
        for ii=1:nFiles
            close all 
        disp(['Analyzing file ',num2str(ii),' of ',num2str(nFiles)])
        fileName = d(ii).name;
        runNumber = str2num(fileName(10));
        dummy = load([dataPathName,fileName]); 
        MEMR_mem = dummy.MEMR_mem;
        FN = [dataPathName,fileName];
        MEMR_Fig4_new(FN)

 FolderName = 'D:\check\';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = [fileName,num2str(get(FigHandle, 'Number'))];
  set(0, 'CurrentFigure', FigHandle);
  saveas(FigHandle,fullfile(FolderName, [FigName '.png'])); %Specify format for the figure
end
   

end