function [] = ARLasVersionInfo()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ARLasVersionInfo
%
% Contains a list of the updates made to ARLas.
% Starting date is November 16, 2016
% 
% Auditory Research Lab, The University of Iowa
% Deptartment of Communication Sciences & Disorders
% The University of Iowa
% Author: Shawn S. Goodman, PhD
% Date: November 16, 2016
% Last Updated: January 7, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2016.11.16
%   Changed handling of initialization defaults.
%       Now uses the desired names and looks for appropriate indices.
%       This change was precipitated by problems with the operating
%       systems frequently changing device parameters. 
%   Fixed bugs in the Abort Experiment button.
%   Fixed bugs in the deleting of initialization gui
% 2016.11.29
%   Changed handling stimTrain.
%       Now stimTrain is a structure, with a field for each output channel.
%       This allows each output to be a vector or a matrix. 
%       If matrices are present, playrecARLas.queue loops over each column.
%       This allows random noise to be played in one or more channels,
%       among other things. 
%       THERE IS NO BACKWARDS COMPATIBILITY FOR THIS UPDATE: All experiment
%       files must be changed so that they load stimTrain using the channel
%       fields. Example of old version: stimTrain = randn(44100,2);
%               Example of new version: stimTrain.Ch1 = randn(44100,1);
%                                       stimTrain.Ch2 = randn(44100,1);
% 2016.12.07
%   Added userInfo as a public property of playrecARLas.
%       This will be saved as part of the header file.
%       This allows users to save important information from within
%       experiment files. 
%               Example: f2 = 1000; % frequency in Hz
%                        obj.objPlayrec.userInfo.f2Freq = 1000; 
%                        obj.objPlayrec.run % run the stimulus                         
%   Added backwards compatibility to stimTrain handling. The old way of
%       using stimTrain (simply make a matrix, with each column being an
%       output channel can still be used. A new property, stimTrainEx
%       (Ex stands for "explicit") can also be used to specify channels
%       explicitly, and to use either vectors or matrices for each channel.
%               Example of old version: stimTrain = randn(44100,2);
%               Example of new version: stimTrainEx.Ch1 = randn(44100,1);
%                                       stimTrainEx.Ch2 = randn(44100,1);
% 2016.12.13
%   Added automatic creation of the Data folder on startup, if it doesn't
%       exist. This will enable versions saved on Github to not include an 
%       empty data folder.
%
% 2017.04.11
%   Changed a line in playrecARLas.queue. Now only nReps + 1 copies of the
%       are added to the queue. (Used to be 2*nReps). Since we no longer do
%       online artifact rejection, the extra overhead is not necessary and
%       can cause unexpected performance in certain conditions.
%
% 2017.04.19
%   Fixed a bug in arlas.m line 814.
%       txt did not exist; added additional line: txt = 'problem with specified channels';
%
% 2017.06.14
%   Fixed a bug in playrecARLas.m line 770.
%       ARlas was failing when users tried to switch output screen views at
%       the very end of a run. This fix simply returns and does not provide
%       the update.
%
% 2017.07.14
%   Major overhaul and partial re-write of playrecARLas.m.
%       Changes include removal of holding pen, addition of zero-phase IIR
%       filtering on the fly, ability to record very short stimulus trains
%       (shorter than the system delay). This fixes a timing glitch,
%       wherein short sequences sometimes showed an apparent sudden shift
%       in delay. Collected data data are now saved when abort button is pressed. 
%       Fixed numbering tag system so that if a value in the sequence is
%       missing, it will not be filled in; rather, the program now looks
%       for the largest value in the sequence and increments by 1.
%       arlas.m now looks to make sure that the other two classes are of
%       matching version numbers, otherwise alerts user to discrepancy.
%
% 2017.07.20
%   Fixed playrecARLas to check for current Matlab version number and use
%       the function strfind instead of contains if < 2016b is being used.
%       in playrecARLas, now suppress all error messages except the first
%       one. In playrecARLas, fix an error in retrieve data that was
%       counting sample delay from delay instead of delay + 1 sample. In
%       playrecARLas, updated to use playrec's system delay estimate if
%       sufficienty close to user's measured delay.
%
% 2017.07.24
%   Fixed arlas to allow for single repetion (nReps = 1). In this case, no
%       correction is applied for system delay. Numerous other small
%       changes made to improve stability. This version suppresses multiple 
%       error messages, hopefully enabling easier debugging.
%
% 2017.07.26
%   Fixed a plotting glitch that was overwriteing axis labels on other open
%       figures. This was occurring in playrecARLas when new plots were
%       made for the waveforms. There is now a visible, short, delay;
%       however, this is preferable to overwritting other figure axis
%       labels.
%
% 2017.08.11
%   Added a helpful error message to users regarding the old (no longer
%       supported) syntax for setting nReps. New message prints the fix on the
%       screen to help users update their code.
%
% 2017.08.14
%   Turned off one of the features implemented on 2017.07.20. No longer
%       uses playrec's system delay estimate. Found it was not providing
%       benefit on windows systems, and was causing problems on Mac.
%
% 2017.11.04
%   playrecARlas now returns the IIR filter information (SOS and G). This
%       allows users to apply the filter post-hoc, or to appy the same
%       filter to stimuli that were applied to the recordings. (This was
%       used to improve the ThevCal fitting.)
%
% 2018.11.14
%   playrecARlas updated to work with sound cards that have a larger number
%       of channels (or virtual channels). Updates around line 595.
%
% 2021.07.07
%    THIS VERSION HAS added hack to switch data writing to D drive, NOT C
%    DRIVE!
%
% 2021.31.07
%    Fixed the hack above more elegantly. Implemented two new files in 
%    Peripheral\sysConfigs:  pathSetup.m and fileNamesSetup.m. The former
%    is used to set things like writing to other locations. The latter is
%    to suggest experiment names and pre-populate the ID field accordingly.
%    This will help maintain consistency across different experimenters.
%    Also made a major change to file structure writing. Now everyting goes
%    as follows:
%       Data\
%             SubjID\
%                    SubjID_31JUL2021\
%                                     SubjID_31JUL2021_testName_run1
%
%    This will be a big help keeping track of multiplie runs withing and
%    across multipl visits from the same subjects. Several core updates
%    were made ot facilitate these changes.
%
% 2021.04.08
%    Had a big problem forgetting to disable ER10X output limiter.
%    Implemented a reminder dialog box when arlas starts. ARLas checks
%    hardwareSetup (with a new field) to see whether limiter is supposed to
%    be being used.
%
% 2022.03.23
%    Now using standardized file structure for calibrations.
%    Added code to arlas so that if thevCal, LTCal, micCal, or calCheck
%    folders are not present under ...Peripheral\calibrations\, they are
%    automatically created.
% 2022.04.01
%    Updated the function enterID in arlas.m. Now adds newly created
%    folders to the path.
%
% 2024.10.04
%   In method "cleanUp" (in playrecARLas), replaced function call fopen with
%   function call openedFiles. The code now tries both functions to retain 
%   backwards compatability.
%
% 2025.01.07
%   In playrecARLas, added new public variable "servedNeat". This works
%   with the audiogram and plays a stripped down version, with no attempt
%   to zero pad. Also does not re-initialize and does not write data. Much
%   faster. Set =0 (default, old standard way) or =1 for this new quick
%   version.

