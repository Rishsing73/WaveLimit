%% this file will read the data from plexon and behavior and pass it throught errasr and then save it to nex for cleaning/sorting
% read all the relevant files

disp('please select the .pl2 file');
filename = internalPL2ResolveFilename(''); % load neural file
pl2 = PL2GetFileIndex(filename);

disp('please select the .mat file corrosponding to the neural recording with correct timestamps');
[file,path] = uigetfile; % load beahvioral file
correctTS = load(fullfile(path,file));

disp('please select the .mat file corrosponding to the neural recording with behavior data');
[fileB,pathB] = uigetfile; % load beahvioral file
bhv = load(fullfile(pathB,fileB));

%%
% get the data of plexon in a matrix and then pick out the stim data for
% errasr
channel_source = 'WB'; % parameter for plexon file
nchannels_per_shank = 16;% channels recorded from

[ad] = PL2AdBySource(filename, channel_source, 1);
dat_plx = NaN(nchannels_per_shank,length(ad.Values)); % initialize an empty matrix
dat_plx(1,:) = ad.Values;

for i = 2:nchannels_per_shank % Number of event channels
    [ad] = PL2AdBySource(filename, channel_source, i);
    dat_plx(i,:) = ad.Values;
end

% dat_plx is the RECORDING

%% Pick the correct time stamp of stimulation and pick the data for errasr

ct = correctTS.data.('sham_stim')(1); % timestamp file of stimulation with correct time
trial_start = correctTS.data.('Trial')(1); % picking the relevant trial

bt = bhv.trials(trial_start:(trial_start+height(correctTS.data)-1)); % select the range

shamstim = [bt.('stimOnTime')];
td = shamstim(1) - ct;

stim_ts = bt{:, timestamp_list('stimOnTime')} - td; % time stamps corrected to plx file
stim_ts = stim_ts(bt{:, timestamp_list('stimTrial')}==1); % stimulation trials
stim_ts = stim_ts.*40000; % converting time into frequency timestamp

pre_window = 10000; % the window to take out the before the stim
post_window = 29999; % the window after stimOnTime


%% grab few trials
stim_ts = stim_ts(1:4);
data_raw = NaN(length(stim_ts),40000,16);

for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    data_raw(i,:,:) = spike_raw(:,takeout:takeafter)';
end
%% Grab spontaneous data as well - If using Becket_pca

% data = NaN(40000,16,length(stim_ts));
% for i = 1:numel(stim_ts)
%     takeout = floor(stim_ts(i)- 50000);
%     takeafter = floor(stim_ts(i) - 10000-1);
%     data(:,:,i) = spike_raw(:,takeout:takeafter)';
% end

% % Plotting parameters`
yOffset = 6;

figure;
for channel = 1:16
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ data_raw(:,channel,1)); %plot to verify that every thing is right
end

xlabel('time');
ylabel('channels');


%% ERRASR params

opts = ERAASR.Parameters();
opts.Fs = 40000; % samples per second
Fms = opts.Fs / 1000; % multiply to convert ms to samples

opts.thresholdHPCornerHz = 250;
opts.thresholdChannel = 8;
opts.thresholdValue = 3;

opts.alignChannel = 1;
opts.alignUpsampleBy = 12;
opts.alignWindowPre = Fms * 0.5;
opts.alignWindowDuration = Fms * 15;

% 60 ms stim, align using 20 ms pre start to 110 post
opts.extractWindowPre = Fms * 30;
opts.extractWindowDuration = Fms * 560;
opts.cleanStartSamplesPreThreshold = Fms * 1;
        
opts.cleanHPCornerHz = 10; % light high pass filtering at the start of cleaning
opts.cleanHPOrder = 4; % high pass filter order 
opts.cleanUpsampleBy = 1; % upsample by this ratio during cleaning
opts.samplesPerPulse = Fms * 5; % 3 ms pulses - changes to 5
opts.nPulses = 100;

opts.nPC_channels = 12;
opts.nPC_trials = 3;
opts.nPC_pulses = 6;

opts.omit_bandwidth_channels = 3;
opts.omit_bandwidth_trials = 1;
opts.omit_bandwidth_pulses = 1;

opts.alignPulsesOverTrain = true; % do a secondary alignment within each train, in case you think there is pulse to pulse jitter. Works best with upsampling
opts.pcaOnlyOmitted = true; % if true, build PCs only from non-omitted channels/trials/pulses. if false, build PCs from all but set coefficients to zero post hoc
opts.lamda = 0.1;
opts.cleanOverChannelsIndividualTrials = true;
opts.cleanOverPulsesIndividualChannels = false;
opts.cleanOverTrialsIndividualChannels = false;

opts.cleanPostStim = false; % clean the post stim window using a single PCR over channels

opts.showFigures = false; % useful for debugging and seeing well how the cleaning works
opts.plotTrials = 1; % which trials to plot in figures, can be vector
opts.plotPulses = 1; % which pulses to plot in figures, can be vector
opts.figurePath = [pwd '/timon_171101_postbug'];% folder to save the figures
% if ~exist(opts.figurePath, 'dir')
%     mkdir(opts.figurePath)
% else
%     disp('The folder does exist.');
% end

opts.saveFigures = false; % whether to save the figures
opts.saveFigureCommand = @(filepath) print('-dpng', '-r300', [filepath '.png']); % specify a custom command to save the figure

%% Do alignment and cleaning procedure

[dataCleaned, extract] = ERAASR.cleanTrials(data_n,data, opts);

%% Create the plot iteratively - to verify that dta has been cleaned
figure(1);
hold on;
yOffset = 5;
for channel = 1:16
    yValues = channel + yOffset;  % Add offset to y-values
    plot(yValues+ dataCleaned(1,:,channel));
end

%% Put the values back in dat_plx

for i = 1:numel(stim_ts)
    takeout = floor(stim_ts(i)-pre_window);
    takeafter = floor(stim_ts(i) + post_window);
    dat_plx(:,takeout:takeafter) = squeeze(dataCleaned(i,:,:))' ;
end

%% write the data into nexfile or should we combine the data

