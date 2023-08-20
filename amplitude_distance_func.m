%===========================================================================================%
% 13/07/2022 - By Amir Sarig
% - Distance between dorsal & ventral stroke reversals in flapping kineatics data
% - Called from 'extract_data_wasp_forward_flight.m'
%===========================================================================================%
function [amplitude,mean_amplitude_distance]=amplitude_distance_func(AllData_coords_storage,cycles_names)
% clear
% The following variable was manually created. It contain trial name &
% Vxy, with categorization (by frame and correlated percentage)of each cycle to 3 phases. base on manually
% observation at the records
s=sprintf('File named ''phases_duration.mat'' is missing.\n It was created manually.\n It is a table containing a row for each wing cycle, The columns are as follow:\n cycle_name, Vxy, downstroke_frames, upstroke_frames, clap_phase _frames, downstroke_prct, upstroke_prct, clap_phase _prct'); %error message
if strcmp(cycles_names(1).name(1:4),'wasp') %if working with wasps' data which I manually prepared phases seperation
    if exist('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\phases_duration.mat', 'file') == 2 %if the file exist
        load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\phases_duration.mat') %load phases 
        % load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\coordinates.mat') %load coordinates
        for trial=1:size(cycle_phases,1) %for each flapping cycle
            ventral_reversal_frame=cycle_phases{trial,3}; %extract frame of end of downstroke
        %     dorsal_reversal_frame=cycle_phases{trial,3}+cycle_phases{trial,4}+1; %begging of clap
            dorsal_reversal_frame=cycle_phases{trial,3}+cycle_phases{trial,4}+1:cycle_phases{trial,3}+cycle_phases{trial,4}+cycle_phases{trial,5}; %extract frames of all the clap phase
            ventral_reversal_coords=(AllData_coords_storage(trial).original(4).data(ventral_reversal_frame,:)+AllData_coords_storage(trial).original(6).data(ventral_reversal_frame,:))./2; %mean of right & left wingtips location at ventral_reversal point
        %     dorsal_reversal_coords=(AllData_coords_storage(trial).original(4).data(dorsal_reversal_frame,:)+AllData_coords_storage(trial).original(6).data(dorsal_reversal_frame,:))./2; %mean of right & left wingtips
            dorsal_reversal_coords=(mean(AllData_coords_storage(trial).original(4).data(dorsal_reversal_frame,:),1)+mean(AllData_coords_storage(trial).original(6).data(dorsal_reversal_frame,:),1))./2; %mean location of wingtips during clap & mean of right-left  
            amplitude(trial,:)=[cycle_phases{trial,1},cycle_phases{trial,2},norm(dorsal_reversal_coords-ventral_reversal_coords),cycle_phases{trial,6},cycle_phases{trial,8}]; %distance between ventral reversal point to mean location during clap
        end
        amplitude=array2table(amplitude,'VariableNames',{'cycle_name','Vxy','amplitude','ventral_stroke_phase','dorsal_stroke_begging_phase'}); %store data
        cycles_affiliation=floor(amplitude{:,1}); %extract insect number that performed each flapping cycle
        individuals_names=unique(cycles_affiliation,'stable'); %save each insect number once. keep the order as appear in the data
        for individual=1:length(individuals_names) %for each individual insect
            current_individual=individuals_names(individual)==cycles_affiliation; %exctract the rows with this individual data
            temp(individual)=mean(amplitude{current_individual,3}); %calc mean amplitude for each individual 
        end
        mean_amplitude_distance=temp'; %rotated mean amplitude array from row to coulumn
    else %if data don't exist
        errordlg(s,'''phases_duration.mat'' file doesn''t exist');
    end
else %if data don't exist
    fprintf('%s\n',s);
    amplitude=nan(size(AllData_coords_storage,2),3);
    mean_amplitude_distance=nan(size(AllData_coords_storage,2),1);
end
end