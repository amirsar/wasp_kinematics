%09/03/2022 Main script for forward flight in E. mundus article
% call the following funcions: (In order of appearance) 
% - 'kine_cycles_frames_validation': validate that chosen frames contain data & extract prior and posterior frames
% - 'extract_coords': extract coordinate of chosen frames from kine solution
% - 'kine_wings': calculate Euler angles from wings coordinates
% - 'mean_V_from_kine_solution': calculate 'Advance ration', 'Frq.', 'Vxyz', 'Vxy', 'Vz','V_Ltip','V_Rtip'
% - 'wing_trajectory_data': Calculate & plot 3D wing trajectory with body possition and chord. Make use of func.: 'wing_trajectory_plot'
% - 'Meanangles': plot Euler angles of flapping cycles with mean
% - 'CFDbody_measurments': Store cycle's data that relevant for CFD
% - 'Vxy_trend_on_angles': Plot wings angles of all input cycles, while color coded by Vxy
% - 'amplitude_distance_func': Calculate distance between dorsal & ventral stroke reversals in flapping kineatics data
% - 'Vxy_trend_on_wingtip_relative2body': Calculate wingtip distance (normalized to body length). Plot: 1)cycle's wingtip trajectory relative wing-base, one over the other, colored by cycle's body's Vxy. 2)max, median & min V cycles individually, colored by to wingtip instentanious V. 3)stroke plane amlitude. 4)Wingtip 3D speed relative to body  5)Body 3D speed
% - 'Vxy_trend_on_wingtip_relative2V': Calculate pitch. Plot: 1)cycle's wingtip trajectory relative to body's velocity vector, one over the other, colored by cycle's body's Vxy. 2)wingtip 3D speed relative to air. 3)mid-cord 3D velocity relative to air at body's velocity vector direction. Make use of func.: 'wingtip_world_viewpoint' which Load & store insantious force data & Plot each cycle's kinematics with force vectors overlay
% - 'Vxy_trend_on_AoA_n_pitch': Calculate Angle-of-Attack & force components. Plot body pitch & angle of attack angles, cycles color coded by Vxy
% - 'AoA_phases': Plot comparison between each cycle's upstroke to downstorke AoA & choed velcity
% - 'wing_length_related_data': Calculate wingtip's traveled distance normalized to body length
clear
close all
%% load functions
addpath('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Code\DigitizingTools\kine analysis');
addpath('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Code'); 
%% extracting data
pth = uigetdir('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Digitizing\Kinematics\WaspKinematics\Kine\fitted','Select The Directory That Contains Solution Directories From Kine'); %ask user for path of the data directories
fprintf('Extract data from: %s\n',pth);
cd(pth); %change active directory to the folder that contains the kine solution files
[file_name,file_path]=uigetfile('*.xlsx','Select ''wing cycles data file''','D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\wing_cycles_data_file'); %ask user to choose file that contain wing cycle names and frames
[frames_index,trial_names]=xlsread(fullfile(file_path,file_name)); %read data. 1st variable is strings from the excel file, the 2nd is the numeric data
draw_plot = strcmp(str2num(questdlg(sprintf('Draw each cycle''s plots?\n (if uncertain - choose the default)'),'Draw separate plots?','Yes','No','No')),'Yes'); %ask user whether to Draw plots
use_filter=questdlg(sprintf('Filter original digitizing data?\n (if uncertain - choose the default)'),'Filter digitization data?','Gaussian','Moving avg.','No','No'); %ask user 
if strcmp(use_filter,'No')
    use_filter=false;
end
repeats=size(frames_index,1); %determain how many wing cycles according to data index file
tcorr = (0.02:0.02:1.0)'; %50 X values between 0 to 1 (breakpoints) that need to be interpolant data (Y value) in the wing cycle
Hz=6000; %recording frame rate
AllData_V_tip_storage=table; CFDbody4sim=[]; %storage variables
for wing_cycle=1:repeats %repeat for each wing_cycle. 
    %% reset variables
    idx=1; %will be used in the following functions to indicate there is only one wing cycle
    continuous='f'; %will be used to indicate following functions that there isn't sequance of wing cycles
    V_data = table; V_tip_storage.body_center.x=[]; pitch=[]; coords_storage=[]; %will be used for 'mean_V_from_kine_solution' function
    rotated_coords=struct(); coords_index=[0,0]; stroke_plane = table;%will be used for 'wing_tip_trajectory' function
    %% calc
    trial=char(trial_names(wing_cycle+1,1)); %recieve solution name
    fprintf('Trial: %s, cycle %d. Switch wings: %s\n',trial,frames_index(wing_cycle,2),char(trial_names(wing_cycle+1,6)))
    frames = [frames_index(wing_cycle,3),frames_index(wing_cycle,4)]; %recieve wing_cycle's frame
    cd(trial) %change active directory to the current trial
    s=what; %recieve *.mat files names in current folder
    load(char(s.mat)); %load trial's solution data
    [frames, prior_frames, posterior_frames] = kine_cycles_frames_validation(data, frames, 1, 1, continuous); %testing in order to make sure there are digitizing data in the cycle, and marking extra data in edges
    if strcmp(trial_names(wing_cycle+1,6),'Yes') %recieve correction type
        switch_wings=2;
    elseif strcmp(trial_names(wing_cycle+1,6),'No')
        switch_wings=1;
    end
    [data, coords] = extract_coords(data,frames(idx),frames(idx+1),prior_frames,posterior_frames,switch_wings,Hz,tcorr,use_filter); %extract coordinates
    angles=kine_wings(data, frames(1:2)); %extract angles
    angles = radtodeg(angles'); %convert from column to row, from radians to degrees
    angles([1,4],:)=angles([1,4],:).*(-1); %reverse direction so that positive will point headwise
    elevation=angles([2,5],:); %exctract elevation angles
    horizontal=mean(angles([1,4],:),1); %exctract flapping/horizontal angles
    incidence=mean(angles([3,6],:),1); %exctract incidence angles
    [V_data,V_tip_storage,pitch,coords_storage]=mean_V_from_kine_solution(coords, data, continuous,frames(idx),frames(idx+1),Hz,V_data,idx,idx,draw_plot,tcorr,V_tip_storage, prior_frames, posterior_frames,elevation,pitch,coords_storage); %calc mean trial velocity
    [rotated_coords,coords_index,stroke_plane]=wing_trajectory_data(coords,idx,idx,rotated_coords,coords_index,draw_plot,tcorr,V_tip_storage,stroke_plane,V_data,degtorad(elevation)); %plot wing tip trajectory
    horizontal=mean(angles([1,4],:),1); %mean flapping/horizontal angles
    [horizontal_minima,horizontal_minima_loc] = findpeaks(mean(horizontal.*(-1),1),'SortStr','descend','MinPeakDistance',round(length(horizontal)/4)); %find local negative peaks while minimal distance between peaks is at least 1/4 wing cycle. returns the peaks in descending order, from the largest to the smallest value.
    horizontal_minima=(-1).*horizontal_minima; %return values to original direction
    cycles_names(wing_cycle).name=trial; %for use in 'Meanangles' & 'Vxy_trend_on_angles' functions
    angles_storage(2*wing_cycle-1).Vxy=V_data{1,6};
    angles_storage(2*wing_cycle-1).advance_ratio=V_data{1,3};
    angles_storage(2*wing_cycle-1).pitch=pitch;
    angles_storage(2*wing_cycle-1).elevation=degtorad(angles(2,:)); angles_storage(2*wing_cycle).elevation=degtorad(angles(5,:)); %for use in the following function
    angles_storage(2*wing_cycle-1).horizontal=degtorad(angles(1,:)); angles_storage(2*wing_cycle).horizontal=degtorad(angles(4,:)); %for use in the following function
    angles_storage(2*wing_cycle-1).incidence=degtorad(angles(3,:)); angles_storage(2*wing_cycle).incidence=degtorad(angles(6,:)); %for use in the following function
    rotated_coords_storage(wing_cycle).data=rotated_coords;
    coords_storage.cycle_name=str2double(sprintf('%d.%d',frames_index(wing_cycle,1), frames_index(wing_cycle,2))); 
    AllData_coords_storage(wing_cycle)=coords_storage;
    V_tip_storage.cycle_name=str2double(sprintf('%d.%d',frames_index(wing_cycle,1), frames_index(wing_cycle,2))); 
    AllData_V_tip_storage_rawData(wing_cycle)=V_tip_storage;
    AllData_V_tip_storage=[AllData_V_tip_storage; V_data]; %add 'V_data' under existing data
    report(wing_cycle,:)=[frames_index(wing_cycle,1), str2double(sprintf('%d.%d',frames_index(wing_cycle,1), frames_index(wing_cycle,2))),angles_storage(2*wing_cycle-1).pitch.ventral_face_down,mean([V_data{1,10},V_data{1,11}]), mean([V_data{1,12},V_data{1,13}]), V_data{1,4}, V_data{1,6}, V_data{1,9}, min(incidence), max(incidence), min(horizontal), max(horizontal), min(mean(elevation,1)), max(mean(elevation(:,round(size(elevation,2)/2):end),1)), horizontal_minima(2), min(pitch.angle), max(pitch.angle), mean(pitch.angle), pitch.stroke_plane_relative2body.abs_flapping_angle_from_horizontal_direction, pitch.stroke_plane_relative2body.flapping_angle_from_horizon, pitch.stroke_plane_relative2body.flapping_angle_from_body, incidence(round(0.7*length(incidence)))]; %save data for report: trial, cycle, body ventral direction, mean R&L wingtip groundV3D, mean R&L wingtip arialV3D, Frq., Vxy, Vz, min(incidence), max(incidence), min(horizontal), max(horizontal), min(elevation), tmp_max(elevation), tmp_min_horizontal, min pitch, max pitch, average pitch, SP from horizontal direction, SP from horizon, SP from body, incidence at 0.7t 
    angles_storage(2*wing_cycle-1).cycle_name=report(wing_cycle,2);
    angles_storage(2*wing_cycle-1).horizontal_tmp_min_relative_loc=(horizontal_minima_loc(2)/length(horizontal))*100; %end of fling relative location during the flapping cycle, according to horizontal angle
%     [clap_length,downstroke_length,upstroke_length] = cycle_phases_func(angles(2,:),angles(5,:),incidence); %calc cycle_phases length
%     cycle_phases_length(wing_cycle,:)=[report(wing_cycle,2),report(wing_cycle,7),clap_length,downstroke_length,upstroke_length]; %store it
    cd(pth); %change active directory to the folder that contains the kine solution files
    %% plot of each cycle individualy
    if draw_plot
%         wingtip_world_viewpoint(wing_cycle,V_data{1,6},AllData_coords_storage);
        mean_cycle = Meanangles(idx*2,cycles_names(wing_cycle),tcorr,angles_storage((2*wing_cycle-1):(2*wing_cycle))); %plot angles
        for i=4:7
            if i==4
                figure(i-1)
            else
                figure(i)
            end
            s=sprintf('%s.%d_%d.jpg',trial(1:13),frames_index(wing_cycle,2),i);
            saveas(gcf,s)
        end
    end
    CFDbody4sim=CFDbody_measurments(wing_cycle,AllData_V_tip_storage_rawData,V_data,pitch,coords,CFDbody4sim);
    clearvars -except CFDbody4sim cycle_phases_length report rotated_coords_storage AllData_coords_storage AllData_V_tip_storage_rawData AllData_V_tip_storage angles_storage V_data cycles_names wing_cycle Hz tcorr repeats use_filter draw_plot trial_names frames_index pth %clear all variables except those needed to run the loop and those who store data
    close all
end
%% plots for article
[Flp, Elv, Inc, CFDwings_angles4sim]=Vxy_trend_on_angles(report(:,7),angles_storage,cycles_names,frames_index(:,2),tcorr);
[amplitude_distance,mean_amplitude_distance]=amplitude_distance_func(AllData_coords_storage,cycles_names);
[I,insect_corners_x,insect_corners_y,wingtip_travel_dist_body]=Vxy_trend_on_wingtip_relative2body(report(:,7),report(:,6),rotated_coords_storage,AllData_V_tip_storage_rawData,amplitude_distance);
[wingtip_travel_dist_air,Vtip_ground_3d,angles_storage,Force]=Vxy_trend_on_wingtip_relative2V(report(:,7),AllData_coords_storage,angles_storage,AllData_V_tip_storage_rawData);
[abs_AoA,AoA,force_components,angles_storage]=Vxy_trend_on_AoA_n_pitch(report(:,7),angles_storage,tcorr,AllData_coords_storage,Vtip_ground_3d);
for trial=1:size(AoA,1) %for each trial
    [AoA_tmp_max(trial),AoA_tmp_max_index(trial)]=max(AoA(trial,1:12)); %extract max AoA between 0 to 0.24 of the flapping cycle
end
report=array2table([report(:,:),AoA_tmp_max',mean(AoA,2), wingtip_travel_dist_air', wingtip_travel_dist_body',amplitude_distance{:,3}],'VariableNames',{'trial', 'cycle', 'body_ventral_direction', 'wingtip_groundV3D', 'wingtip_arialV3D', 'Frq', 'Vxy', 'Vz', 'min_incidence', 'max_incidence', 'min_horizontal', 'max_horizontal', 'min_elevation', 'tmp_max_elevation', 'tmp_min_horizontal', 'min_pitch', 'max_pitch', 'average_pitch', 'SP_from_horizontal_direction', 'SP_from_horizon', 'SP_from_body', 'incidence_at_7of10t','AoA_tmp_max','mean_AoA', 'wingtip_travel_dist_relative_to_air', 'wingtip_travel_dist_relative_to_body','amplitude_distance'});
% cycle_phases_length=table(cycle_phases_length(:,1),cycle_phases_length(:,2),cycle_phases_length(:,3),cycle_phases_length(:,4),cycle_phases_length(:,5),'VariableNames',{'cycle','Vxy','clap_length','downstroke_length','upstroke_length'});
AoA_phases;
normalized_wingtip_distance=wing_length_related_data(AllData_coords_storage);
wing_drawing %plot supplementary figure of wing components. Independend script
pitch2Fdirection_figure %plot pitch to force-direction figure. Independend script