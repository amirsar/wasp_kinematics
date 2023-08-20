% plot wingtip trajectory leveled to the body's velocity direction. Sample of 2 following flapping cycles from 3 trials with different Vxy.
clear
close all
%% load functions
addpath('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Code\DigitizingTools\kine analysis');
addpath('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Code'); 
%% extracting data
pth = uigetdir('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Digitizing\Kinematics\WaspKinematics\Kine\fitted','Select The Directory That Contains Solution Directories From Kine'); %ask user for path of the data directories
fprintf('Extract data from: %s\n',pth);
cd(pth); %change active directory to the folder that contains the kine solution files
[file_name,file_path]=uigetfile('*.xlsx','Select wing cycles file','D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article'); %ask user to choose file that contain wing cycle names and frames
[frames_index,trial_names]=xlsread(fullfile(file_path,file_name)); %read control table data. 1st variable is strings from the excel file, the 2nd is the numeric data
% leveled=figure; hold on
% centered=figure; hold on
%% manually enter data section
frames(1,1:3)=[frames_index(1,3),frames_index(2,4),2]; %7.14-7.15 - slow
frames(2,1:3)=[frames_index(5,3),frames_index(4,4),1]; %6.2-6.3 - medium
frames(3,1:3)=[frames_index(7,3),frames_index(6,4),1]; %8.1-8.2 - fast
% frames(4,1:3)=[frames_index(10,3),frames_index(11,4),2]; %11.2-11.3 - fastest, dorasl body side facing ground
trials_folders=[trial_names(1+1,1);trial_names(1+4,1);trial_names(1+6,1);trial_names(1+10,1)];
%% 
for trial=1:size(frames,1)
    second_clap_frame(trial)=round((frames(trial,2)-frames(trial,1))/2);
    switch_wings(trial)=frames(trial,3); %swith wings status for the above trials
end
relative_Vxy=0:1/(trial-1):1; %will be use to color trajectory in plot
tcorr = (0.02:0.02:1.0)'; %50 X values between 0 to 1 (breakpoints) that need to be interpolant data (Y value) in the wing cycle
Hz=6000; %recording frame rate
continuous='f'; %will be used to indicate following functions that there isn't sequance of wing cycles
for trial=1:size(frames,1) %repeat for min,median,max Vxy
    color=[relative_Vxy(trial) 0 1-relative_Vxy(trial)]; %trials uniqe color
    cd(char(trials_folders(trial))) %change active directory to the current trial
    fprintf('Trial: %s\n',char(trials_folders(trial)));
    s=what; %recieve *.mat files names in current folder
    load(char(s.mat)); %load trial's solution data
    [trial_frames, prior_frames, posterior_frames] = kine_cycles_frames_validation(data, frames(trial,1:2), 1, 1, continuous); %testing in order to make sure there are digitizing data in the cycle, and marking extra data in edges
    coords = extract_coords(data,trial_frames(1),trial_frames(2),prior_frames,posterior_frames,switch_wings(trial),Hz,tcorr); %extract coordinates
    for i=1:size(coords,2)
        trial_coords(i).data=coords(i).original';
%         trial_coords(i).data=coords(i).interpolation;
    end
%% calc normalize wingtip trajectory
    body_center=(trial_coords(1).data(:,1)+trial_coords(2).data(:,1))./2; %calc body center in the 1st frame
    for digitizing_point=1:size(trial_coords,2) % start each cycle at (0,0,0), repeat for each digiting point. order: 1)head;2)abd;3)Lbase;4)Ltip;5)Rbase;6)Rtip;7)Llead;8)Rlead;9)Ltrail;10)Rtrail;
        normalize(digitizing_point).coords=(trial_coords(digitizing_point).data-repmat(body_center,1,size(trial_coords(digitizing_point).data,2)))'; %(all frames of each digitizing_point minus body center at 1st frame) 
    end
    % plot data
    figure
    tip=(normalize(4).coords'+normalize(6).coords')./2;
    if tip(1,1)>tip(1,end) %if the last frame is to the left on the X axis from the 1st frame
        tip(1,:)=-tip(1,:); %horizontaly flip the coordinates
    end
    temp=plot(tip(1,:),tip(3,:),'Color',color,'LineWidth',2); %plot 2D trajectory. color it according to Vxy
%     lgnd2(trial).data=temp; %store for legend
    xlim([-0.001 0.003])
    axis equal
    set(gca,'FontSize',20);
    xlabel('Horizontal axis (m)', 'FontSize',24)%Body velocity axis
    ylabel('Vertical axis (m)', 'FontSize',24)%Vertical perpendicular axis
%% calc rotated trajectory
    body_center=(normalize(1).coords(:,:)+normalize(2).coords(:,:))./2; %calc body center of each frame: (head+abd)/2
    velocity_vector=body_center(end,:)-body_center(1,:); %cycle's velocity_vector is BC's last frame minus 1st frame
    horizontal_angle=2*pi-atan2(velocity_vector(2),velocity_vector(1)); %calculate cycle's velocity_vector projection on XY plain
    z_rotated=([cos(horizontal_angle) -sin(horizontal_angle) 0; sin(horizontal_angle) cos(horizontal_angle) 0; 0 0 1]*velocity_vector')';% rotated around vertical axis in order to have unified horizontal direction on the XZ plane. So that X is the horizontal direction of body velocity.
    vertical_angle=2*pi-acos(dot(z_rotated([1 3]),[1,0])./(norm(z_rotated([1 3]))*norm([1,0]))); %calculate cycle's velocity_vector Vertical_angle from it's horizontal direction
    fprintf('Rotation - horizontal angle: %1.3f degrees, vertical angle: %1.3f degrees\n',rad2deg(horizontal_angle),radtodeg(vertical_angle))
    for dim=1:3 %calc mean of left&right wing. repeat for each axis
        mean_L_R(1).coords(dim,:)=mean([normalize(4).coords(:,dim), normalize(6).coords(:,dim)],2); %wingtip mean is calculated along columns
        mean_L_R(2).coords(dim,:)=mean([normalize(7).coords(:,dim), normalize(8).coords(:,dim)],2); %leading edge mean is calculated along columns
        mean_L_R(3).coords(dim,:)=mean([normalize(9).coords(:,dim), normalize(10).coords(:,dim)],2); %trailing edge mean is calculated along columns
        mean_L_R(4).coords(dim,:)=normalize(1).coords(:,dim); %head
        mean_L_R(5).coords(dim,:)=normalize(2).coords(:,dim); %abdomen
    end
    for object=1:size(mean_L_R,2) %repeat for each object. order: 1)wingtip 2)lead 3)trail 4)head 5)abdomen 6)V_mid_chord
        for frame=1:size(mean_L_R(object).coords,2) %rotate each frame's wingtip coordinate by horizontal & vertical angles
            z_rotated=([cos(horizontal_angle) -sin(horizontal_angle) 0; sin(horizontal_angle) cos(horizontal_angle) 0; 0 0 1]*mean_L_R(object).coords(:,frame))';% rotated in order to oriente XZ plain
            y_rotated=([cos(-vertical_angle) 0 sin(-vertical_angle); 0 1 0; -sin(-vertical_angle) 0 cos(-vertical_angle)]*z_rotated')';% rotated in order to oriente XY plain
            object_relative2V(object).coords(:,frame)=y_rotated'; %store in a way that each axis is a row
        end
    end
    body_center=((object_relative2V(4).coords+object_relative2V(5).coords)./2); % rotated body center
    body_vector=(object_relative2V(4).coords-body_center)-(object_relative2V(5).coords-body_center); % normalized rotated head-abdomen
    wingtip=object_relative2V(1).coords-body_center; % normalized rotated wingtip
    clap_frame=round(0.95*size(object_relative2V(1).coords,2)); %middle of clap should be around 0.9 flapping cycle length, so 0.95 is for the 0.9 of the 2nd cycle 
    transverse=cross(wingtip(:,clap_frame),body_vector(:,clap_frame)); %body transverse axis, from right to left of the body: cross(wings length,body_vector) 
    dorsal_norm=cross(body_vector(:,clap_frame),transverse); %dorsal perpendicular to the plane of the body
    roll_angle=2*pi-acos(dot(dorsal_norm([1 3]),[0,1])./(norm(dorsal_norm([1 3]))*norm([0,1]))); %calculate cycle's roll_angle between Z axis to it's rotated body dorsal direction at 2nd clap
    for object=1:size(mean_L_R,2) %repeat for each object. order: 1)wingtip 2)lead 3)trail ...
        for frame=1:size(object_relative2V(object).coords,2) %rotate each frame's wingtip coordinate
            shifted=object_relative2V(object).coords(:,frame)-body_center(:,frame); %set body center as (0,0) because rotation is around the beginning of the axes
            x_rotated=([1 0 0; 0 cos(roll_angle) -sin(roll_angle); 0 sin(roll_angle) cos(roll_angle)]*shifted)';% rotated by roll_angle
            object_relative2V(object).coords(:,frame)=x_rotated'+body_center(:,frame); %store in a way that each axis is a row, return to original position
        end
    end
    if (360-rad2deg(roll_angle))<=90
        line_style='-'; %plot continous curve if ventral_face_down
    else
        line_style='--'; %plot dotted curve if ventral_face_up
    end
    for dim=1:3 %calc distance between wingtip to middle of chord. repeat for each axis
        mean_L_R_chord_shift(dim,:)=object_relative2V(1).coords(dim,:)-mean([object_relative2V(2).coords(dim,:); object_relative2V(3).coords(dim,:)],1); %distance between middle of the cord to wingtip. mean is calculated along columns
    end
%% plot data
    figure
    hold on
    
%     temp=plot(object_relative2V(1).coords(1,:),object_relative2V(1).coords(3,:),'Color',color,'LineWidth',2,'LineStyle',line_style); %plot 2D trajectory. color it according to Vxy
    
    for i=1:numel(object_relative2V(1).coords(3,:))
        if i<=numel(object_relative2V(1).coords(3,:))/2
            k=i;
        else
            k=i-numel(object_relative2V(1).coords(3,:))/2;
        end
        phase_color=[k/(numel(object_relative2V(1).coords(3,:))/2) 1-k/(numel(object_relative2V(1).coords(3,:))/2) 1];
        if i>1
            plot([object_relative2V(1).coords(1,i-1) object_relative2V(1).coords(1,i)],[object_relative2V(1).coords(3,i-1) object_relative2V(1).coords(3,i)],'Color',phase_color,'LineWidth',2,'LineStyle',line_style); %plot 2D trajectory. color it according to Vxy
        end
%         scatter(body_center(1,i),body_center(3,i),25,phase_color,'filled')
        G(i)=scatter(object_relative2V(1).coords(1,i),object_relative2V(1).coords(3,i),35,'filled','kd'); %plot 2D trajectory. color it according to Vxy
        uistack(G(i), 'top')
    end
%     scatter(object_relative2V(1).coords(1,1),object_relative2V(1).coords(3,1),'h','MarkerEdgeColor',color,'MarkerFaceColor',color) %mark 1st cycle beginging  
%     scatter(object_relative2V(1).coords(1,second_clap_frame(trial)),object_relative2V(1).coords(3,second_clap_frame(trial)),'h','MarkerEdgeColor',color,'MarkerFaceColor',color) %mark 2nd cycle beginging

%     lgnd1(trial).data=temp; %store for legend
    xlim([-0.001 0.003])
    axis equal
    set(gca,'FontSize',20);
    xlabel('$$\hat{x}$$ (m)','Interpreter','Latex', 'FontSize',24)%Body velocity axis
    ylabel('$$\hat{y}$$ (m)','Interpreter','Latex', 'FontSize',24)%Vertical perpendicular axis
    cd(pth); %return active directory to the folder that contains the kine solution files
    clearvars -except pth relative_Vxy frames second_clap_frame trials_folders switch_wings tcorr Hz continuous lgnd1 lgnd2 leveled centered%delete variables from this trial
end
% figure(leveled)
% axis equal
% legend([lgnd1(:).data],['Slow'],['Median'],['Fast'],['Fastest']); %create legend that differentiate between ventral face up to ventral face down
% set(gca,'FontSize',20);
% xlabel('$$\hat{x}$$ (m)','Interpreter','Latex', 'FontSize',24)%Body velocity axis
% ylabel('$$\hat{y}$$ (m)','Interpreter','Latex', 'FontSize',24)%Vertical perpendicular axis
% figure(centered)
% axis equal
% legend([lgnd2(:).data],['Slow'],['Median'],['Fast'],['Fastest']); %create legend that differentiate between ventral face up to ventral face down
% set(gca,'FontSize',20);
% xlabel('x (m)', 'FontSize',24)%Body velocity axis
% ylabel('y (m)', 'FontSize',24)%Vertical perpendicular axis

% print(gcf, 'amplitude', '-dtiff', '-r600');