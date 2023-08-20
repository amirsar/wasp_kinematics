clear;
close all;
current_path=pwd; %recieve current path and save it in order to get back to it at the end of the function.

%% ask user for information
flapping_cycles_counter=0;
trials_counter=0;
climb_range_of_interest=[0,10]; %set whice data will be stored
file_path=uigetdir('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Digitizing\Wind tunnel\Horizontal\MSc_E.mundus\withFlappingPlane','Choose the folder containing tracks data:'); %ask the user to enter the path of the data.
cd(file_path); %change the folder that MATLAB is using to the folder with the files.
speeds=dir; %recieve data of the items in the folder.
num_WindSpeeds=sum([speeds(~ismember({speeds.name},{'.','..'})).isdir]); %count the number of sub-folders in the folder.
Hz=5000; %setting the frame rate of the record
angles=[]; %reset index row for data summary table
%% extract data
for WindSpeed=1:num_WindSpeeds %repeat as many times as items in the folder (wind speeds)
    wind_path=fullfile(file_path,speeds(WindSpeed+2).name); %create the full path of the wind speed folder.
    cd(wind_path); %change the folder that MATLAB is using to the folder with the data of trials in this wind speed.
    num_files=size(dir); %count the number of files in the folder.
    num_files=num_files(1)-2; %make a correction to the previous count, becuse the two first objects are not folders.
    folder=dir; %recieve information about the items in the folder.
    folder=folder(3:end,:); %delete the first two lines, because they are not folders.
    %calculate and store 3D points (p1-posterior, p2-head, p3-clap, p4-ventral reversal point)
    for trial=1:num_files %repeat as many times as items in the folder (trials)
        negative_pitch_flap=0; %will be used to mark trials with negative pitch
        path=fullfile(file_path,speeds(WindSpeed+2).name,'\',folder(trial).name); %create the full path of the folder.
        cd(path); %change the folder that MATLAB is using to the folder with the files.
        data_folder=dir; %recieve information about the files in the folder.
        xyz=xlsread(data_folder(5).name); %load to "file" the file with coordinates matrix.
        [r,c]=size(xyz);  
        xyz=xyz./100; %converting the input data from units is cm to meter
        body_center=((xyz(:,1:3)+xyz(:,4:6))./2);
        flight=(velocity(body_center,Hz));
        p1_normalize=xyz(:,1:3)-body_center;
        p1=xyz(:,1:3);
        p2_normalize=xyz(:,4:6)-body_center;
        p2=xyz(:,4:6);
        current_speed=(str2num(speeds(WindSpeed+2).name))*(-1);
        flight(:,4)=(flight(:,1))-current_speed; %convert ground speed to movement speed
        flight(:,5) = sqrt(flight(:,1).^2+flight(:,2).^2); %find the diagonal which is the Horizontal speed (in xy plain), based on ground speed
        flight(:,6) = sqrt(flight(:,4).^2+flight(:,2).^2); %find the diagonal which is the Horizontal speed (in xy plain), based on movement speed
%% calc body pitch        
        clap_frames=find(~isnan(xyz(:,7))); %all frames with clap
        for idx=1:r %repeat for each row (frame) in the file (trial)
            if ~isnan(xyz(idx,7)) %calculate only in frames of wings clap
                p3=xyz(idx,7:9); % dorsal clap
                p3_normalize=xyz(idx,7:9)-body_center(idx,:); %normalized dorsal clap
                transverse=cross(p3_normalize,(p2_normalize(idx,:)-p1_normalize(idx,:))); %body transverse axis, from right to left of the body: cross(clap,(head-abd)) 
                dorsal_norm=cross((p2_normalize(idx,:)-p1_normalize(idx,:)),transverse); %dorsal perpendicular to the plane of the body: cross((head-abd), transverse) 
                if dorsal_norm(3)>=0 %if body's dorsal perpendicular is above body center
                    ventral_face_ground=1; %mark as true
                    pitch_sign=(1); %will leave pitch angle to range of 0 till 90
                else %if dorsal wings' clap is bellow body center
                    ventral_face_ground=0; %mark as false
                    pitch_sign=(-1); %will pass pitch angle to range of 90 till 180 
                end
                dorsal_clap_frame=idx;
            end
            if ~isnan(xyz(idx,10)) %calculate only in frames of wings at ventral position
                if clap_frames((dorsal_clap_frame==clap_frames))~=clap_frames(end) %if it's not the last clap
                    flapping_cycles_counter=flapping_cycles_counter+1;
                    next_clap=clap_frames(find(dorsal_clap_frame==clap_frames)+1); %extract frame number of next dorsal clap
                    cycle_frq=1/((next_clap-dorsal_clap_frame)/Hz);
                    Vxy_ground=mean(flight(dorsal_clap_frame:next_clap,5)); %mean of value during flapping cycle
                    Vxy_arial=mean(flight(dorsal_clap_frame:next_clap,6)); %mean of value during flapping cycle
                    Vz=mean(flight(dorsal_clap_frame:next_clap,3)); %mean of value during flapping cycle
                    climb=atan2d(Vz,Vxy_arial); %mean vertical velocity angle
                    p4=xyz(idx,10:12); %ventral clap
                    p4_normalize=xyz(idx,10:12)-body_center(idx,:); %normalized ventral clap
                    flapping_plane=p3-p4; %p3-clap, p4-ventral
                    normalized_flapping_plane=p3_normalize-p4_normalize; %p3-clap, p4-ventral
                    amplitude=norm(normalized_flapping_plane);
                    body=p2_normalize(dorsal_clap_frame:next_clap,:)-p1_normalize(dorsal_clap_frame:next_clap,:); %normalized mean body vector
                    pitch=atan2d(body(:,3),pitch_sign*sqrt((body(:,1).^2)+(body(:,2).^2)));
                    for frame=1:length(pitch) %for each frame
                        if pitch(frame)<0 %if pitch is negative
                            negative_pitch_flap=1; %will mark this trial
                            pitch(frame)=pitch(frame)+360; %add 360. will shift range from [-180 180] to [0 360], in order to keep contiousity
                        end
                    end
%% calc Stroke Plane                    
% SP relative to cameras
                    mean_body=mean(p2(dorsal_clap_frame:next_clap,:)-p1(dorsal_clap_frame:next_clap,:),1);
                    transverse=cross(flapping_plane,mean_body); %body transverse axis, from right to left of the body: cross((clap-ventral_reversal),mean(head-abd)) 
                    dorsal_norm=cross(flapping_plane,transverse); %dorsal perpendicular to the flapping_plane: cross(flapping_plane, transverse) 
                    if dorsal_norm(3)>=(mean_body(3)/2) %if body's dorsal perpendicular is above body center
                        flapping_plane_sign=(1); %will leave pitch angle to range of 0 till 90
                    else %if dorsal wings' clap is bellow body center
                        flapping_plane_sign=(-1); %will pass pitch angle to range of 90 till 180 
                    end
                    flapping_angle_from_horizon=atan2d(flapping_plane(3),flapping_plane_sign*sqrt((flapping_plane(1)^2)+(flapping_plane(2)^2))); 
                    for frame=1:length(flapping_angle_from_horizon)
                        if flapping_angle_from_horizon(frame)<0 %if SP is negative
                            flapping_angle_from_horizon(frame)=flapping_angle_from_horizon(frame)+360; %add 360. will shift range from [-180 180] to [0 360], in order to keep contiousity
                        end
                    end
                    abs_flapping_angle_from_horizon=abs(flapping_angle_from_horizon-180);
                    flapping_angle_from_body=flapping_angle_from_horizon-mean(pitch);
% SP relative to body                  
                    normalize_flapping_angle_from_body=acosd(dot(mean(body,1),normalized_flapping_plane)./(norm(mean(body,1))*norm(normalized_flapping_plane))); %angle between mean body position to flapping palne: acos(dot(body,flapping_plane)/(|body|*|flapping_plane|))
                    normalize_flapping_angle_from_horizontal_direction=mean(pitch)+normalize_flapping_angle_from_body;
                    if normalize_flapping_angle_from_horizontal_direction>180 %change angles range from [0 360] to absolute[-180 180]
                        abs_normalize_flapping_angle_from_horizon=abs(normalize_flapping_angle_from_horizontal_direction-360); 
                    else
                        abs_normalize_flapping_angle_from_horizon=normalize_flapping_angle_from_horizontal_direction;
                    end
% %olde SP relative to body - angle is in the opposite direction (>90)
%                     normalize_flapping_angle_from_horizon=atan2d(normalized_flapping_plane(3),sqrt((normalized_flapping_plane(1)^2)+(normalized_flapping_plane(2)^2))); %wrong clac: 1) horizontal component is size not vector, therefor the angle is limited up to |90|. 2) angles beyong |90| are meaningless because the axes aren't oriented to the body orientations
%                     normalize_flapping_angle_from_body=normalize_flapping_angle_from_horizon-mean(pitch);
%% report data                    
                    if climb>=climb_range_of_interest(1) && climb<=climb_range_of_interest(2) %store data only if climb is in range of interest
                        if mean(pitch)<=90 %store data only if mean pitch is up to 90 degree (only ventral body side facing ground cycles)
                            if negative_pitch_flap %if final sample include trials that their negative pitch angle was increased by 360 degrees
                                fprintf('In trial %d, a cycle had negative pitch that was increased by 360 degrees\n',str2double(folder(trial).name(5:end))) %print trial name
                            end
                            if ~ventral_face_ground %if final sample include trial with ventral body side face upwards
                                fprintf('%d\n',str2double(folder(trial).name(5:end))) %print trial name
                            end
                            angles(1+size(angles,1),:)=[str2double(speeds(WindSpeed+2).name)*(-1),str2double(folder(trial).name(5:end)),idx,ventral_face_ground,Vxy_ground,Vxy_arial,Vz,climb,mean(pitch),median(pitch),normalize_flapping_angle_from_body,abs_normalize_flapping_angle_from_horizon,normalize_flapping_angle_from_horizontal_direction,cycle_frq,amplitude]; %store data
                        end
                        clearvars pitch
                    end
                end
            end
        end
    end
    trials_counter=trials_counter+trial;
end
table = array2table(angles,'VariableNames',{'wind_speed', 'trial', 'frame','ventral_face_down','mean_Vxy_ground','mean_Vxy_movement','mean_Vz','mean_vertical_velocity_angle','mean_body_pitch','median_body_pitch','normalize_flapping_angle_from_mean_body_pitch','abs_normalize_flapping_angle_from_horizon','normalize_flapping_angle_from_horizon','cycle_frq','amplitude'});
cd(current_path);
%% plot
variable=[9,13,11];
variable_name={'mean Body Pitch (degree)'; 'Stroke Plane from Horizon(degree)'; 'Stroke Plane from mean Body Pitch(degree)'};
for i=1:length(variable)
    figure
    plot(table{:,variable(i)},table{:,6},'ko','markerfacecolor','k')
    ylabel('Horizontal aerial speed (m/s)','fontsize',20)
    xlabel(variable_name(i),'fontsize',20)
    xlim([0 200])
end
fprintf('Out of %d flapping cycle, from %d trial - Only %d flapping cycle from %d trials were enter to the sample. i.e. were with positive climb angle of up to 10 degrees and ventral body side facing ground\n',flapping_cycles_counter,trials_counter, size(angles,1), length(unique(table{:,2})))
%%
trials=unique(table{:,2}); %extract trials names
c=1:size(table,2); %extract data from all columns
c(3)=[]; %except frames
for trial=1:length(trials) %for each trial
    mean_data(trial,:)=mean(table{trials(trial)==table{:,2},c},1); %calc mean for each variable
end
mean_data = array2table(mean_data,'VariableNames',{'wind_speed', 'trial','ventral_face_down','mean_Vxy_ground','mean_Vxy_movement','mean_Vz','mean_vertical_velocity_angle','mean_body_pitch','median_body_pitch','normalize_flapping_angle_from_mean_body_pitch','abs_normalize_flapping_angle_from_horizon','normalize_flapping_angle_from_horizon','cycle_frq','amplitude'});