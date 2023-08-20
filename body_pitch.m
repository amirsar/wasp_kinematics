% Called from mean_V_from_kine_solution function

function pitch=body_pitch(pitch,trial,cycle,index,head,abd,Ltip,Rtip,elevation,tcorr,draw_plot)
%% pitch calculation
pitch(index).trial=trial; %store trial
pitch(index).cycle=cycle; %store cycle
body_center=((head+abd)./2); % body center
normalize_body_vector=(head-body_center)-(abd-body_center); % normalized body vector
normalize_wingtip=(((Ltip-body_center)+(Rtip-body_center))./2); % normalized mean left&right wingtip
frame=round(0.9*size(head,2)); %frame during clap
transverse=cross(normalize_wingtip(:,frame),normalize_body_vector(:,frame)); %body transverse axis, from right to left of the body: cross(wings length,body_vector) 
dorsal_norm=cross(normalize_body_vector(:,frame),transverse); %dorsal perpendicular to the plane of the body
if dorsal_norm(3)>=0 %if body's dorsal perpendicular is above body center
    pitch_sign=(1); %will leave pitch angle to range of -90 till 90. meaning ventral_face_ground
    pitch(index).ventral_face_down=1; %define ventral side as facing down
else %if dorsal perpendicular is bellow body center
    pitch_sign=(-1); %will pass pitch angle to range of |90| till |180|. meaning dorsal_face_ground
    pitch(index).ventral_face_down=0; %define ventral side as facing up for this trial
end
% angle=atan2d(normalize_body_vector(3,:),pitch_sign*sqrt((normalize_body_vector(1,:).^2)+(normalize_body_vector(2,:).^2))); %pitch calculation
angle=atan2d(normalize_body_vector(3,:),sqrt((normalize_body_vector(1,:).^2)+(normalize_body_vector(2,:).^2))); %pitch calculation
% for frame=1:length(angle)
%     if angle(frame)<0 %if pitch is negative
%         angle(frame)=angle(frame)+360; %add 360. will shift range from [-180 180] to [0 360], in order to keep contiousity
%     end
% end

%% Stroke plae
% SP relative to body
normalized_flapping_plane=normalize_wingtip(:,round(size(normalize_wingtip,2)*0.9))-normalize_wingtip(:,round(size(normalize_wingtip,2)*0.5)); %wingtip position at 90% of the cycle (clap) minus position at 50% (ventral reversal point)
normalize_flapping_angle_from_body=acosd(dot(mean(normalize_body_vector,2),normalized_flapping_plane)./(norm(mean(normalize_body_vector,2))*norm(normalized_flapping_plane))); %angle between mean body position to flapping palne: acos(dot(mean(body),flapping_plane)/(|mean(body)|*|flapping_plane|))
normalize_flapping_angle_from_horizon=mean(angle)+normalize_flapping_angle_from_body;    
if normalize_flapping_angle_from_horizon>180 %change angles range from [0 360] to absolute[-180 180]
    abs_normalize_flapping_angle_from_horizontal_direction=abs(normalize_flapping_angle_from_horizon-360); 
else
    abs_normalize_flapping_angle_from_horizontal_direction=normalize_flapping_angle_from_horizon;
end
if abs(mean(angle))>90 %meaning ventral facing up, and horizontal dirction is towards 180 not 0 degrees
    abs_normalize_flapping_angle_from_horizontal_direction=180-abs_normalize_flapping_angle_from_horizontal_direction; %reverse direction
end
pitch(index).stroke_plane_relative2body.flapping_angle_from_horizon=normalize_flapping_angle_from_horizon;
pitch(index).stroke_plane_relative2body.abs_flapping_angle_from_horizontal_direction=abs_normalize_flapping_angle_from_horizontal_direction;
pitch(index).stroke_plane_relative2body.flapping_angle_from_body=normalize_flapping_angle_from_body;

% SP relative to cameras
% mean_body=mean(head-abd,2);
% wingtip=(Ltip+Rtip)./2;
% flapping_plane=wingtip(:,round(size(wingtip,2)*0.9))-wingtip(:,round(size(wingtip,2)*0.5)); %wingtip position at 90% of the cycle minus position at 50%
% transverse=cross(flapping_plane,mean_body); %body transverse axis, from right to left of the body: cross((clap-ventral_reversal),mean(head-abd)) 
% dorsal_norm=cross(flapping_plane,transverse); %dorsal perpendicular to the flapping_plane: cross(flapping_plane, transverse) 
% if dorsal_norm(3)>=(mean_body(3)/2) %if body's dorsal perpendicular is above body center
%     flapping_plane_sign=(1); %will leave pitch angle to range of 0 till 90
% else %if dorsal wings' clap is bellow body center
%     flapping_plane_sign=(-1); %will pass pitch angle to range of 90 till 180 
% end
% flapping_angle_from_horizon=atan2d(flapping_plane(3),flapping_plane_sign*sqrt((flapping_plane(1)^2)+(flapping_plane(2)^2))); 
% for frame=1:length(flapping_angle_from_horizon)
%     if flapping_angle_from_horizon(frame)<0 %if SP is negative
%         pitch(index).stroke_plane.flapping_angle_from_horizon(frame)=flapping_angle_from_horizon(frame)+360; %add 360. will shift range from [-180 180] to [0 360], in order to keep contiousity
%     end
% end
% pitch(index).stroke_plane.abs_flapping_angle_from_horizon=abs(flapping_angle_from_horizon-180);
% pitch(index).stroke_plane.flapping_angle_from_body=flapping_angle_from_horizon-mean(angle);


%% extend wing cycle to unified length
rr = length(angle);
t = (1./rr:1./rr:1.0)'; %X values (breakpoints) between 0 to 1 with data (Y value) in the wing cycle
cs=spline(t,[angle(1) angle(1,:) angle(end)]); % piecewise polynomial form of the cubic spline interpolant of cycle_data for later use with the first and last value in cycle_data are used as the endslopes for the cubic spline.
interpolant_angle = ppval(cs,tcorr); %Spline Interpolation with Specified Endpoint Slopes
interpolant_angle(1)=angle(1); interpolant_angle(end)=angle(end); %ensure that first and last value in cycle_data are used
pitch(index).angle=interpolant_angle';

elevation=mean(elevation,1); %mean of both wings
rr = length(elevation);
t = (1./rr:1./rr:1.0)'; %X values (breakpoints) between 0 to 1 with data (Y value) in the wing cycle
cs=spline(t,[elevation(1) elevation(:)' elevation(end)]); % piecewise polynomial form of the cubic spline interpolant of cycle_data for later use with the first and last value in cycle_data are used as the endslopes for the cubic spline.
interpolant_elevation = ppval(cs,tcorr)'; %Spline Interpolation with Specified Endpoint Slopes
interpolant_elevation(1)=elevation(1); interpolant_elevation(end)=elevation(end); %ensure that first and last value in cycle_data are used

%% downstroke & upstroke classification
[val_min,loc_min]=min(interpolant_elevation); %downstroke minimum point & it's location
[val_max,loc_max]=max(interpolant_elevation); %upstroke maximum point & it's location 
%test for common mistakes in stroke direction classification
if loc_max<loc_min || loc_max<length(interpolant_elevation)*0.666 || loc_min<length(interpolant_elevation)*0.25
    %upstroke max point comes before downstroke min point || upstroke
    %max point comes before the end of the cycle || min point comes too early
    s=sprintf('upstroke max point comes before downstroke min point, \nor upstroke max point comes before the last third of the cycle, \nor min point comes before the second quarter of the cycle');
    choice = questdlg(s, 'Min\Max Error', 'Manually correct','Leave it that way','Leave it that way');
    if strcmp(choice,'Manually correct') %if user choose to manually mark min & max points
        fprintf('Error in elevation angle of trial %d cycle %d.\n',pitch(index).trial,pitch(index).cycle); %print to command window
        figure
        plot(interpolant_elevation) %plot elevation
        tit=sprintf('Trial %d cycle %d.\n Red: Min. Green: Max',pitch(index).trial,pitch(index).cycle); %title 
        title(tit);
        hold on
        plot(loc_min,interpolant_elevation(loc_min),'o','color','r') %plot rejected locations
        plot(loc_max,interpolant_elevation(loc_max),'o','color','g') %plot rejected locations
        msg=sprintf('Error in elevation angle.\nChoose new min point and then max point.\n Consider to remove this cycle from sample.');
        uiwait(msgbox(msg));
        [x,y,button] = ginput(2); %recieve new location coordinates from user
        close %close insect image
        loc_min=round(x(1)); loc_max=round(x(2)); %store user choise
    end
end
pitch(index).wing_cycle=[loc_min loc_max]; %1-downstroke 2-upstroke

%% plot
if draw_plot %if it is the last cycle
    %create super-title for the following subplot
    f=figure;
    p = uipanel('Parent',f,'BorderType','none'); 
    tit=sprintf('Periodic changes in body pitch (Gray is downstroke)');
    p.Title = tit; 
    p.TitlePosition = 'centertop'; 
    p.FontSize = 12;
    p.FontWeight = 'bold';
    trials=unique(pitch(1).trial(:)); %extract trials numbers
    for i=1:length(trials) %repeat for each trial
        if length(trials)>5
            if length(trials)>10
                number_of_columns=3;
            else
                number_of_columns=2;
            end
        else
            number_of_columns=1;
        end
        subplot(ceil(length(trials)/number_of_columns),number_of_columns,i,'Parent',p); %slote for each trial
        cycles=trials(i)==pitch(1).trial(:); %extract row naumbers of cycles in this trial
        cycles=find(cycles); %switch to index numbers
        pitch_angle=[]; %reset variable
        for j=1:length(cycles) %repeat for each wing cycle
            pitch_angle((length(pitch_angle)+1):(length(pitch_angle)+length(pitch(cycles(j)).angle)))=pitch(cycles(j)).angle; %extract trial's pitch angles, cycle after cycle
        end
        sum_previous_length=0; %will be used to place background relative to it's wing cycle
        plot(pitch_angle) %plot in order to recieve axis limits
        axis tight
        hold on
        for j=1:length(cycles) %set gray backround while in downstroke
            patch(([0 pitch(cycles(j)).wing_cycle(1) pitch(cycles(j)).wing_cycle(1) 0]+sum_previous_length),[max(ylim) max(ylim) min(ylim) min(ylim)],[0.7 0.7 0.7]) %wing cycle begins with downstroke
            patch(([pitch(cycles(j)).wing_cycle(2) length(pitch(cycles(j)).angle) length(pitch(cycles(j)).angle) pitch(cycles(j)).wing_cycle(2)]+sum_previous_length),[max(ylim) max(ylim) min(ylim) min(ylim)],[0.7 0.7 0.7]) %wing cycle ends with downstroke
            sum_previous_length=sum_previous_length+length(pitch(cycles(j)).angle); %will be used to place next iteration background
        end
        plot(pitch_angle) %plot again because the background cover the previous line
        tit=sprintf('Trial %d',i);
        title(tit);
        set(gca,'XTick',0:length(tcorr):length(cycles)*length(tcorr))
        set(gca,'XTickLabel',0:length(cycles))
        if any(i==(ceil(length(trials)/number_of_columns)*number_of_columns)-number_of_columns+1:ceil(length(trials)/number_of_columns)*number_of_columns) %add xlabel for the buttom of each column in the subplot
            xlabel('Wing cycles')
        end
        if any(i==(0:1:ceil(length(trials)/number_of_columns))) %add ylabel for the left figure of each row in the subplot
        ylabel('Pitch (degree)')
        end
    end
end
