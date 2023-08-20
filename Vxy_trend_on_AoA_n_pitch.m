%===========================================================================================%
% By Amir Sarig 16/03/2022
% - Calculate Angle-of-Attack & force components
% - Plot body pitch & angle of attack angles of all input cycles, while color coded by Vxy
% Called from 'extract_data_wasp_forward_flight' 
%===========================================================================================%

function [abs_AoA,AoA_with_direction,force_components,angles_storage] = Vxy_trend_on_AoA_n_pitch(Vxy,angles_storage,interpolant_flap_cycle_fraction,coords,Vtip_3d)
relative_Vxy=Vxy-min(Vxy); %normalize body horizonta velocity so that minimal value is zero
relative_Vxy=relative_Vxy./max(relative_Vxy); %normalize body horizonta velocity so that maximal value is one
pitch_plot=figure; hold on;
Aoa_plot=figure; hold on;
force_components=[];
% figure
% hold on
legend_in_next_cycle=0; legend_activated=0;
for wing_cycle=1:1:(length(angles_storage)/2)
    %% calc AoA
    for dig_point=1:size(coords(wing_cycle).interpolation,2)
        interpolation(dig_point).interpolation=coords(wing_cycle).interpolation(dig_point).data; %extract interpolated data
    end
    for side=1:2 %once for right, one for left
        if side==1 %left, exctract data
            base=interpolation(3).interpolation(:,:); tip=interpolation(4).interpolation(:,:); lead=interpolation(7).interpolation(:,:); trail=interpolation(9).interpolation(:,:);
            V_mid_chord=coords(wing_cycle).V_mid_chord.left';
        else %right, exctract data
            base=interpolation(5).interpolation(:,:); tip=interpolation(6).interpolation(:,:); lead=interpolation(8).interpolation(:,:); trail=interpolation(10).interpolation(:,:);
            V_mid_chord=coords(wing_cycle).V_mid_chord.right';
        end
        mid_chord=(lead'+trail')./2; %calc mid chord at ~0.7 wing length
        chord_vector=lead'-trail';
        lead_vector=lead'-base';
        trail_vector=trail'-base';
        span_vector=tip'-base';
        body_vector=interpolation(2).interpolation(:,:)'-interpolation(1).interpolation(:,:)';
        for frame=1:length(mid_chord) %calc each frame separately
            chord_angle(side)=atan2d(chord_vector(frame,3),norm(chord_vector(frame,1:2))); %store wing chord angle relative to ground for each wing
            chord_size(side)=sqrt((chord_vector(frame,1)^2)+(chord_vector(frame,2)^2)+(chord_vector(frame,3)^2));
            span_size(side)=sqrt((span_vector(frame,1)^2)+(span_vector(frame,2)^2)+(span_vector(frame,3)^2));
            body_size(side)=sqrt((body_vector(frame,1)^2)+(body_vector(frame,2)^2)+(body_vector(frame,3)^2));
            if side==1 %left
                wing_normal=cross(lead_vector(frame,:),trail_vector(frame,:)); %dorsal perpendicular to wing plane
                wing_n_velocity_normal=cross(V_mid_chord(frame,:)./norm(V_mid_chord(frame,:)),span_vector(frame,:)); %dorsal perpendicular to the plane of the wing and the wing velocity direction
            else %right
                wing_normal=cross(trail_vector(frame,:),lead_vector(frame,:)); %dorsal perpendicular to wing plane
                wing_n_velocity_normal=cross(span_vector(frame,:),V_mid_chord(frame,:)./norm(V_mid_chord(frame,:))); %dorsal perpendicular to the plane of the wing and the wing velocity direction
            end
            if abs(90-acosd(dot(chord_vector(frame,:),span_vector(frame,:))./(norm(chord_vector(frame,:))*norm(span_vector(frame,:)))))<5 && abs(90-acosd(dot(chord_vector(frame,:),wing_normal)./(norm(chord_vector(frame,:))*norm(wing_normal)))<5) %test if vector composing the rotation matrix are orthogonal to each other (plus/minus 5 degrees)
                rotation_matrix=[chord_vector(frame,:)./norm(chord_vector(frame,:));
                    span_vector(frame,:)./norm(span_vector(frame,:));
                    wing_normal./norm(wing_normal)];
            else
                errordlg('Vector composing the rotation matrix are not orthogonal to each other. *Divietion of mote then 5 degrees')
            end
            rotated_wing_n_velocity_normal=rotation_matrix*wing_n_velocity_normal';
            if side==1 %store left wing AoA
                abs_AoA(wing_cycle,frame)=acosd(dot(wing_normal,wing_n_velocity_normal)./(norm(wing_normal)*norm(wing_n_velocity_normal))); %AoA=acos((n1*n2)/(|n1|*|n2|))
                AoA_with_direction(wing_cycle,frame)=atan2d(rotated_wing_n_velocity_normal(1),rotated_wing_n_velocity_normal(3)); %angle between wing_normal to wing_n_velocity_normal
            else % calc mean of both wing's AoA
                abs_AoA(wing_cycle,frame)=mean([abs_AoA(wing_cycle,frame),acosd(dot(wing_normal,wing_n_velocity_normal)./(norm(wing_normal)*norm(wing_n_velocity_normal)))]); 
                AoA_with_direction(wing_cycle,frame)=mean([AoA_with_direction(wing_cycle,frame),atan2d(rotated_wing_n_velocity_normal(1),rotated_wing_n_velocity_normal(3))]);
                if any(frame==(0:10:length(mid_chord))) %for each 10th frame
                    force_components(size(force_components,1)+1,:)=[wing_cycle,frame,Vtip_3d(wing_cycle,frame),mean(chord_angle),AoA_with_direction(wing_cycle,frame),mean(chord_size),mean(span_size),mean(body_size)]; %store: cycle, frame, wingtip ground speed, chord angle, AoA, chrod size, wing length, body length
                end
            end
        end 
        angles_storage(2*wing_cycle-1).AoA=AoA_with_direction(wing_cycle,:);
        angles_storage(2*wing_cycle-1).abs_AoA=abs_AoA(wing_cycle,:);
    end
%% plot 
%     subplot(1,2,2) %body pitch
    figure(pitch_plot);
    if angles_storage(2*wing_cycle-1).pitch.ventral_face_down
        line_style='-'; %plot continous curve if ventral_face_down
        pitch_data=angles_storage(2*wing_cycle-1).pitch.angle; %exctract pitch data
    else
        line_style=':'; %plot dotted curve if ventral_face_up
%         pitch_data=angles_storage(2*wing_cycle-1).pitch.angle.*(-1); %reverse angle to be aligned with data from trials with ventral side face down
        pitch_data=angles_storage(2*wing_cycle-1).pitch.angle; %
    end
    temp=plot(interpolant_flap_cycle_fraction,pitch_data,'color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineStyle',line_style,'LineWidth',2);
    clear pitch_data
%     hold on
    lgnd(wing_cycle).data=temp;
    if ~legend_activated %if legend haven't been activated yet
        if ~angles_storage(2*wing_cycle-1).pitch.ventral_face_down && wing_cycle>1 %if ventral face down in this cycle and it isn't the 1st cycle
           legend([lgnd((wing_cycle-1):wing_cycle).data],['ventral face down'],['ventral face up']); %create legend that differentiate between ventral face up to ventral face down
           legend_activated=1;
        elseif ~angles_storage(2*wing_cycle-1).pitch.ventral_face_down && wing_cycle==1 %if ventral face down in this cycle but it is the 1st cycle
           legend_in_next_cycle=wing_cycle; %store flag to create legend one there will be cycle that ventral face up
        elseif legend_in_next_cycle && angles_storage(2*wing_cycle-1).pitch.ventral_face_down %if ventral face down was previously found and in this cycle ventral face up
            legend([lgnd([legend_in_next_cycle,wing_cycle]).data],['ventral face up'],['ventral face down']); %create legend that differentiate between ventral face up to ventral face down
            legend_activated=1;
        end
    end
%     subplot(1,3,2) %AoA
%     plot(interpolant_flap_cycle_fraction,AoA(wing_cycle,:),'color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineStyle',line_style,'LineWidth',2);
%     hold on
%     subplot(1,2,1) %AoA with +/- direction
    figure(Aoa_plot);
    plot(interpolant_flap_cycle_fraction,AoA_with_direction(wing_cycle,:),'color',[relative_Vxy(wing_cycle) 0 1-relative_Vxy(wing_cycle)],'LineStyle',line_style,'LineWidth',2);
    hold on
end
%% plot annotations
% subplot(1,2,2)
figure(pitch_plot)
title('Body Pitch');
xlabel('Flap cycle fraction') %X axis label 
ylabel('Angle (degree)') %Y axis label for 1st figure in the plot implay for all the row
% ylim([-90 90])
% subplot(1,2,1)
figure(Aoa_plot)
title('Angle of Attack');
xlabel('Flap cycle fraction') %X axis label 
ylabel('Angle (degree)') %Y axis label for 1st figure in the plot implay for all the row
% pos = get(subplot(1,2,2),'Position');
% cbh =colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)  0.02  pos(4)]); %Create Colorbar
figure
cbh =colorbar;
cbh.Ticks = 0: 1/10: 1;  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.01:1); zeros(1,101); (1:-0.01:0)]'; %create matrix of RGB colormap: Low-Blue, High-Red
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(min(Vxy), max(Vxy), 11) ; %create labels 
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.3f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
cycles=unique(force_components(:,1)); %extract cycles names
for wing_cycle=1:length(cycles)
    i=force_components(:,1)==cycles(wing_cycle); %extract rows index of all data from current wing cycle
    mean_force_components(wing_cycle,:)=mean(force_components(i,:),1); %calc wing cycle's mean value of each parameter
end
mean_force_components(:,2)=[]; %delete 2nd column  that contain frame number data
force_components=array2table(mean_force_components,'VariableNames',{'cycle', 'wingtip_ground_speed', 'chord_angle', 'AoA', 'chrod_size', 'wing_length', 'body_length'}); %convert to table variable
end