function Force=wingtip_world_viewpoint(wing_cycle,Vxy,coords,object_relative2V,Force)
load('D:\GalR10\Desktop\OneDrive - Tel-Aviv University\Publications\Wasp kinematics article\Data\CFD\Thomas'' data\instantaneous_force\only force_data\CFD_instantaneous_F.mat')
figure
hold on
if wing_cycle==11 %execute only for the last cycle, in order to return it to intiall orientation of dorsal face ground
    body_center=((object_relative2V(4).coords+object_relative2V(5).coords)./2); % rotated body center
    for object=1:size(object_relative2V,2) %repeat for each object. order: 1)wingtip 2)lead 3)trail ...
        for frame=1:size(object_relative2V(object).coords,2) %rotate each frame's wingtip coordinate
            shifted=object_relative2V(object).coords(:,frame)-body_center(:,frame); %set body center as (0,0) because rotation is around the beginning of the axes
            x_rotated=([1 0 0; 0 cos(pi) -sin(pi); 0 sin(pi) cos(pi)]*shifted)';% rotated by 180 degrees, from dorsal face up to dorsal face ground
            object_relative2V(object).coords(:,frame)=x_rotated'+body_center(:,frame); %store in a way that each axis is a row, return to original position
        end
    end
end
% V_mid_chord=(coords(wing_cycle).V_mid_chord.left'+coords(wing_cycle).V_mid_chord.right')./2; %mid-chord vector speed, L&R wings mean 
V_mid_chord=object_relative2V(6).coords';
%% calc normalize wing speed
% body_center=(coords(wing_cycle).interpolation(1).data(:,:)+coords(wing_cycle).interpolation(2).data(:,:))./2; %calc body center 
for frame=1:length(V_mid_chord) %for each frame
    V3D(frame,:)=norm(V_mid_chord(frame,:)); % calc 3D wing speed
    V3D_unit_vector(frame,:)=V_mid_chord(frame,:)./V3D(frame); %create unit vector of wing speed
end
relative_V3D=V3D/max(V3D); %normalize 3D wing speed to be between 0 to 1
chords_size=mean([norm(coords(wing_cycle).interpolation(7).data(:,1)-coords(wing_cycle).interpolation(9).data(:,1)); norm(coords(wing_cycle).interpolation(8).data(:,1)-coords(wing_cycle).interpolation(10).data(:,1))]); %size of chord at 1st frame (L&R wings mean)
V_mid_chord=V3D_unit_vector.*repmat(relative_V3D*chords_size,1,3); %normalize wing speed size according to chord size and relative wing speed

%% normalize coordinates
% if body_center(1,1)>body_center(1,end) %if flight direction is towards negative X axis
%     direction=-1;
% else
%     direction=1;
% end
% V_mid_chord(:,1)=V_mid_chord(:,1).*direction; %normalize chord velocity direction towards possitive X axis
% for digitizing_point=1:size(coords(wing_cycle).interpolation,2) % repeat for each digiting point. order: 1)head;2)abd;3)Lbase;4)Ltip;5)Rbase;6)Rtip;7)Llead;8)Rlead;9)Ltrail;10)Rtrail;
%     normalize(digitizing_point).coords=(coords(wing_cycle).interpolation(digitizing_point).data-repmat(body_center(:,1),1,50))'; % start each cycle at (0,0,0), (all frames of each digitizing_point minus body center at 1st frame) 
%     normalize(digitizing_point).coords(:,1)=normalize(digitizing_point).coords(:,1).*direction; %normalize flight direction towards possitive X axis
% end
% for dim=1:3 %store data. repeat for each axis
%     data(1).coords(dim,:)=mean([normalize(4).coords(:,dim), normalize(6).coords(:,dim)],2); %wingtip left&right wing mean. is calculated along columns
%     data(2).coords(dim,:)=mean([normalize(7).coords(:,dim), normalize(8).coords(:,dim)],2); %leading edge. left&right wingmean is calculated along columns
%     data(3).coords(dim,:)=mean([normalize(9).coords(:,dim), normalize(10).coords(:,dim)],2); %trailing edge. left&right wingmean is calculated along columns
%     data(4).coords(dim,:)=normalize(1).coords(:,dim); %head
%     data(5).coords(dim,:)=normalize(2).coords(:,dim); %abdomen
%     data(6).coords(dim,:)=mean([normalize(1).coords(:,dim), normalize(2).coords(:,dim)],2); %body center
% end
%% extract data points
for dim=1:3 %store data. repeat for each axis
    data(1).coords(dim,:)=object_relative2V(1).coords(dim,:); %wingtip left&right wing mean. is calculated along columns
    data(2).coords(dim,:)=object_relative2V(2).coords(dim,:); %leading edge. left&right wingmean is calculated along columns
    data(3).coords(dim,:)=object_relative2V(3).coords(dim,:); %trailing edge. left&right wingmean is calculated along columns
    data(4).coords(dim,:)=object_relative2V(4).coords(dim,:); %head
    data(5).coords(dim,:)=object_relative2V(5).coords(dim,:); %abdomen
    data(6).coords(dim,:)=mean([data(4).coords(dim,:), data(5).coords(dim,:)],2); %body center
end
%% calc viewpoint rotation & chord position shift
% velocity_vector=data(6).coords(:,end)-data(6).coords(:,1); %cycle's velocity_vector is body center's last frame minus 1st frame
% horizontal_angle=atan2d(velocity_vector(2),velocity_vector(1))-90; %calculate cycle's velocity_vector projection on XY plain minus 90 degrees
% if velocity_vector(1)*velocity_vector(2)<0
%     horizontal_angle=horizontal_angle-180; %perform horizontal angle correction
% end
% vertical_angle=atan2d(velocity_vector(3),sqrt((velocity_vector(1)^2)+(velocity_vector(2)^2))); %calculate cycle's velocity_vector Vertical_angle
for dim=1:3 %calc distance between wingtip to middle of chord. repeat for each axis
    chord_shift(dim,:)=data(1).coords(dim,:)-mean([data(2).coords(dim,:); data(3).coords(dim,:)],1); %distance between middle of the cord to wingtip. mean is calculated along columns
end

%% plot
% steps=numel(data(5).coords(1,:))/10:numel(data(5).coords(1,:))/10:numel(data(5).coords(1,:)); %divid flap cycle to 10 steps
% for i=2:1:numel(data(5).coords(1,:)) %for each phase of wing cycle
%     plot3(data(1).coords(1,i-1:i),data(1).coords(2,i-1:i),data(1).coords(3,i-1:i),'Color',[i/numel(data(5).coords(1,:)) 1-i/numel(data(5).coords(1,:)) 1],'LineWidth',2) %plot mean 3D trajectory betweem this phase to the previous one. color it according to time along flapping cycle
%     if any(i==steps) %for each of the 10th fractions of the flap cycle 
%         plot3([data(4).coords(1,i), data(5).coords(1,i)],[data(4).coords(2,i), data(5).coords(2,i)],[data(4).coords(3,i), data(5).coords(3,i)],'Color',[0.5 0.5 0.5],'LineWidth',1) %plot body location.  
%         plot3(data(4).coords(1,i),data(4).coords(2,i),data(4).coords(3,i),'d','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]) %plot head location as full diamond. 
%         plot3([data(2).coords(1,i), data(3).coords(1,i)]+chord_shift(1,i),[data(2).coords(2,i), data(3).coords(2,i)]+chord_shift(2,i),[data(2).coords(3,i), data(3).coords(3,i)]+chord_shift(3,i),'Color','k','LineWidth',2) %plot chord location.  
%         plot3(data(2).coords(1,i)+chord_shift(1,i),data(2).coords(2,i)+chord_shift(2,i),data(2).coords(3,i)+chord_shift(3,i),'o','Color','k','MarkerFaceColor','k') %plot leading edge location as full cycle. 
%         quiver3(data(1).coords(1,i),data(1).coords(2,i),data(1).coords(3,i),V_mid_chord(i,1),V_mid_chord(i,2),V_mid_chord(i,3),'color',[0 0.8 0],'LineWidth',3,'autoscale','off','markersize',3,'maxheadsize',3); %plot wing velocity vector
%     end
% end 
% plot3(mean([data(4).coords(1,i), data(5).coords(1,i)]),mean([data(4).coords(2,i), data(5).coords(2,i)]),mean([data(4).coords(3,i), data(5).coords(3,i)]),'kh','MarkerFaceColor',[0.5 0.5 0.5]) %plot last body center
% view(-horizontal_angle,-vertical_angle) %set viewpoint to horizontal and perpendicular

%% plot
steps_amount=20; %set amount of steps
steps=0:numel(data(5).coords(1,:))/steps_amount:numel(data(5).coords(1,:)); %divid flap cycle to 20 segments
steps(1)=1; %set 1st step to 1st frame
scatter(mean([data(4).coords(1,1), data(5).coords(1,1)]),mean([data(4).coords(3,1), data(5).coords(3,1)]),100,'filled','kh') %mark 1st body center with hexagone
for i=1:1:numel(data(5).coords(1,:)) %for each data point of wing cycle
    if i>1 %skip 1st frame
%         plot([data(1).coords(1,i-1),data(1).coords(1,i)],[data(1).coords(3,i-1),data(1).coords(3,i)],'-.','Color',[i/numel(data(5).coords(1,:)) 1-i/numel(data(5).coords(1,:)) 1],'MarkerFaceColor',[i/numel(data(5).coords(1,:)) 1-i/numel(data(5).coords(1,:)) 1],'LineWidth',3) %plot wingtip trajectory: sagittal 2D location of this phase. color it according to time along flapping cycle
        plot([mean([mean(data(2).coords(1,i-1)), mean(data(3).coords(1,i-1))]),mean([mean(data(2).coords(1,i)), mean(data(3).coords(1,i))])],[mean([mean(data(2).coords(3,i-1)), mean(data(3).coords(3,i-1))]),mean([mean(data(2).coords(3,i)), mean(data(3).coords(3,i))])],'-.','Color',[i/numel(data(5).coords(1,:)) 1-i/numel(data(5).coords(1,:)) 1],'MarkerFaceColor',[i/numel(data(5).coords(1,:)) 1-i/numel(data(5).coords(1,:)) 1],'LineWidth',3) %plot mid-chord trajectory: sagittal 2D location of this phase. color it according to time along flapping cycle
    end
%     if i==31 %draw just selected frame on the trajectory, for easy comparision
%         k=i;
    if any(i==round(steps)) %for each of segemnt of the flap cycle (check if data point is defined as new segemnt) 
        if round(steps(i==round(steps)))~=steps(i==round(steps)) %if the step which defines new segemnt is a fraction between two data points
            k=[(i-1),i]; %take data from both points
        else %if step is round number
            k=i; %take data from this step along
        end
        G(i).trajectory=plot([mean(data(4).coords(1,k)), mean(data(5).coords(1,k))],[mean(data(4).coords(3,k)), mean(data(5).coords(3,k))],'Color',[0.5 0.5 0.5],'LineWidth',1); %plot body location.  
        plot(mean(data(4).coords(1,k)),mean(data(4).coords(3,k)),'d','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5]) %plot head location as full diamond. 
%         plot([mean(data(2).coords(1,k)), mean(data(3).coords(1,k))]+mean(chord_shift(1,k)),[mean(data(2).coords(3,k)), mean(data(3).coords(3,k))]+mean(chord_shift(3,k)),'Color','k','LineWidth',2) %plot chord orientation over wingtip trajectory.  
        plot([mean(data(2).coords(1,k)), mean(data(3).coords(1,k))],[mean(data(2).coords(3,k)), mean(data(3).coords(3,k))],'Color','k','LineWidth',2) %plot chord location.  
%         plot(mean(data(2).coords(1,k))+mean(chord_shift(1,k)),mean(data(2).coords(3,k))+mean(chord_shift(3,k)),'o','Color','k','MarkerFaceColor','k') %plot leading edge over chord orientation as full cycle. 
        plot(mean(data(2).coords(1,k)),mean(data(2).coords(3,k)),'o','Color','k','MarkerFaceColor','k') %plot leading edge location as full cycle. 
%         G(i).velocity_vector=quiver(mean(data(1).coords(1,k)),mean(data(1).coords(3,k)),mean(V_mid_chord(k,1)),mean(V_mid_chord(k,3)),'color',[0 0.8 0],'LineWidth',2,'autoscale','off','markersize',3,'maxheadsize',3); %plot wing velocity vector
        step_num=find(i==round(steps))-1; %current step index number ('steps' variable contain zero step, so it needs to be reduced by 1)
        if i==steps(end) %last segment display value which is mean of 1st segment and last
            Force(wing_cycle).Fx(step_num)=mean(CFD(wing_cycle).Fx([1:10,end-10:end]))/1000000; %calc Fx in mN
            Force(wing_cycle).Fz(step_num)=mean(CFD(wing_cycle).Fz([1:10,end-10:end]))/1000000; %calc Fz in mN
        elseif i~=1 %display force of segment, calc from mean data around the arrow base
            F_index=round(length(CFD(wing_cycle).Fx)*(step_num/steps_amount)); %find index of F that correlates to the current segment
            Force(wing_cycle).Fx(step_num)=mean(CFD(wing_cycle).Fx(F_index-10:F_index+10))/1000000; %calc Fx in mN
            Force(wing_cycle).Fz(step_num)=mean(CFD(wing_cycle).Fz(F_index-10:F_index+10))/1000000; %calc Fz in mN
        end
        if step_num %skip 1st frame
%             G(i).force_vector=quiver(mean(data(1).coords(1,k)),mean(data(1).coords(3,k)),Force(wing_cycle).Fx(step_num),Force(wing_cycle).Fz(step_num),'color',[0.6 0.8 0],'LineWidth',3,'autoscale','off','markersize',3,'maxheadsize',3,'linestyle','-'); %plot wing CFD's force vector, from wingtip
            G(i).force_vector=quiver(mean([mean(data(2).coords(1,k)), mean(data(3).coords(1,k))]),mean([mean(data(2).coords(3,k)), mean(data(3).coords(3,k))]),Force(wing_cycle).Fx(step_num),Force(wing_cycle).Fz(step_num),'color',[0.6 0.8 0],'LineWidth',3,'autoscale','off','markersize',3,'maxheadsize',3,'linestyle','-'); %plot wing CFD's force vector, from mid-chord
        end
    end
end
for i=1:1:numel(data(5).coords(1,:)) %for each phase of wing cycle
    uistack(G(i).trajectory, 'bottom') %move body to the bottom of the drawing
%     uistack(G(i).velocity_vector, 'top') %move wing velocity vector to the top of the drawing
    uistack(G(i).force_vector, 'top') %move wing force vector to the top of the drawing
end
%% add annotation
t=sprintf('Wing cycle %1.2f - Body''s horizontal speed: %1.2f m/s',coords(wing_cycle).cycle_name,Vxy(wing_cycle)); title(t,'fontsize',22)
% t=sprintf('Body''s horizontal speed: %1.2f m/s\n Viewpoint - azimuth:%3.0f, elevation:%3.0f degrees',Vxy(wing_cycle),-horizontal_angle,-vertical_angle); title(t)
axis equal
xlim([-0.0004 0.0016])
ylim([-9 11].*10^-4)
xlabel('Horizontal axis (m)', 'FontSize',22) %X axis
ylabel('Vertical axis (m)', 'FontSize',22) %Y axis
legend(G(i).force_vector,'Force (mN)')
set(gca,'FontSize',20)
% xlabel('Horizontal axis (X) (m)', 'FontSize',12) %X axis
% ylabel('Horizontal axis (Y) (m)', 'FontSize',12) %Y axis
% zlabel('Vertical axis (m)', 'FontSize',12) %Z axis
cbh = colorbar ; %Create Colorbar
cbh.Ticks = 0: 1/10: 1;  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.1:1); (1:-0.1:0); ones(1,11)]'; %create matrix of RGB colormap: Low-Cyan, High-Magneta
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(0, 1, 11) ; %create labels
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.1f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
% s=sprintf('Trial_%d.jpg',wing_cycle);
% title('')
% saveas(gcf,s) %save figure as image named after flapping cycle number
end