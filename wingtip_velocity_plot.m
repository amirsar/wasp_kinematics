function wingtip_velocity_plot(mean_L_R_wingtip,V_tip_storage,I,insect_corners_x,insect_corners_y)
%% exctract data
% if any(pts2plot==4) %if user chose to plot left wing
%     V_tip=V_tip_storage.left;
%     mean_x=normalize_wingtip_trajectory.left(1,:);
%     mean_z=normalize_wingtip_trajectory.left(3,:);
% elseif any(pts2plot==6) %if user chose to plot right wing
%     V_tip=V_tip_storage.right;
%     mean_x=normalize_wingtip_trajectory.right(1,:);
%     mean_z=normalize_wingtip_trajectory.right(3,:);
% end
V_tip.x=mean([V_tip_storage.body.left.x; V_tip_storage.body.right.x],1);
V_tip.y=mean([V_tip_storage.body.left.y; V_tip_storage.body.right.y],1);
V_tip.z=mean([V_tip_storage.body.left.z; V_tip_storage.body.right.z],1);
mean_x=mean_L_R_wingtip(1,:);
mean_z=mean_L_R_wingtip(3,:);
    
Vtip=((median(V_tip.x,1).^2)+(median(V_tip.y,1).^2)+(median(V_tip.z,1).^2)).^(1/2);%calc 3D velocity
Vtip_range=[min(Vtip) max(Vtip)]; %store minimal and maximal wingtip velocity 

Vtip=Vtip-min(Vtip); %normalize wingtip velocity so that minimal value is zero
Vtip=Vtip./max(Vtip); %normalize wingtip velocity so that maximal value is one

%% plot data
figure
for i=1:1:numel(Vtip)-1 %for each phase of wing cycle (exept the last)
    %alternative 2 color scale of blue-red [(0:0.01:1); repmat(0,11,1); (1:-0.01:0)]
    plot(mean_x(i:i+1),mean_z(i:i+1),'Color',[Vtip(i),sin((Vtip(i))*3.14),1-Vtip(i)],'LineWidth',4) %plot mean 2D trajectory betweem this phase to the next one. color it according to wingtip velocity
    hold on 
end
plot(mean_x(i+1:1),mean_z(i+1:1),'Color',[Vtip(i+1),sin((Vtip(i))*3.14),1-Vtip(i+1)],'LineWidth',4) %plot mean 2D trajectory betweem last phase to the first one. color it according to wingtip velocity

%% plot annotations
h = image(insect_corners_x,-insect_corners_y,I); %insert image to the figure
uistack(h,'bottom') %set flyer image as background
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0, 101, 11);  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.01:1); sin((0:0.01:1)*3.14); (1:-0.01:0)]'; %create matrix of RGB colormap: Low-Blue, Middle-Green, High-Red
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(Vtip_range(1), Vtip_range(2), 11) ; %create labels
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.3f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
xlabel('Longitudinal axis (body lengths)', 'FontSize',24)
ylabel('Dorsoventral axis (body lengths)', 'FontSize',24)
xlim([-1.5 1.5])
ylim([-2 2])
grid on
set(gca,'FontSize',20);
title('Wingtip aerial velocity (m/s)');
end