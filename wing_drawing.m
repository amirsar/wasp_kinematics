load('D:\My Documents\My Desktop\Kine 2010\wur_flo_insect_tracker\B_manualtracking_kine\01A_kine_SA-3\models\wing_eret1.mat'); %load wing's outline model
%% plot wing's components
figure
plot3(coords(1,:),coords(2,:),zeros(41),'linewidth',2) %draw wing outline
axis equal
hold on
plot3([coords(1,1),coords(1,13)],[coords(2,1),coords(2,13)],[0 0],'linewidth',2) %leading vector
plot3([coords(1,1),coords(1,21)],[coords(2,1),coords(2,21)],[0 0],'linewidth',2) %wing's length
plot3([coords(1,1),coords(1,31)],[coords(2,1),coords(2,31)],[0 0],'linewidth',2) %trailing vector
plot3([coords(1,13),coords(1,31)],[coords(2,13),coords(2,31)],[0 0],'linewidth',2) %wing's span
%% calc perpendiculars
coords(3,1)=0; %add z dimension to the 2D wing outline
N1=cross(coords(:,31),coords(:,13)); %dorsal perpendicular to the wing
N1=N1./norm(N1); %unit vector of the above
quiver3(0,0,0,N1(1),N1(2),N1(3),'linewidth',2,'autoscale','off','maxheadsize',0.3) %wing's dorsal normal 
span=coords(:,13)-coords(:,31); %span vector
V=[-0.4,0.6,-0.2];  %imagenary velocity vector
spanV=V-(span./2)'; %velocity vector from middle of span
N2=cross(span,spanV); %dorsal perpendicular to the wing and velocity
N2=N2./norm(N2); %unit vector of the above
quiver3((coords(1,13)+coords(1,31))/2,0,0,N2(1),N2(2),N2(3),'linewidth',2,'autoscale','off','maxheadsize',0.3) %wing's dorsal normal 
quiver3((coords(1,13)+coords(1,31))/2,0,0,V(1),V(2),V(3),'linewidth',2,'linestyle','-.','autoscale','off','maxheadsize',0.3,'color','k') %imagenary velocity vector
%% remove axes & plot costum axes
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])
set(gca,'visible','off')
quiver3(-0.2,-0.2,-0.2,0.3,0,0,'linewidth',4,'color','k','autoscale','off','maxheadsize',0.3) %X axis
quiver3(-0.2,-0.2,-0.2,0,0.3,0,'linewidth',4,'color','k','autoscale','off','maxheadsize',0.3) %Y axis
quiver3(-0.2,-0.2,-0.2,0,0,0.3,'linewidth',4,'color','k','autoscale','off','maxheadsize',0.3) %Z axis