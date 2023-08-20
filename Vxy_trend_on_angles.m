%===========================================================================================%
% By Amir Sarig 06/03/2022
% - Recieve number of wing cycles (lines) and their names (cycles_names)
% - Plot R & L wings mean angles of all input cycles, while color coded by Vxy
%
% Called from 'extract_data_wasp_forward_flight' 
%===========================================================================================%

function [Flp, Elv, Inc, angles4sim] = Vxy_trend_on_angles(Vxy,angles_storage,trials_names,cycles_names,interpolant_flap_cycle_fraction)
relative_Vxy=Vxy-min(Vxy); %normalize body horizonta velocity so that minimal value is zero
relative_Vxy=relative_Vxy./max(relative_Vxy); %normalize body horizonta velocity so that maximal value is one
side_name={'Left','Right'};
% figure;
% hold on
% flapping / horizontal
% subplot(1,3,1)
horizontal=figure;
hold on
for wing_cycle=1:2:length(angles_storage)
    for side=0:1
        angles4sim(wing_cycle+side).cycle_name=angles_storage(wing_cycle).cycle_name; %record wing cycle name
        angles4sim(wing_cycle+side).side=side_name(side+1); %right/left
        ind=round(wing_cycle/2);
%         data = radtodeg(mean([angles_storage(wing_cycle).horizontal; angles_storage(wing_cycle+1).horizontal],1)); %read data, calc mean of R & L wings, and convert from radians to degrees
        data = radtodeg(angles_storage(wing_cycle+side).horizontal); %read data and convert from radians to degrees
        rr = size(data,2);
        flap_cycle_fraction = (1./rr:1./rr:1.0)'; %X values (breakpoints) of data (Y value) in the wing cycle
        cs=spline(flap_cycle_fraction,[data(1) data data(end)]); % piecewise polynomial form of the cubic spline interpolant of data with the first and last value in data are used as the endslopes for the cubic spline
        angles4sim(wing_cycle+side).horizontal=ppval(cs,interpolant_flap_cycle_fraction); %Spline Interpolation with Specified Endpoint Slopes
    end
        Flp(:,ind) = mean([angles4sim(wing_cycle).horizontal,angles4sim(wing_cycle+1).horizontal],2);
        plot(interpolant_flap_cycle_fraction,Flp(:,ind),'color',[relative_Vxy(ind) 0 1-relative_Vxy(ind)],'LineWidth',2);
        hold on
end
title('Horizontal'); %Flapping
xlabel('Flap cycle fraction') %X axis label 
ylabel('Angle (degree)') %Y axis label for 1st figure in the plot implay for all the row
ylim([-60 60])
% elevation
% subplot(1,3,2)
elevation=figure;
hold on
for wing_cycle=1:2:length(angles_storage)
    for side=0:1
        ind=round(wing_cycle/2);
%         data = radtodeg(mean([angles_storage(wing_cycle).elevation; angles_storage(wing_cycle+1).elevation],1)); %read data, calc mean of R & L wings, and convert from radians to degrees
        data = radtodeg(angles_storage(wing_cycle+side).elevation); %read data and convert from radians to degrees
        rr = size(data,2);
        flap_cycle_fraction = (1./rr:1./rr:1.0)'; %X values (breakpoints) of data (Y value) in the wing cycle
        cs=spline(flap_cycle_fraction,[data(1) data data(end)]); % piecewise polynomial form of the cubic spline interpolant of data for later use with the first and last value in data are used as the endslopes for the cubic spline. Use tcorr from above
        angles4sim(wing_cycle+side).elevation=ppval(cs,interpolant_flap_cycle_fraction); %Spline Interpolation with Specified Endpoint Slopes
    end
    Elv(:,ind) = mean([angles4sim(wing_cycle).elevation,angles4sim(wing_cycle+1).elevation],2);
    plot(interpolant_flap_cycle_fraction,Elv(:,ind),'color',[relative_Vxy(ind) 0 1-relative_Vxy(ind)],'LineWidth',2);
    hold on
end
title('Elevation');
xlabel('Flap cycle fraction') %X axis label 
ylabel('Angle (degree)') %Y axis label for 1st figure in the plot implay for all the row
ylim([-100 100])
% incidence
% subplot(1,3,3)
incidence=figure;
hold on
for wing_cycle=1:2:length(angles_storage)
    for side=0:1
        ind=round(wing_cycle/2);
%         data = radtodeg(mean([angles_storage(wing_cycle).incidence; angles_storage(wing_cycle+1).incidence],1)); %read data, calc mean of R & L wings, and convert from radians to degrees
        data = radtodeg(angles_storage(wing_cycle+side).incidence); %read data and convert from radians to degrees
        rr = size(data,2);
        flap_cycle_fraction = (1./rr:1./rr:1.0)'; %X values (breakpoints) of data (Y value) in the wing cycle
        cs=spline(flap_cycle_fraction,[data(1) data data(end)]); % piecewise polynomial form of the cubic spline interpolant of data with the first and last value in data are used as the endslopes for the cubic spline. Use tcorr from above
        angles4sim(wing_cycle+side).incidence=ppval(cs,interpolant_flap_cycle_fraction); %Spline Interpolation with Specified Endpoint Slopes
    end
    Inc(:,ind) = mean([angles4sim(wing_cycle).incidence,angles4sim(wing_cycle+1).incidence],2);
    temp=plot(interpolant_flap_cycle_fraction,Inc(:,ind),'color',[relative_Vxy(ind) 0 1-relative_Vxy(ind)],'LineWidth',2);
    hold on
    lgnd(ind).data=temp;
%     lgnd(ind).name=sprintf('%s.%d',trials_names(ind).name(end-4:end-3),cycles_names(ind)); %record cycle name
    lgnd(ind).name=sprintf('%1.3f',Vxy(ind)); %record cycle Vxy
end
title('Incidence');
xlabel('Flap cycle fraction') %X axis label 
ylabel('Angle (degree)') %Y axis label for 1st figure in the plot implay for all the row
ylim([-145 145])
% legend
figure
cbh =colorbar;
% pos = get(subplot(1,3,3),'Position');
% cbh =colorbar('Position', [pos(1)+pos(3)+0.02  pos(2)  0.02  pos(4)]); %Create Colorbar
cbh.Ticks = [0: 1/10: 1];  %Create & update 11 ticks from minimal wingtip velocity to maximal
c=[(0:0.01:1); zeros(1,101); (1:-0.01:0)]'; %create matrix of RGB colormap: Low-Blue, High-Red
colormap(c); %apply colormap that fit to the plot range
TickLabels = linspace(min(Vxy), max(Vxy), 11) ; %create labels
for i=1:11
    TickLabels_string(i,:) = num2str(TickLabels(i),'%1.3f'); %modify labels values to string
end
cbh.TickLabels=TickLabels_string; %update labels in figure
% leg=legend([lgnd(p(:)).data],lgnd(p(:)).name); %display legend while cycles are sorted by Vxy
% title(leg,'Vxy (m/s)'); %works only with Matlab 2017+
end