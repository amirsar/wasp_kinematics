%===========================================================================================%
% By Amir Sarig
% - Recieve number of wing cycles (lines) and their names (cycles_names)
% - Read angles from file 'Elevation3.xls' that should be in the active directory
% - Plot angles with mean
% - Store mean cycle
%
% Called from 'AnglesExcel_Amir' 
%
% Updates:
%27.02.2022
%Add constant Y limits
%21.07.2021
%Changed the output from 'count' to 'mean_cycle'
%20.04.2021
%change the spline calculation from 'Elv(:,c) = spline(t,A,tcorr);' to the
%current calculation that specify first and last endpoints of the plot.
%'tcorr' is now calculated once. it was calculated in each loop interation in all 3 loops
%07.01.2020
%changed the max length of wing cycle from 26 or 27 frames (Z or AA in excel) to 52
%frames (AZ in excel)
%24.07.2019 
%adjustments over Gal's version which was made for Eyal's data at 2018 
%===========================================================================================%

function [mean_cycle] = Meanangles(lines,cycles_names,tcorr,angles_storage)
%clear;
%close all

%count=input('Total amount of wings cycles?   ')*2; %ask user for amount of wing cycle that will be in the all the analysis, double it because each cycle has right and left wings
%file_path=uigetdir('','Choose the folder containing excel data file:'); %ask the user to enter the path of the data.
%cd(file_path)

count=lines; %amount of rows in the spreadsheet - number of wing cycles times two (1st for right wing, 2nd for left wing)


%assuming 1 empty line in the excel sheet between elevation to flapping and between flapping to incidence


%storing graphical data for ploting
which_wing='--'; %left wing
figure
hold on
for ind=1:(count/2)
    colors(ind).bar=get(plot(ind:100),'color');
end
close
figure
lgnd_cell=struct2cell(cycles_names);
% elevation
for c=1:count
%     range = strcat('A',num2str(c),':AZ',num2str(c));
%     A = radtodeg(xlsread('Elevation3.xls',1,range)); %read data and convert from radians to degrees
    A = radtodeg(angles_storage(c).elevation); %read data and convert from radians to degrees
    [r,rr] = size(A);
    t = (1./rr:1./rr:1.0)'; %X values (breakpoints) with data (Y value) in the wing cycle
    cs=spline(t,[A(1) A A(end)]); % piecewise polynomial form of the cubic spline interpolant of A for later use with the first and last value in A are used as the endslopes for the cubic spline.
    Elv(:,c) = ppval(cs,tcorr); %Spline Interpolation with Specified Endpoint Slopes
    subplot(1,4,2)
    temp=plot(tcorr,Elv(:,c),which_wing);
    set(temp,'color',colors(round(c/2)).bar)
    lgnd(c).data=temp;
    if strcmp('-.',which_wing)
        which_wing='--'; %left wing
    else
        which_wing='-.'; %right wing
    end
    hold on
end
title('Elevation');
Elv(:,c+1) = mean(Elv(:,1:count)');
plot(tcorr,Elv(:,c+1), 'k', 'LineWidth', 2.5)
xlabel('Flap cycle fraction') %X axis label for middle figure in the plot implay for all the row
ylim([-135 135])

% flapping / horizontal
subplot(1,4,1)
for c=(count+2):(2*count+1)
%     range = strcat('A',num2str(c),':AZ',num2str(c));
%     A = radtodeg(xlsread('Elevation3.xls',1,range)); %read data and convert from radians to degrees
    A = radtodeg(angles_storage(c-(count+1)).horizontal); %read data and convert from radians to degrees
    [r,rr] = size(A);   
    t = (1./rr:1./rr:1.0)'; %X values (breakpoints) with data (Y value) in the wing cycle
    cs=spline(t,[A(1) A A(end)]); % piecewise polynomial form of the cubic spline interpolant of A for later use with the first and last value in A are used as the endslopes for the cubic spline.
    Flp(:,c-(count+1)) = ppval(cs,tcorr); %Spline Interpolation with Specified Endpoint Slopes
    temp=plot(tcorr,Flp(:,c-(count+1)),which_wing);
    set(temp,'color',colors(round((c-(count+2)+1)/2)).bar)
    if strcmp('-.',which_wing)
        which_wing='--';
    else
        which_wing='-.';
    end
    hold on
end
title('Horizontal'); %Flapping
Flp(:,c-count) = mean(Flp(:,1:count)');
plot(tcorr,Flp(:,c-count), 'k', 'LineWidth', 2.5)
ylabel('Angle (degree)') %Y axis label for 1st figure in the plot implay for all the row
ylim([-135 135])

% incidence
subplot(1,4,3)
for c=(2*count+3):(3*count+2)
%     range = strcat('A',num2str(c),':AZ',num2str(c));
%     A = radtodeg(xlsread('Elevation3.xls',1,range)); %read data and convert from radians to degrees
    A = radtodeg(angles_storage(c-(2*count+2)).incidence); %read data and convert from radians to degrees
    [r,rr] = size(A);
    t = (1./rr:1./rr:1.0)'; %X values (breakpoints) with data (Y value) in the wing cycle
    cs=spline(t,[A(1) A A(end)]); % piecewise polynomial form of the cubic spline interpolant of A for later use with the first and last value in A are used as the endslopes for the cubic spline.
    Inc(:,c-(2*count+2)) = ppval(cs,tcorr); %Spline Interpolation with Specified Endpoint Slopes
    temp=plot(tcorr,Inc(:,c-(2*count+2)),which_wing);
    set(temp,'color',colors(round((c-(2*count+3)+1)/2)).bar)
    if strcmp('-.',which_wing)
        which_wing='--';
    else
        which_wing='-.';
    end
    hold on
end
title('Incidence');Inc(:,c-(2*count+1)) = mean(Inc(:,1:count)');
plot(tcorr,Inc(:,c-(2*count+1)), 'k', 'LineWidth', 2.5)
new_lgnd=lgnd(1:2:end);
leg=legend([new_lgnd(:).data],lgnd_cell(1,1,:));
ylim([-180 180])

%%Store data
mean_cycle.Elv=Elv(:,end)';
mean_cycle.Flp=Flp(:,end)';
mean_cycle.Inc=Inc(:,end)';
%Super title, work only in Matlab 2018b +
%sgtitle({'Dash-dot: Right wing'; 'Dashed : Left wing'}); 

end