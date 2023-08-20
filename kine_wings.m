%07.01.2020
%changed 'frames' to be function input rather than user manually input

function [angles] = kine_wings(data, frames)

%% To begin, open the solution of the kine (the variable data should appear)

%clearvars -except data;
% close all;


% Obtaining body & wing points
frames = frames(1):frames(2);
%frames = input('start frame?  '):input('end frame?  ');
head = data.kine.body.data.coords(1,:,frames); head = reshape(head,3,size(head,3)); %reshape so each column is a single dot
abd = data.kine.body.data.coords(2,:,frames); abd = reshape(abd,3,size(abd,3));
Lbase = data.kine.left_wing.data.coords(1,:,frames); Lbase = reshape(Lbase,3,size(Lbase,3));
Ltip = data.kine.left_wing.data.coords(2,:,frames); Ltip = reshape(Ltip,3,size(Ltip,3));
Rbase = data.kine.right_wing.data.coords(1,:,frames); Rbase = reshape(Rbase,3,size(Rbase,3));
Rtip = data.kine.right_wing.data.coords(2,:,frames); Rtip = reshape(Rtip,3,size(Rtip,3));
Llead = data.kine.left_wing.data.model_coords(:,13,frames); Llead = reshape(Llead,3,size(Llead,3));  %third point on the wing's leading edge
Rlead = data.kine.right_wing.data.model_coords(:,13,frames); Rlead = reshape(Rlead,3,size(Rlead,3));
Ltrail = data.kine.left_wing.data.model_coords(:,31,frames); Ltrail = reshape(Ltrail,3,size(Ltrail,3));  %third point on the wing's trailing edge
Rtrail = data.kine.right_wing.data.model_coords(:,31,frames); Rtrail = reshape(Rtrail,3,size(Rtrail,3));

% Defining the axes of the body

x = (head-abd)./repmat(sqrt((head(1,:)-abd(1,:)).^2+(head(2,:)-abd(2,:)).^2+(head(3,:)-abd(3,:)).^2),3,1); %body length unit vector
yy = (Lbase-Rbase)./repmat(sqrt((Lbase(1,:)-Rbase(1,:)).^2+(Lbase(2,:)-Rbase(2,:)).^2+(Lbase(3,:)-Rbase(3,:)).^2),3,1); %body width unit vector
tet = acos(dot(yy,x)); %angle between x & yy (both are already unit vectors)
zz = cross(x,yy); %body hight vector
z = zz./repmat(sin(tet),3,1);  %body hight unit vector (xXyy=1*1*sin(tet))
y = cross(z,x);  %corrected body width unit vector (so it is orthogonal to x)
TM = permute(reshape([x;y;z],3,3,length(frames)),[2,1,3]);  %rotation matrix for each frame(each 2D matrix is [x1,x2,x3;y1,y2,y3;z1,z2,z3])
AllD = reshape([head;abd;Lbase;Ltip;Rbase;Rtip;Llead;Rlead;Ltrail;Rtrail],3,10,length(frames));  %All Data
AllD = AllD-repmat((AllD(:,1,:)+AllD(:,2,:))/2,1,10);

% Rotation

P = zeros(size(AllD));
for n = 1:length(frames) %for each frame
    P(:,:,n) = TM(:,:,n)*AllD(:,:,n); %rotate all digitizing points to body cordinate system
end
% NormP = P./repmat(P(1,1,:)-P(1,2,:),3,10);
% LwingTraj(:,:) = NormP(:,4,:);
% RwingTraj(:,:) = NormP(:,6,:);

% New axes of the body

Lwing = reshape(P(:,4,:)-P(:,3,:),3,length(frames)); %left wing length (tip-base) vector in body coordinate system, in all frames
Rwing = reshape(P(:,6,:)-P(:,5,:),3,length(frames)); %right wing length (tip-base) vector in body coordinate system, in all frames

% New calculations for the old angles. Used in figure 6 of the paper

Lflap = atan2(Lwing(2,:),Lwing(1,:));
Rflap = atan2(-Rwing(2,:),Rwing(1,:));
Lelev = atan2(Lwing(3,:),sqrt(Lwing(1,:).^2+Lwing(2,:).^2));
Relev = atan2(Rwing(3,:),sqrt(Rwing(1,:).^2+Rwing(2,:).^2));

% Collecting data for figure 5b

LeftTip = [reshape(P(1,4,:),1,length(frames))./(2*reshape(P(1,1,:),1,length(frames)));reshape(P(3,4,:),1,length(frames))./(2*reshape(P(1,1,:),1,length(frames)))]; %x,z normalized by body length
RightTip = [reshape(P(1,6,:),1,length(frames))./(2*reshape(P(1,1,:),1,length(frames)));reshape(P(3,6,:),1,length(frames))./(2*reshape(P(1,1,:),1,length(frames)))]; %x,z normalized by body length
Lchord = reshape(P(:,7,:)-P(:,9,:),3,length(frames)); %wing chord vectors
Rchord = reshape(P(:,8,:)-P(:,10,:),3,length(frames)); %wing chord vectors
Lleng = sqrt(Lchord(1,:).^2+Lchord(3,:).^2)./(2*reshape(P(1,1,:),1,length(frames))); %normalized by body length
Rleng = sqrt(Rchord(1,:).^2+Rchord(3,:).^2)./(2*reshape(P(1,1,:),1,length(frames))); %normalized by body length
L5ang = atan2(Lchord(3,:),Lchord(1,:));
R5ang = atan2(Rchord(3,:),Rchord(1,:));

% Calculating the incidence angle: Rotating to xy, then to y

% LEFT - wing to positive y
Lyz = atan2(Lwing(3,:),Lwing(2,:)); %elevation only in the body yz plane, not "true" elevation
PLeft = P-repmat(P(:,3,:),[1,10,1]); %all rotated digitizing points minus rotated left wing base
Lrot1 = zeros(3,4,length(frames));
for n = 1:length(frames) %for each frame
    Lrot1(:,:,n) = [1,0,0;0,cos(Lyz(n)),sin(Lyz(n));0,-sin(Lyz(n)),cos(Lyz(n))]*PLeft(:,[3,4,7,9],n); %rotation of left wing's base,tip,lead,trail around bod's X axis to body xy plane 
end
Lwing2 = reshape(Lrot1(:,2,:)-Lrot1(:,1,:),3,length(frames));  %left wing span vectors of all frames
Lxy = atan2(Lwing2(1,:),Lwing2(2,:)); %angle of the wing after rotation to xy plane, relative to the y axis (positive y is 0, head is +90)
Lrot2 = zeros(3,4,length(frames));
for n = 1:length(frames) %for each frame
    Lrot2(:,:,n) = [cos(Lxy(n)),-sin(Lxy(n)),0;sin(Lxy(n)),cos(Lxy(n)),0;0,0,1]*Lrot1(:,:,n); %rotation of left wing's base,tip,lead,trail around body's Z axis to body positive y
end
Linci(1,:) = atan2((Lrot2(3,3,:)-Lrot2(3,4,:)),(Lrot2(1,3,:)-Lrot2(1,4,:))); %angle of wing chord (lead-trail) on the xz plane, after rotating the wing to the y axis (fling is negative)
% LeftAng(1,:) = atan2(Lrot2(3,3,:),Lrot2(1,3,:));
% LeftAng(2,:) = atan2(-Lrot2(3,4,:),-Lrot2(1,4,:));
% LeftAng(3,:) = atan2((Lrot2(3,3,:)-Lrot2(3,4,:)),(Lrot2(1,3,:)-Lrot2(1,4,:)));

% RIGHT - wing to negative y
Ryz = pi+atan2(Rwing(3,:),Rwing(2,:)); %elevation only in the yz plane, not "true" elevation
PRight = P-repmat(P(:,5,:),[1,10,1]);
Rrot1 = zeros(3,4,length(frames));
for n = 1:length(frames) %rotation of wing to xy plane
    Rrot1(:,:,n) = [1,0,0;0,cos(Ryz(n)),sin(Ryz(n));0,-sin(Ryz(n)),cos(Ryz(n))]*PRight(:,[5,6,8,10],n); %base,tip,lead,trail
end
Rwing2 = reshape(Rrot1(:,2,:)-Rrot1(:,1,:),3,length(frames));  %wing length vectors
Rxy = pi+atan2(Rwing2(1,:),Rwing2(2,:)); %angle of the wing after rotation to xy plane, relative to the y axis (positive y is 0, head is +90)
Rrot2 = zeros(3,4,length(frames));
for n = 1:length(frames) %rotation of to negative y
    Rrot2(:,:,n) = [cos(Rxy(n)),-sin(Rxy(n)),0;sin(Rxy(n)),cos(Rxy(n)),0;0,0,1]*Rrot1(:,:,n); %base,tip,lead,trail
end
Rinci(1,:) = atan2((Rrot2(3,3,:)-Rrot2(3,4,:)),(Rrot2(1,3,:)-Rrot2(1,4,:))); %angle of wing chord (lead-trail) on the xz plane, after rotating the wing to the y axis (fling is negative)

% NEW - save data for Gal (for moment calculation)
Lflap(Lflap<-pi/2) = Lflap(Lflap<-pi/2)+pi*2;
Rflap(Rflap<-pi/2) = Rflap(Rflap<-pi/2)+pi*2;
Ryz = atan2(Rwing(3,:),-Rwing(2,:));
Rxy = atan2(Rwing2(1,:),-Rwing2(2,:));
reshapeMat = reshape(1:numel(AllD),30 ,size(AllD,3))'; %rearrange so rows are frames and columns are xyz
reshapeMat = reshapeMat(:,[1:12,19:21,25:27,13:18,22:24,28:30]); %rearrange to the order of dataOrder
rawData = AllD(reshapeMat);
rotatedData = P(reshapeMat);
oldAngles = [Lflap;Lyz;Linci;Rflap;Ryz;Rinci]';
angles = [-Lxy;Lyz;-Linci;-Rxy;Ryz;-Rinci]'; %(-) changes xy angles so head is -90, and inci angles so fling is positive
dataOrder = {'Head','Abdomen','Left base','Left tip','Left leading','Left trailing','Right base','Right tip','Right leading','Right trailing'};
oldAngleOrder = {'Left sweep','Left elevation','Left incidence','Right sweep','Right elevation','Right incidence'};
angleOrder = {'Left Horizontal','Left elevation','Left incidence','Right Horizontal','Right elevation','Right incidence'};

% figure(1); plot(rad2deg((Lflap+Rflap)/2)); figure(2); plot(rad2deg((Lelev+Relev)/2)); figure(3); plot(rad2deg((Linci+Rinci)/2)); 

% save('D:\My Documents\My Desktop\Er11.mat','angleOrder','angles','dataOrder','oldAngleOrder','oldAngles','rawData','rotatedData')


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@old junk 

% RightAng(1,:) = atan2(Rrot2(3,3,:),Rrot2(1,3,:));
% RightAng(2,:) = atan2(-Rrot2(3,4,:),-Rrot2(1,4,:));
% RightAng(3,:) = atan2((Rrot2(3,3,:)-Rrot2(3,4,:)),(Rrot2(1,3,:)-Rrot2(1,4,:)));

% Test = sum((round(LeftAng(1,:)*1000)==round(LeftAng(2,:)*1000))+(round(LeftAng(2,:)*1000)==round(LeftAng(3,:)*1000))+(round(RightAng(1,:)*1000)==round(RightAng(2,:)*1000))+(round(RightAng(2,:)*1000)==round(RightAng(3,:)*1000)))==4*size(LeftAng,2)
% the test always succeeded

% % NO GOOD - rotating to yz, then to y
% 
% % LEFT - wing to positive y
% NoLxy = atan(Lwing(1,:)./Lwing(2,:))+pi*(Lwing(2,:)<0);
% NoPLeft = P-repmat(P(:,3,:),[1,10,1]);
% NoLrot1 = zeros(3,4,length(frames));
% for n = 1:length(frames)
%     NoLrot1(:,:,n) = [cos(NoLxy(n)),-sin(NoLxy(n)),0;sin(NoLxy(n)),cos(NoLxy(n)),0;0,0,1]*NoPLeft(:,[3,4,7,9],n); %base,tip,lead,trail
% end
% NoLwing2 = reshape(NoLrot1(:,2,:)-NoLrot1(:,1,:),3,length(frames));  %wing length vectors
% NoLyz = atan(NoLwing2(3,:)./NoLwing2(2,:))+pi*(NoLwing2(2,:)<0);
% NoLrot2 = zeros(3,4,length(frames));
% for n = 1:length(frames)
%     NoLrot2(:,:,n) = [1,0,0;0,cos(NoLyz(n)),sin(NoLyz(n));0,-sin(NoLyz(n)),cos(NoLyz(n))]*NoLrot1(:,:,n); %base,tip,lead,trail
% end
% NoLeftAng(1,:) = atan(NoLrot2(3,3,:)./NoLrot2(1,3,:))+pi*(NoLrot2(1,3,:)<0);
% NoLeftAng(2,:) = atan(NoLrot2(3,4,:)./NoLrot2(1,4,:))+pi*(NoLrot2(1,4,:)>0);
% NoLeftAng(3,:) = atan((NoLrot2(3,3,:)-NoLrot2(3,4,:))./(NoLrot2(1,3,:)-NoLrot2(1,4,:)))+pi*((NoLrot2(1,3,:)-NoLrot2(1,4,:))<0);
% 
% % RIGHT - wing to negative y
% NoRxy = atan(Rwing(1,:)./Rwing(2,:))+pi*(Rwing(2,:)>0);
% NoPRight = P-repmat(P(:,5,:),[1,10,1]);
% NoRrot1 = zeros(3,4,length(frames));
% for n = 1:length(frames)
%     NoRrot1(:,:,n) = [cos(NoRxy(n)),-sin(NoRxy(n)),0;sin(NoRxy(n)),cos(NoRxy(n)),0;0,0,1]*NoPRight(:,[5,6,8,10],n); %base,tip,lead,trail
% end
% NoRwing2 = reshape(NoRrot1(:,2,:)-NoRrot1(:,1,:),3,length(frames));  %wing length vectors
% NoRyz = atan(NoRwing2(3,:)./NoRwing2(2,:))+pi*(NoRwing2(2,:)>0);
% NoRrot2 = zeros(3,4,length(frames));
% for n = 1:length(frames)
%     NoRrot2(:,:,n) = [1,0,0;0,cos(NoRyz(n)),sin(NoRyz(n));0,-sin(NoRyz(n)),cos(NoRyz(n))]*NoRrot1(:,:,n); %base,tip,lead,trail
% end
% NoRightAng(1,:) = atan(NoRrot2(3,3,:)./NoRrot2(1,3,:))+pi*(NoRrot2(1,3,:)<0);
% NoRightAng(2,:) = atan(NoRrot2(3,4,:)./NoRrot2(1,4,:))+pi*(NoRrot2(1,4,:)>0);
% NoRightAng(3,:) = atan((NoRrot2(3,3,:)-NoRrot2(3,4,:))./(NoRrot2(1,3,:)-NoRrot2(1,4,:)))+pi*((NoRrot2(1,3,:)-NoRrot2(1,4,:))<0);
% NoTest = sum((round(NoLeftAng(1,:)*1000)==round(NoLeftAng(2,:)*1000))+(round(NoLeftAng(2,:)*1000)==round(NoLeftAng(3,:)*1000))+(round(NoRightAng(1,:)*1000)==round(NoRightAng(2,:)*1000))+(round(NoRightAng(2,:)*1000)==round(NoRightAng(3,:)*1000)))==4*size(NoLeftAng,2)
% 
% %Lelev = acos((Lwing(1,:)^2+Lwing(2,:)^2)./(Lwing(1,:)^2+Lwing(2,:)^2+Lwing(3,:)^2));
end