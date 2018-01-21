%--------------------------Overall SubjectMapAnalysis----------------------
%this script needs output from the script "PositionAnalysis" which analyzed
%the individual subjects
PartList={3755,6876};
savepath = 'C:/Users/vivia/Dropbox/Project Seahaven/Tracking/Position/Results/';
%--------------------------------------------------------------------------
%% Show overlaid maps
Number = length(PartList);
map = imread('map5.png'); 
map = imresize(map,[500 450]);
for ii = 1: Number
    e = cell2mat(PartList(ii));
    x = load(['path_VP_' num2str(e) '.mat']);
    color = randi([0 255],1,3);
    len = size(x.path,2);
    for a=1:len-1
        map(int64(x.path(1,a)),int64(x.path(2,a)),1) = color(1);
        map(int64(x.path(1,a)),int64(x.path(2,a)),2) = color(2);
        map(int64(x.path(1,a)),int64(x.path(2,a)),3) = color(3);
    end
end
image(map)%save as jpg when displayed
imwrite(map,fullfile(savepath,'OverlaidMaps.jpeg'));
%% Show heatmap
newp=zeros([51 46]);
for x=1:500
    for y =1:450
        newp(floor(x/10)+1,floor(y/10)+1) = newp(floor(x/10)+1,floor(y/10)+1)+pos(x,y);
    end
end
h=pcolor(newp/norm(newp));colorbar;
set(h, 'EdgeColor', 'none');
%% show single maps
for ii = 1: Number
    e = cell2mat(PartList(ii));
    x = load(['map_VP_' num2str(e) '.mat']);
    figure;
    imshow(x.map);
    hold on;
end
%% Individual North
%shows graph of individual north for all subjects
% true north = rotation of 270
r = 1; % Radius
for ii = 1: Number
    e = cell2mat(PartList(ii));
    n = load(['North_VP_' num2str(e) '.mat']);
    t = cell2mat(n.north(3))-180; % Angle in degrees, -180 to have north on top
    [x,y] = pol2cart(t/180*pi,r);
    hold on;
    plot([0 x],[0,y])
    legendInfo{ii} = ['Subject = ' num2str(e)];
end
t = 90;%true north at 270 degrees -> -180 = 90
[x,y] = pol2cart(t/180*pi,r);
hold on;
plot([0 x],[0,y])
legendInfo{Number+1} = ['True North'];
legend(legendInfo)
saveas(gcf,fullfile(savepath,'IndividualNorth.jpeg'));
%clear all;