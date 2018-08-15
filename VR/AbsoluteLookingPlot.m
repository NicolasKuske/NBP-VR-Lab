
%% plot average looking time of correct decided houses vs incorrect decided houses for the absolute task in VR and Belt over all sessions 

%working on Laptop or in Office?  
WorkPlace = str2double(input('Enter 1 for Laptop or 2 for Office: ','s'));

%which Dataset?
DataSet = input('Enter A for Dataset A and B for Dataset B: ', 's');

OfficePath='/net/home/student/n/nkuske/Dropbox/NKuskePhDProjects/AlignmentStudy/';
LaptopPath='C:\Users\nkuske\Dropbox\NKuskePhDProjects\AlignmentStudy\';

if WorkPlace==2
%PC
sjnums= readtable(strcat(OfficePath,'AlignmentAnalysis/Familiarity/Data/Set', DataSet, '/measurementList.txt'));
else
sjnums= readtable(strcat(LaptopPath,'AlignmentAnalysis\Familiarity\Data\Set', DataSet, '\measurementList.txt'));
end

 
M1Num = [];
M2Num = [];
M3Num = [];
for s=1:length(sjnums.Var1)
    if ~isnan(sjnums.Var3(s))
        M1Num = [M1Num sjnums.Var1(s)];
        M2Num = [M2Num sjnums.Var2(s)];
        M3Num = [M3Num sjnums.Var3(s)];
    end
end
Measurements = [M1Num;M2Num;M3Num];

Number = length(M1Num);

%matrices with one vector for each sessions
%where all correct and wrong decisions are stored
%each decision entry is the looking time (house information is lost)
Good=struct(); Good.Session1=[]; Good.Session2=[]; Good.Session3=[]; 
Bad=struct(); Bad.Session1=[]; Bad.Session2=[]; Bad.Session3=[]; 


%read in participants of every session
for iii=1:3
%for every entry in MiiiNum read in NumViewsD file.
for n = 1:Number
    disp(n);
    
    
%for every house in this participants 
%Output.Absolute.Trial_3s or Output.Absolute.Trial_Inf... 
if WorkPlace==2
%PC
Subj_HouseTime=load(strcat(OfficePath,'AlignmentAnalysis/Familiarity/Data/Set', DataSet,...
    '/NumViewsD_VP_',num2str(Measurements(iii,n)), '.mat'));

Subj_Perf=load(strcat(OfficePath,'AlignmentVR_AllPCsData/tracking_experiments_AlignmentVR/',...
    'Results/AlignmentVR_SubjNo_', num2str(Measurements(iii,n)), '.mat'));

else
Subj_HouseTime=load(strcat(LaptopPath,'\AlignmentAnalysis\Familiarity\Data\Set', DataSet, '\NumViewsD_','VP_', num2str(Measurements(iii,n)), '.mat'));
Subj_Perf=load(strcat(LaptopPath,'\AlignmentVR_AllPCsData\tracking_experiments_AlignmentVR\Results\AlignmentVR_SubjNo_', num2str(Measurements(iii,n)), '.mat'));

end   


    %convert house numbers to numbers
    %delete last N(ot)H(ouse) row
Subj_HouseTime.NumViews(height(Subj_HouseTime.NumViews),:) = [];
Subj_HouseTime.NumViews.House=str2double(extractBefore(Subj_HouseTime.NumViews.House,'_'));
        
    for ii=1:36
        
        %extract house number and from house number extract looking time
        %if both are the same, add time to good otherwise add time to bad
        if strcmp(Subj_Perf.Output.Absolute.Trial_3s(ii).Correct,Subj_Perf.Output.Absolute.Trial_3s(ii).Decision)
        Good.(strcat('Session',num2str(iii)))=[Good.(strcat('Session',num2str(iii))) Subj_HouseTime.NumViews.occ(Subj_HouseTime.NumViews.House == Subj_Perf.Output.Absolute.Trial_3s(ii).House_Nr)];
        else
        Bad.(strcat('Session',num2str(iii)))=[Bad.(strcat('Session',num2str(iii))) Subj_HouseTime.NumViews.occ(Subj_HouseTime.NumViews.House == Subj_Perf.Output.Absolute.Trial_3s(ii).House_Nr)];
        end
        
        if strcmp(Subj_Perf.Output.Absolute.Trial_Inf(ii).Correct,Subj_Perf.Output.Absolute.Trial_Inf(ii).Decision)
        Good.(strcat('Session',num2str(iii)))=[Good.(strcat('Session',num2str(iii))) Subj_HouseTime.NumViews.occ(Subj_HouseTime.NumViews.House == Subj_Perf.Output.Absolute.Trial_Inf(ii).House_Nr)];
        else
        Bad.(strcat('Session',num2str(iii)))=[Bad.(strcat('Session',num2str(iii))) Subj_HouseTime.NumViews.occ(Subj_HouseTime.NumViews.House == Subj_Perf.Output.Absolute.Trial_Inf(ii).House_Nr)];
        end
        
    end

end
end


%plot mean and scatter
%% PLOT
fL=figure('Name','AbsoluteLooking','NumberTitle','off');
axL=axes('Parent', fL);  

%plot scattered datapoints first so that the bar plot overlays on them

%create a xG value for each good datapoint (different color for plotting)
xG=[]; xB=[];

%good first session
for i=1:length(Good.Session1)
    
    xG= [xG , 0.85];
    
end    

%bad first session
for i=1:length(Bad.Session1)
    xB= [xB , 1.15];
end 


for i=1:length(Good.Session2)
    xG= [xG , 1.85];
end    

for i=1:length(Bad.Session2)
    xB= [xB , 2.15];
end 

for i=1:length(Good.Session3)
    xG= [xG , 2.85];
end    

for i=1:length(Bad.Session3)
    xB= [xB , 3.15];
end 


yG = [Good.Session1, Good.Session2, Good.Session3];
sz=52;
Gcol=[.0 .4 .0];

pL=scatter(axL,xG,yG,sz, Gcol, 'filled', 'jitter','on', 'jitterAmount',0.05);
pL.MarkerEdgeAlpha = 0.2;
pL.MarkerFaceAlpha = 0.2;

hold on

yB = [Bad.Session1, Bad.Session2, Bad.Session3];
sz=52;
Bcol=[.5 .0 .0];

pL=scatter(axL,xB,yB,sz, Bcol, 'filled', 'jitter','on', 'jitterAmount',0.05);
pL.MarkerEdgeAlpha = 0.2;
pL.MarkerFaceAlpha = 0.2;


%% Boxplots
%create groups
grp = [zeros(1,length(Good.Session1))+.85,ones(1,length(Bad.Session1))+.15,...
    ones(1,length(Good.Session2))+.85, ones(1,length(Bad.Session2))+1.15,...
    ones(1,length(Good.Session3))+1.85, ones(1,length(Bad.Session3))+2.15];
%plot the boxes
bx=boxplot(axL,[Good.Session1, Bad.Session1, Good.Session2, Bad.Session2,...
    Good.Session3, Bad.Session3], grp,'Positions',[0.85 1.15 1.85 2.15 2.85 3.15],...
    'Colors','k', 'Widths',0.25, 'Symbol','');

%linewidth&style
set(bx,'linew',1.5);
set(findobj(axL,'LineStyle','--'),'LineStyle','-');

lines = axL.Children; % get handles to the lines in the HGGroup object
uw = findobj(lines, 'tag', 'Upper Whisker');           % get handle to "Upper Whisker" line
uav = findobj(lines, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw = findobj(lines, 'tag', 'Lower Whisker');           % get handle to "Lower Whisker" line
lav = findobj(lines, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line

uw(1).YData(1,2) = quantile(Bad.Session3,.975);
uw(2).YData(1,2) = quantile(Good.Session3,.975);
uw(3).YData(1,2) = quantile(Bad.Session2,.975);
uw(4).YData(1,2) = quantile(Good.Session2,.975);
uw(5).YData(1,2) = quantile(Bad.Session1,.975);
uw(6).YData(1,2) = quantile(Good.Session1,.975);

uav(1).YData(:) = quantile(Bad.Session3,.975);
uav(2).YData(:) = quantile(Good.Session3,.975);
uav(3).YData(:) = quantile(Bad.Session2,.975);
uav(4).YData(:) = quantile(Good.Session2,.975);
uav(5).YData(:) = quantile(Bad.Session1,.975);
uav(6).YData(:) = quantile(Good.Session1,.975);

lw(1).YData(1,1) = quantile(Bad.Session3,.025);
lw(2).YData(1,1) = quantile(Good.Session3,.025);
lw(3).YData(1,1) = quantile(Bad.Session2,.025);
lw(4).YData(1,1) = quantile(Good.Session2,.025);
lw(5).YData(1,1) = quantile(Bad.Session1,.025);
lw(6).YData(1,1) = quantile(Good.Session1,.025);

lav(1).YData(:) = quantile(Bad.Session3,.025);
lav(2).YData(:) = quantile(Good.Session3,.025);
lav(3).YData(:) = quantile(Bad.Session2,.025);
lav(4).YData(:) = quantile(Good.Session2,.025);
lav(5).YData(:) = quantile(Bad.Session1,.025);
lav(6).YData(:) = quantile(Good.Session1,.025);





hold off



axL.FontSize=14;
axL.TickLength=[0.02 0.02];
% xlim(axS,[0 36]);
% axS.XTick=0:18:36;
ylim(axL,[0.0 75]);
axL.YTick=0:30:60;
% axS.XGrid='on';
%axL.YGrid='on';
xticks(axL,[1 2 3]);

%Naming each of the bar groups
xticklabels(axL,{ 'First', 'Second', 'Third'});

xlabel(axL,'Session')
ylabel(axL,'Fixation Time')

%title(axM,'All Subjects: Marginal Sessions ');

box(axL, 'off')

%create a legend object for the figure axis
lL=legend(axL, {'Correct Houses' 'Wrong Houses'});



%% zoomed plot
fLZ=figure('Name','AbsoluteZoomed','NumberTitle','off');
axLZ=axes('Parent', fLZ);  

pL=scatter(axLZ,xG,yG,sz, Gcol, 'filled', 'jitter','on', 'jitterAmount',0.05);
pL.MarkerEdgeAlpha = 0.15;
pL.MarkerFaceAlpha = 0.15;

hold on

pL=scatter(axLZ,xB,yB,sz, Bcol, 'filled', 'jitter','on', 'jitterAmount',0.05);
pL.MarkerEdgeAlpha = 0.15;
pL.MarkerFaceAlpha = 0.15;


%plot the boxes
bx=boxplot(axLZ,[Good.Session1, Bad.Session1, Good.Session2, Bad.Session2,...
    Good.Session3, Bad.Session3], grp,'Positions',[0.85 1.15 1.85 2.15 2.85 3.15],...
    'Colors','k', 'Widths',0.25, 'Symbol','');

%linewidth
set(bx,'linew',1.5);
set(findobj(axLZ,'LineStyle','--'),'LineStyle','-');
%here same function as above but generally more specific: set(findobj(gcf,'-regexp','Tag','\w*Whisker'),'LineStyle','-')

lines = axLZ.Children; % get handles to the lines in the HGGroup object
uw = findobj(lines, 'tag', 'Upper Whisker');           % get handle to "Upper Whisker" line
uav = findobj(lines, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
lw = findobj(lines, 'tag', 'Lower Whisker');           % get handle to "Lower Whisker" line
lav = findobj(lines, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line

uw(1).YData(1,2) = quantile(Bad.Session3,.975);
uw(2).YData(1,2) = quantile(Good.Session3,.975);
uw(3).YData(1,2) = quantile(Bad.Session2,.975);
uw(4).YData(1,2) = quantile(Good.Session2,.975);
uw(5).YData(1,2) = quantile(Bad.Session1,.975);
uw(6).YData(1,2) = quantile(Good.Session1,.975);

uav(1).YData(:) = quantile(Bad.Session3,.975);
uav(2).YData(:) = quantile(Good.Session3,.975);
uav(3).YData(:) = quantile(Bad.Session2,.975);
uav(4).YData(:) = quantile(Good.Session2,.975);
uav(5).YData(:) = quantile(Bad.Session1,.975);
uav(6).YData(:) = quantile(Good.Session1,.975);

lw(1).YData(1,1) = quantile(Bad.Session3,.025);
lw(2).YData(1,1) = quantile(Good.Session3,.025);
lw(3).YData(1,1) = quantile(Bad.Session2,.025);
lw(4).YData(1,1) = quantile(Good.Session2,.025);
lw(5).YData(1,1) = quantile(Bad.Session1,.025);
lw(6).YData(1,1) = quantile(Good.Session1,.025);

lav(1).YData(:) = quantile(Bad.Session3,.025);
lav(2).YData(:) = quantile(Good.Session3,.025);
lav(3).YData(:) = quantile(Bad.Session2,.025);
lav(4).YData(:) = quantile(Good.Session2,.025);
lav(5).YData(:) = quantile(Bad.Session1,.025);
lav(6).YData(:) = quantile(Good.Session1,.025);




hold off
axLZ.FontSize=14;
axLZ.TickLength=[0.02 0.02];
% xlim(axS,[0 36]);
% axS.XTick=0:18:36;
ylim(axLZ,[0.0 15]);
axLZ.YTick=0:5:15;
% axS.XGrid='on';
%axL.YGrid='on';
xticks(axLZ,[1 2 3]);

%Naming each of the bar groups
xticklabels(axLZ,{ 'First', 'Second', 'Third'});

xlabel(axLZ,'Session')
ylabel(axLZ,'Fixation Time')

%title(axM,'All Subjects: Marginal Sessions ');

box(axLZ, 'off')

%create a legend object for the figure axis
lLZ=legend(axLZ, {'Correct Houses' 'Wrong Houses'});


