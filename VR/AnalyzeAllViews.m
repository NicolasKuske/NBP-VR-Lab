
%% ----------------Analyze Raw ViewedHouses Files (1st Level)----------------
% Author is Viviane Kakerbeck, slight modifications from Nicolas Kuske
%
% Transforms ViewedHouses a txt list giving the number of house looked at every time 
% into NumViewD mat table that holds the accumulated time (and average distance) each house was looked at  


%working on Laptop or in Office?  
WorkPlace = str2double(input('Enter 1 for Laptop or 2 for Office: ','s'));

%which Dataset?
DataSet = input('Enter A for Dataset A and B for Dataset B: ', 's');


if WorkPlace==2
%PC
savepath = strcat('/net/home/student/n/nkuske/Dropbox/NKuskePhDProjects/AlignmentStudy/AlignmentAnalysis/Familiarity/Data/Set', DataSet);
files = dir(strcat(savepath,'/','ViewedHouses_VP*.txt'));%Analyzes all subjectfiles in that directory
else
savepath = strcat('C:\Users\nkuske\Dropbox\NKuskePhDProjects\AlignmentStudy\AlignmentAnalysis\Familiarity\Data\Set', DataSet);
files = dir(strcat(savepath,'\','ViewedHouses_VP*.txt'));
end

%--------------------------------------------------------------------------

Number = length(files);
avgdist = cell(1,Number);
for ii = 5:5%Number
    
    
    if WorkPlace==2
%PC
data = fopen(strcat(files(ii).folder,'/', files(ii).name));
else
data = fopen(strcat(files(ii).folder,'\', files(ii).name));
    end
    
    suj_num = files(ii).name(16:19) %leave to see output for check

    data = textscan(data,'%s','delimiter', '\n');
    data = data{1};
    data = table2array(cell2table(data));
    %initialize fields
    len = int32(length(data));
    houses = cell(1,len);
    distance = zeros(1,len);
    timestamps = zeros(1,len);
    for a = 1:double(len)
        line = textscan(data{a},'%s','delimiter', ',');line = line{1};
        houses{a} = char(line{1});
        distance(a) = str2num(cell2mat(line(2)));
        timestamps(a) = str2num(cell2mat(line(3)));
    end
    avgdist{Number} = mean(distance);
    clear data;
    %calculate how often one house was looked at:
    [uniqueX, ~, J]=unique(cellstr(houses)); %uniqueX = which elements exist in houses, J = to which of the elements in uniqueX does the element in houses correspond 
    occ = histc(J, 1:numel(uniqueX))/30; %histogram bincounts to count the # of occurances of each house
    NumViews = table(uniqueX',occ);
    
    
    
%     %% ------------------------Make Timeline---------------------------------
%     housenumbers = cell(1,len);
%     for a = 1:len-1%convert data
%        if houses{a}(1:2)=='NH'
%            housenumbers{a} = NaN; 
%        else
%            housenumbers{a} = str2num(houses{a}(1:3));
%        end
%     end
%     h = cell2mat(housenumbers);
%     clear a;clear h;
% 
%     list = {,};
%     first = 1;
%     last = 1;
%     label = 'NaN';
%     for e = 2:length(housenumbers)-1
%         if not((isnan(housenumbers{e})&&isnan(housenumbers{e-1}))||housenumbers{e}== housenumbers{e-1})
%             last = e-1;
%             list{end+1,1}=houses{e-1};
%             list{end,2}=[first,last];
%             first = e;
%         end
%     end
%     list{end+1,1}=houses{e-1};
%     list{end,2}=[first,length(housenumbers)-1];
%     clear e; clear first; clear last;
%     %make sure the list is ordered in reverse chronological order (to help identify the last label of a row)
%     [~, order] = sortrows(vertcat(list{:, 2}), [-1 2]);
%     list = list(order, :);
%     %identify unique label and generate vertical positioning
%     [labels, idx, ylabels] = unique(list(:, 1), 'stable');
%     ypatches = max(ylabels) + 1 - [ylabels, ylabels, ylabels-1, ylabels-1]'; 
%     ylabels = max(ylabels) + 1 - ylabels(idx);
%     %generate horizonal positioning
%     xpatches = [vertcat(list{:, 2}), fliplr(vertcat(list{:, 2}))]';
%     xlabels = xpatches(2, idx);
%     %plot
%     figure;
%     color = parula(size(list, 1)); %color distribution
%     patch(xpatches, ypatches, reshape(color, 1, [], 3)); 
%     %text(xlabels-(xlabels-1200), ylabels+0.5, labels, 'fontsize', 10); 
%     text(xlabels+5, ylabels+0.5, labels, 'fontsize', 10);
%     xlabel('Time (Minutes)');
%     set(gca,'XTickLabel',{1 double(len/6/1980) len/6*2/1980 len/6*3/1980 len/6*4/1980 len/6*5/1980 len/1980},'XTick',[1 len/6 len/6*2 len/6*3 len/6*4 len/6*5 len]);%!!only for 30 Mintes Sessions!!
%     grid on
%     saveas(gcf,fullfile([savepath 'Results/'],['Timeline_VP' suj_num '.jpeg']));
%     clear color; clear idx; clear labels; clear order; clear xlabels; clear ylabels; clear xpatches; clear ypatches;




    %% NumView creation
    %Calculate average distance from which each house was looked at and
    %variance in the distance:
    lenHouses = length(uniqueX);
    distances = cell(lenHouses,3);
    distances(1:lenHouses) = uniqueX;
    ix=cellfun('isempty',distances);
    distances(ix)={[0]};
    for a = 1:len
        houseN = J(a);
        distances{houseN,2}=[distances{houseN,2},distance(a)];
    end
    for a = 1:lenHouses
        distances{a,2}=distances{a,2}(distances{a,2}~=0);
        distances{a,3}=var(distances{a,2});
        distances{a,2}=mean(distances{a,2});
    end
    %Save everything in NumViews:------------------------------------------
    NumViews = [NumViews distances];
    NumViews(:,[3])=[];
    remove = isnan(NumViews.Var4);%remove houses that we're 'seen' from further away than the far clip plane
    NumViews(remove,:)=[];
    NumViews.Properties.VariableNames{'Var1'}='House';NumViews.Properties.VariableNames{'Var4'}='DistanceMean';NumViews.Properties.VariableNames{'Var5'}='DistanceVariance';
    %clear uniqueX;clear occ;clear J;
    %Save NumViews as a matlab table:
    
    if WorkPlace==2
    %PC
    current_name = strcat(savepath,'/','NumViewsD_','VP_',num2str(suj_num),'.mat');
    else
    current_name = strcat(savepath,'\','NumViewsD_','VP_',num2str(suj_num),'.mat');
    end

   save(current_name,'NumViews')
end
%clear all;

% %shows how often each number appeared
%    A = occ;
%    x = unique(A);
%    N = numel(x);
% 
%    count = zeros(N,1);
% 
%    for k = 1:N
%       count(k) = sum(A==x(k));
%    end
% 
%    disp([ x(:) count ]);