%PartList = {1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,18,19,20,21,22,24,25,26,28,29};
PartList = {34};
Number = length(PartList);
avgdist = cell(1,Number);
for ii = 1:Number
    suj_num = cell2mat(PartList(ii));
    if suj_num < 10
        file = strcat('ViewedHouses_VP',num2str(0),num2str(suj_num),'.txt');
    else
        file = strcat('ViewedHouses_VP',num2str(suj_num),'.txt');
    end 
    disp(file);
    data = fopen(file);
    data = textscan(data,'%s','delimiter', '\n');
    data = data{1};
    data = table2array(cell2table(data));
    %initialize fields
    len = int16(length(data));
    houses = cell(1,len);
    distance = zeros(1,len);
    timestamps = zeros(1,len);
    for a = 1:double(len)
        line = textscan(data{a},'%s','delimiter', ',');line = line{1};
        houses{a} = char(line{1});
        distance(a) = str2num(cell2mat(line(2)));
        timestamps(a) = str2num(cell2mat(line(3)));
    end
    avgdist = mean(distance);
    clear data;
    %calculate how often one house was looked at:
    [uniqueX, ~, J]=unique(cellstr(houses)); %uniqueX = which elements exist in houses, J = to which of the elements in uniqueX does the element in houses correspond 
    occ = histc(J, 1:numel(uniqueX)); %histogram bincounts to count the # of occurances of each house
    NumViews = table(uniqueX',occ);

    %make timeline---------------------------------------------------------
    housenumbers = cell(1,len);
    for a = 1:len-1%convert data
       if houses{a}(1:2)=='NH'
           housenumbers{a} = NaN; 
       else
           housenumbers{a} = str2num(houses{a}(1:3));
       end
    end
    h = cell2mat(housenumbers);
    clear a;clear h;

    list = {,};
    first = 1;
    last = 1;
    label = 'NaN';
    for e = 2:length(housenumbers)-1
        if not((isnan(housenumbers{e})&&isnan(housenumbers{e-1}))||housenumbers{e}== housenumbers{e-1})
            last = e-1;
            list{end+1,1}=houses{e-1};
            list{end,2}=[first,last];
            first = e;
        end
    end
    list{end+1,1}=houses{e-1};
    list{end,2}=[first,length(housenumbers)-1];
    clear e; clear first; clear last;
    %%Create timeline
    %make sure the list is ordered in reverse chronological order (to help identify the last label of a row)
    [~, order] = sortrows(vertcat(list{:, 2}), [-1 2]);
    list = list(order, :);
    %identify unique label and generate vertical positioning
    [labels, idx, ylabels] = unique(list(:, 1), 'stable');
    ypatches = max(ylabels) + 1 - [ylabels, ylabels, ylabels-1, ylabels-1]'; 
    ylabels = max(ylabels) + 1 - ylabels(idx);
    %generate horizonal positioning
    xpatches = [vertcat(list{:, 2}), fliplr(vertcat(list{:, 2}))]';
    xlabels = xpatches(2, idx);
    %plot
    figure;
    color = parula(size(list, 1)); %color distribution
    patch(xpatches, ypatches, reshape(color, 1, [], 3)); 
    %text(xlabels-(xlabels-1200), ylabels+0.5, labels, 'fontsize', 10); 
    text(xlabels+5, ylabels+0.5, labels, 'fontsize', 10);
    xlabel('Time (seconds/30)');
    grid on
    clear color; clear idx; clear labels; clear order; clear xlabels; clear ylabels; clear xpatches; clear ypatches;
    %----------------------------------------------------------------------
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
    %save everything in NumViews:------------------------------------------
    NumViews = [NumViews distances];
    NumViews(:,[3])=[];
    remove = isnan(NumViews.Var4);%remove houses that we're 'seen' from further away than the far clip plane
    NumViews(remove,:)=[];
    NumViews.Properties.VariableNames{'Var1'}='House';NumViews.Properties.VariableNames{'Var4'}='DistanceMean';NumViews.Properties.VariableNames{'Var5'}='DistanceVariance';
    clear uniqueX;clear occ;clear J;
    %Save NumViews as a matlab table:
    if suj_num < 10
        current_name = strcat('D:/v.kakerbeck/Tracking/ViewedHouses/','NumViewsD_','VP_',num2str(0),num2str(suj_num),'.mat');
    else
        current_name = strcat('D:/v.kakerbeck/Tracking/ViewedHouses/','NumViewsD_','VP_',num2str(suj_num),'.mat');
    end 
    save(current_name,'NumViews')
end