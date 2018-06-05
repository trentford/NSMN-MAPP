clear
clc

%CPC spatial information
% R = georasterref;
% R.RasterSize = [32,64];
% R.Latlim = [20.5, 52.5];
% R.Lonlim = [-130.5, -66.5];
% R.ColumnsStartFrom = 'north';
% R.RowsStartFrom = 'west';
% -----------------------------
%NLDAS-2 spatial information
% R = georasterref;
% R.RasterSize = [224,464];
% R.Latlim = [25.0625, 52.9375];
% R.Lonlim = [-124.9375, -67.0625];
% R.ColumnsStartFrom = 'south';
% R.RowsStartFrom = 'west';
%---------------------------------
%SMAP L4 spatial information
% R = georasterref;
% R.RasterSize = [302,643];
% R.Latlim = [24, 51];
% R.Lonlim = [-126, -66];
% R.ColumnsStartFrom = 'north';
% R.RowsStartFrom = 'west';
%-------------------------------
%ECV spatial information
R = georasterref;
R.RasterSize = [121,281];
R.Latlim = [20, 50];
R.Lonlim = [-130, -60];
R.ColumnsStartFrom = 'north';
R.RowsStartFrom = 'west';

%Go Salukis!
yell = 'Go Salukis'

soilm = load('G:\Research-Backup\MAPP Drought\Remote Sensing\ECV\ECV_VWC_1998_2016.mat');
soilm = soilm.soilm;

impacts = importdata('All_Drought_Impacts.xlsx');
impacts = impacts.data;

time = importdata('G:\Research-Backup\MAPP Drought\Remote Sensing\ECV\ECV_Dates.csv');
time = time.data;
clear r1 r2

data = importdata('G:\Research-Backup\MAPP Drought\In Situ\SCAN\SCAN_MS_VWC.csv');
data = data.data;
stations = unique(data(:,1));
%---------------------------------------------------------------
% Vertically interpolate observation data
% data(:,6) = (data(:,7).*(1/6))+(data(:,8).*(5/6));
% data(:,6) = data(:,7);
%-----------------------------------------------------------------
depth = 10;
impacts = impacts(impacts(:,1) == 3, :);

% lat = load('G:\Research-Backup\MAPP Drought\Remote Sensing\SMOS\SMOS_lat.mat');
% lat = lat.lat;
% 
% lon = load('G:\Research-Backup\MAPP Drought\Remote Sensing\SMOS\SMOS_lon.mat');
% lon = lon.lon;
%Average in situ data over the study reigon
avgData = [];
years = unique(data(:,2));
for i = 1:length(years)
    subset = data(data(:,2) == years(i,1), :);
    doys = unique(subset(:,5));
    for ii = 1:length(doys)
        sub = subset(subset(:,5) == ii, :);
        avgData = [avgData; sub(1,2:5) nanmean(sub(:,depth),1)];
        clear sub
    end
    clear subset doys
end

%Compute percentiles from the daily, averaged VWC
[obsPerc] = InSituPercentiles(avgData);

%Match stations with grid cells - Average model data over all grid cells in
%the study region
for i = 1:length(stations)
    subset = data(data(:,1) == stations(i), :);
    [row,col] = latlon2pix(R,subset(1,end-1),subset(1,end));
    gcells(i,1) = stations(i);
    gcells(i,2) = round(row);
    gcells(i,3) = round(col);
    
%     [c,ind] = nanmin(abs(subset(1,end-1)-lat));
%     gcells(i,1) = stations(i);
%     gcells(i,2) = ind(1); clear c ind
%     [c,ind] = (nanmin(abs(subset(1,end)-lon)));
%     gcells(i,3) = ind(1); clear c ind
    
    
    clear row col subset c ind
end

[C,ix,ic] = unique(gcells(:,2:3),'rows','stable');
for i = 1:length(C(:,1))
    avg(:,i) = soilm(C(i,1),C(i,2),:);
end
clear C ix ic
  
model = [time nanmedian(avg,2)];

%Now compute percentiles from daily, averaged model VWC
[modelPerc] = ModelPercentiles(model); 

clear avg gcells model rows stations
%----------------------------------------------------------
%Find the date of each impact and composite 5-day and 7-day VWC percentiles
%from model and in situ lead up to impact day
leads(:,1) = [55; 49; 42; 35; 28; 21; 14; 7];
leads(:,2) = [49; 42; 35; 28; 21; 14; 7; 1];

count = 1;
for i = 1:length(impacts)
    [C,ia,ib] = intersect(obsPerc(:,1:3),impacts(i,2:4),'rows');
    if isempty(ia) == 0
        for ii = 1:length(leads)
            impComp(count,ii) = nanmean(obsPerc(ia-leads(ii,1):ia-leads(ii,2),5));            
        end
        count = count + 1;
    end
    clear C ia ib
end

for i = 1:length(impComp(1,:))
    tpRate(i,1) = (length(impComp(impComp(:,i) <= 0.3, i)))/(length(impComp(:,i)));
    avg(i,1) = nanmean(impComp(impComp(:,i) <= 0.3, i));
end

tpRate'
avg'
trapz(tpRate)

