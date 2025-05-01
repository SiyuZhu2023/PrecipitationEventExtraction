clear all
clc

%% load the data
MEPPpath = '.\MPEE_V1';
addpath(MEPPpath);


DataPath = '.\TestData';
DataName1 = 'PMall_Vtest.mat';
DataName2 = fullfile(DataPath,DataName1);
load(DataName2,'PMall');

%% core function
PMinput = PMall;
Eventlist = MPEE_V1(PMinput);

% For more details please check Function file
%   Eventlist = MPEE_V1(PMinput, ...
%       'SeedThreshold',      1.5,    ... 
%       'PeripheralThreshold',0.1,    ... 
%       'FilterKernel',       [5 5 5],... 
%       'DilationIterations', 50,     ... 
%       'DilationRadius',     1,      ... 
%       'MinEventSize',       200);       


%% --- Print field explanations ------------------------------------------
fprintf('\nEventlist contains the following fields:\n');
fprintf('  id : sequential identifier for the event\n');
fprintf('  MS : row indices (latitude grid) of pixels in the event\n');
fprintf('  NS : column indices (longitude grid) of pixels in the event\n');
fprintf('  TS : time indices (e.g., hourly) of voxels in the event\n');
fprintf('  PS : precipitation intensity (mm h^-1) for each voxel\n');

%% --- Display the first five events as examples --------------------------
nShow = min(5, numel(Eventlist));
fprintf('\n----- Example events (first %d) -----\n', nShow);

for k = 1:nShow
    ev = Eventlist(k);

    % 统计信息
    nPix   = numel(ev.PS);
    maxP   = max(ev.PS);
    meanP  = mean(ev.PS);
    durT   = numel(unique(ev.TS));   % 时长（时间层数）
    minLat = min(ev.MS);
    maxLat = max(ev.MS);
    minLon = min(ev.NS);
    maxLon = max(ev.NS);

    fprintf(['Event %3d | pixels=%5d | maxP=%6.2f mm/h | meanP=%5.2f mm/h | ', ...
             'duration=%3d steps | MS=[%d–%d] | NS=[%d–%d]\n'], ...
             ev.id, nPix, maxP, meanP, durT, minLat, maxLat, minLon, maxLon);
end




