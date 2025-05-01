function Eventlist = MPEE_V1(PMinput, varargin)
% EXTRACTEXTREMEPRECIPITATION Detects 3‑D extreme precipitation events in a
% precipitation data cube (M × N × T).
%
%   EVENTS = EXTRACTEXTREMEPRECIPITATION(PMALL) extracts extreme
%   precipitation events using default parameters.
%
%   EVENTS = EXTRACTEXTREMEPRECIPITATION(PMALL,'PARAM',VALUE,...) overrides
%   the defaults.  Accepted parameters are:
%       'SeedThreshold'        – High threshold (mm h‑1) for precipitation
%                                cores (default 1.5).
%       'PeripheralThreshold'  – Lower threshold (mm h‑1) for peripheral
%                                growth (default 0.1).
%       'FilterKernel'         – 3‑element vector with the size of the 3‑D
%                                mean filter (default [5 5 5]).
%       'DilationIterations'   – Number of 3‑D dilations used to grow
%                                events (default 50).
%       'DilationRadius'       – Radius (pixels) of the spherical structuring
%                                element for dilation (default 1).
%       'MinEventSize'         – Minimum number of voxels required to keep
%                                an event (default 200).
%
%   Output
%       EVENTS – array of structures with fields:
%           id  – sequential identifier
%           MS  – row indices of voxels (1‑based)
%           NS  – column indices of voxels (1‑based)
%           TS  – time indices of voxels (1‑based)
%           PS  – precipitation values (mm h‑1) for the voxels
%
%   Example
%       load('CONUS404_1980_subset.mat','PM');
%       ev = extractExtremePrecipitation(PM,'SeedThreshold',2);
%
%   Siyu Zhu, 01‑May‑2025
%   University of Oklahoma
%  
%   Please cite:
%   XXX
% ---------------------------------------------------------------------

% --------------------------- Input parsing ---------------------------
assert(ndims(PMinput)==3,'PMall must be a 3‑D array.')
p = inputParser;
addRequired(p,'PMall');
addParameter(p,'SeedThreshold',1.5,@isnumeric);
addParameter(p,'PeripheralThreshold',0.1,@isnumeric);
addParameter(p,'FilterKernel',[5 5 5],@(x) isnumeric(x)&&numel(x)==3);
addParameter(p,'DilationIterations',50,@isscalar);
addParameter(p,'DilationRadius',1,@isscalar);
addParameter(p,'MinEventSize',200,@isscalar);
parse(p,PMinput,varargin{:});
S = p.Results;

% -------------------------- Pre‑processing ---------------------------
% Smooth field with a 3‑D mean filter to suppress single‑pixel noise
kernel = ones(S.FilterKernel)/prod(S.FilterKernel);
PMs = imfilter(double(PMinput),kernel,'symmetric');

% ------------------------ Threshold masking -------------------------
bwPeripheral = PMs >= S.PeripheralThreshold;  % broad candidate mask
bwCore       = PMs >= S.SeedThreshold;        % seed (core) mask

% ----------------------- Seed‑based grouping ------------------------
labelM = bwlabeln(bwCore,6);                  % 6‑connectivity in 3‑D

% ------------------------- Region growing ---------------------------
SE3d = strel('sphere',S.DilationRadius);
for iIter = 1:S.DilationIterations
    currMask  = labelM > 0;
    dilated   = imdilate(labelM,SE3d);
    growMask  = dilated > 0 & ~currMask & bwPeripheral;
    labelM(growMask) = dilated(growMask);     % inherit labels
end

% --------------------- Extract event properties ---------------------
props  = regionprops(labelM,'PixelList','Area');
keepID = find([props.Area] >= S.MinEventSize);

nKeep  = numel(keepID);
Eventlist = struct('id',cell(1,nKeep),...
                'MS',cell(1,nKeep),...
                'NS',cell(1,nKeep),...
                'TS',cell(1,nKeep),...
                'PS',cell(1,nKeep));

for k = 1:nKeep
    idx = keepID(k);
    vox = props(idx).PixelList;
    ms  = vox(:,2);
    ns  = vox(:,1);
    ts  = vox(:,3);
    lin = sub2ind(size(PMinput),ms,ns,ts);

    Eventlist(k).id = k;
    Eventlist(k).MS = ms;
    Eventlist(k).NS = ns;
    Eventlist(k).TS = ts;
    Eventlist(k).PS = PMinput(lin);
end

end