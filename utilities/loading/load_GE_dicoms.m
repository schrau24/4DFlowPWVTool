function [directory, nframes, res, fov, pixdim, timeres, v, MAG, magWeightVel, angio, vMean, VENC] = ...
    load_GE_dicoms()
% Get and load input directory

[directory] = uigetdir('Select Reconstructed Data');

warning('off','all');
%% grab each parrec and save corresponding data
disp('Loading data')

% the anatomy directory
tmpMAG = dir([directory '\ANATOMY\*']); tmpMAG(1:2) = [];
tmpV1 = dir([directory '\AP FLOW\*']); tmpV1(1:2) = [];
tmpV2 = dir([directory '\LR FLOW\*']); tmpV2(1:2) = [];
tmpV3 = dir([directory '\SI FLOW\*']); tmpV3(1:2) = [];
% read in first header to get some image info
info = dicominfo(fullfile(tmpMAG(end).folder,tmpMAG(end).name));
nframes = info.CardiacNumberOfImages;                                       % number of reconstructed frames
timeres = (60/info.HeartRate*1000-info.TriggerWindow)/nframes;              % temporal resolution, in ms
nslices = length(tmpMAG)/nframes;
fov = double(info.AcquisitionMatrix);                                       % Field of view in cm
fov = squeeze(fov(fov>0)); fov(3) = nslices*info.SpacingBetweenSlices;
fov = fov/10;
VENC = double(info.Private_0019_10cc);                                      % venc, in mm/s

res = double([info.Width info.Height nslices]);                             % number of pixels in row,col,slices
pixdim = [10*fov(1)/res(1) 10*fov(2)/res(2) info.SpacingBetweenSlices];     % the reconstructed resolution

% now loop over all magnitude and velocities and load in dicoms, takes
% awhile
% outer loop is slices, inner is cardiac frames
mag = zeros([res nframes],'single');
vx = zeros([res nframes],'single');
vy = zeros([res nframes],'single');
vz = zeros([res nframes],'single');
count = 0;
h = waitbar(0,'Reading in all dicoms...');
for slc = 1:nslices
    for cf = 1:nframes
        count = count+1;
        mag(:,:,slc,cf) = single(dicomread(fullfile(tmpMAG(count).folder, tmpMAG(count).name)));
        vx(:,:,slc,cf) = single(dicomread(fullfile(tmpV1(count).folder, tmpV1(count).name)));
        vy(:,:,slc,cf) = single(dicomread(fullfile(tmpV2(count).folder, tmpV2(count).name)));
        vz(:,:,slc,cf) = single(dicomread(fullfile(tmpV3(count).folder, tmpV3(count).name)));
        
        if mod(count,100) == 0
            waitbar (count/(nslices*nframes), h);
        end
    end
end
close(h)

vy = -vy;       % flip vy!
v = cat(5,vx,vy,vz); v = permute(v, [1 2 3 5 4]);
clear  vx vy vz

% take the means
vMean = mean(v,5);
MAG = mag./max(mag(:)); clear mag;

% Calculate a Polynomial
[poly_fitx,poly_fity, poly_fitz] = background_phase_correction(mean(MAG,4),vMean(:,:,:,1),vMean(:,:,:,2),vMean(:,:,:,3),VENC);
%%
disp('Correcting Data with Polynomial');
xrange = single( linspace(-1,1,size(MAG,1)));
yrange = single( linspace(-1,1,size(MAG,2)));
zrange = single( linspace(-1,1,size(MAG,3)));
[y,x,z] = meshgrid( yrange,xrange,zrange);

disp('   Vx');
back = evaluate_poly(x,y,z,poly_fitx);
back = single(back);
vMean(:,:,:,1) = vMean(:,:,:,1) - back;
for m = 0 : nframes - 1
    v(:,:,:,1,m+1) = v(:,:,:,1,m+1) - back;
end

disp('   Vy');
back = evaluate_poly(x,y,z,poly_fity);
back = single(back);
vMean(:,:,:,2) = vMean(:,:,:,2) - back;
for m = 0 : nframes - 1
    v(:,:,:,2,m+1) = v(:,:,:,2,m+1) - back;
end

disp('   Vz');
back = evaluate_poly(x,y,z,poly_fitz);
back = single(back);
vMean(:,:,:,3) = vMean(:,:,:,3) - back;
for m = 0 : nframes - 1
    v(:,:,:,3,m+1) = v(:,:,:,3,m+1) - back;
end

[magWeightVel, angio] = calc_angio(MAG, v, VENC);

disp('Load Data finished');
return