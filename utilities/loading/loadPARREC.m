function [directory, nframes, res, fov, pixdim, timeres, v, MAG, magWeightVel, angio, vMean, VENC, ori] = ...
    loadPARREC()
% Get and load input directory

[filename,directory] = uigetfile('*.rec','Select Reconstructed Data');
fBase = filename(1:end-5);

warning('off','all');
%% grab each parrec and save corresponding data
disp('Loading data')
PARRECFILE = fullfile(directory,[fBase, '1.rec']);
[IMG1,~] = readrec_V4_2(PARRECFILE, 'noscale');
IMG1 = double(IMG1);
vx = squeeze(IMG1(:,:,:,:,:,2,:))-2048;
mag1 = squeeze(IMG1(:,:,:,:,:,1,:));

PARRECFILE = fullfile(directory,[fBase, '2.rec']);
[IMG2,~] = readrec_V4_2(PARRECFILE,'noscale');
IMG2 = double(IMG2);
vy = squeeze(IMG2(:,:,:,:,:,2,:))-2048;
mag2 = squeeze(IMG2(:,:,:,:,:,1,:));

PARRECFILE = fullfile(directory,[fBase, '3.rec']);
[IMG3,header] = readrec_V4_2(PARRECFILE, 'noscale');
IMG3 = double(IMG3);
vz = squeeze(IMG3(:,:,:,:,:,2,:))-2048;
mag3 = squeeze(IMG3(:,:,:,:,:,1,:));
warning('on','all');

MAG = mean(cat(5,mag1,mag2,mag3),5);

vy = -vy;       % flip vy!
v = cat(5,vx,vy,vz); v = permute(v, [1 2 3 5 4]);
clear mag1 mag2 mag3 IMG1 IMG2 IMG3 vx vy vz

nframes = header.nphases;                                       % number of reconstructed frames
timeres = max(header.tbl(:,header.tblcols.ttime))/nframes;      % temporal resolution, in ms
fov = header.fov;                                               % Field of view in cm
res = round([header.nrows header.ncols header.nslices]);        % number of pixels in row,col,slices
VENC = max(header.pevelocity)*10;                               % venc, in mm/s
pixdim = header.pixdim;                                         % the reconstructed resolution
ori = header.tbl(1,26);                                         % orientation number (1 - axial, 2 - sagittal, 3 - coronal)

% % load the mask .mat file, if available. then check to see if some of the
% % trailing frames have high acceleration factors (R>30). Remove these
% tmp = dir([directory '/*mask.mat']);
% if ~isempty(tmp)
%     disp('Checking for cardiac frames with high acceleration (R>30)')
%     load(fullfile(directory,tmp(1).name));
%     mask = squeeze(mask(:,:,:,1,:,:,:,:,:,1));  % cardiac phase is 4th dim
%     
%     % loop through cardiac phases and calculate the R for each
%     for cp = 1:nframes
%         m = mask(:,:,:,cp);
%         R(cp) = numel(m)/sum(m(:)) /4*pi;
%     end
%     
%     idxToRemove = find(R>30);
%     
%     % remove these
%     MAG(:,:,:,idxToRemove) = [];
%     v(:,:,:,:,idxToRemove) = [];
%     
%     % update other parameters
%     nframes = nframes - numel(idxToRemove);
%     disp(['Removing the last ' num2str(numel(idxToRemove)) ' cardiac frames'])
% end
    
% scale velocity
v = (v./2048)*VENC;

% take the means
vMean = mean(v,5);
MAG = MAG./max(MAG(:));

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