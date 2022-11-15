
function [magWeightVel, angio] = calc_angio(MAG, v, Venc)

% time resolved (for contouring later)
Vmag = squeeze(sqrt( sum( v.^2,4)));
Vmag(Vmag > Venc) = Venc;
% 
magWeightVel = 32000*MAG.*sin( pi/2*Vmag / Venc);
mm = mean(MAG,4); 


% time averaged (looks good)
vMean = mean(v,5);
Vmag = squeeze(sqrt( sum( vMean.^2,4)));
Vmag(Vmag > Venc) = Venc;

angio = 32000*mm.*sin( pi/2*Vmag / Venc);


% % calculate softmax angio image
% scaling = 0.05*pi;
% v2 = v/Venc*pi;
% 
% vx = squeeze(v2(:,:,:,1,:));
% vx_sm = sum(exp(vx/scaling).*vx,4)./(sum(exp(vx/scaling),4));
% 
% vy = squeeze(v2(:,:,:,2,:));
% vy_sm = sum(exp(vy/scaling).*vy,4)./(sum(exp(vy/scaling),4));
% 
% vz = squeeze(v2(:,:,:,3,:));
% vz_sm = sum(exp(vz/scaling).*vz,4)./(sum(exp(vz/scaling),4));
% 
% Vmag = squeeze(sqrt(vx_sm.^2 + vy_sm.^2 + vz_sm.^2));
% Vmag(Vmag > pi) = pi;
% angio = 32000*mm.*sin( pi/2*Vmag / pi);
return