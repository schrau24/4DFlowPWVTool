function [D, fitObject, dist_total] = calc_pwv(waveforms,dist_total, timeres, PWVcalctype, scale)

% PWV calc type: 1 is cross correlation, 2 is TTF, 3 is Wavelet

% to interp to 1/scale ms, we need how many frames?
nFrames = size(waveforms,2);
nFrames_interp = floor(timeres/scale*nFrames);

% interp params
x = 1:nFrames;
xq = linspace(1,nFrames,nFrames_interp);
METHOD = 'pchip';
plot_steps = false;

% grab first flow waveform
flow_sl1 = smooth(circshift(interp1(x,waveforms(1,:),xq,METHOD),50),round(nFrames_interp/15));
if min(flow_sl1) < 0
    flow_sl1 = flow_sl1 + abs(min(flow_sl1));
end
if min(flow_sl1) > 0
    flow_sl1 = flow_sl1 - (min(flow_sl1));
end

% normalize
flow_sl1 = flow_sl1./max(abs(flow_sl1));

% find the systolic upslope, TTF and wavelet
[maxFl, indMax1] = max(flow_sl1);
pts1 = find(flow_sl1 >= 0.2*maxFl & flow_sl1 <= 0.8*maxFl);
pts1 = pts1(pts1<indMax1);

if PWVcalctype == 2     % TTF
    fitObject = polyfit(pts1,flow_sl1(pts1),1);
    D(1) = -fitObject(2) / fitObject(1);
end

% for wavelet
[y1,PERIOD,SCALE,COI,DJ, PARAMOUT, K] = contwt(flow_sl1,.001,1,.125,[],[],'dog',4);
f1 = 1./PERIOD;
ind = [];
% loop over all flow curves
for slice = 2:size(waveforms,1)

    flow_sl2 = smooth(circshift(interp1(x,waveforms(slice,:),xq,METHOD),50),round(nFrames_interp/15));
    
    if min(flow_sl2) < 0
        flow_sl2 = flow_sl2 + abs(min(flow_sl2));
    end
    if min(flow_sl2) > 0
        flow_sl2 = flow_sl2 - (min(flow_sl2));
    end

%     normalize
    flow_sl2 = flow_sl2./max(abs(flow_sl2));
    
    % TTF and wavelet
    % find the systolic upslope (20% < x < 80%) at or after maxInd of
    % flow_sl1
    [maxFl, indMax] = max(flow_sl2(indMax1:end));
    pts2 = find(flow_sl2 >= 0.2*maxFl & flow_sl2 <= 0.8*maxFl);
    pts2 = pts2(pts2<indMax+indMax1);
    
    switch PWVcalctype
        case 1      % cross-correlation
%             tempDelay = abs(finddelay(flow_sl1,flow_sl2)) * scale;
            % delay calculated on systolic upstroke only as per Saintounge
            % paperclc; clear all; close all;

            sl1 = zeros(size(flow_sl1));sl1(pts1) = flow_sl1(pts1);
            sl2 = zeros(size(flow_sl2));sl2(pts2) = flow_sl2(pts2);
            tempDelay = abs(finddelay(sl1,sl2)) * scale;
            if (dist_total(slice-1) > 20 && tempDelay < 1) || ...
                    (dist_total(slice-1) > 40 && tempDelay < 2) || ...
                    (dist_total(slice-1) > 60 && tempDelay < 3) || ...
                    (dist_total(slice-1) > 80 && tempDelay < 4) || ...
                    (dist_total(slice-1) > 100 && tempDelay < 6)
                % force the delays to make sense along vessel
                ind = cat(1,ind,slice);     % to remove at the end for fitting
            end
            D(slice-1) = tempDelay * scale;
            
        case 2      % TTF
            fitObject = polyfit(pts2,flow_sl2(pts2),1);
            D(slice) = -fitObject(2) / fitObject(1);
            
        case 3      % Wavelet time delay estimation
            [y2,~,~,~,~,~,~] = contwt(flow_sl2,.001,1,.125,[],[],'dog',4);
            % the complex cross-spectrum
            y = y1.*conj(y2);
            
            % find frequencies 1:10 Hz and keep only those
            ind_f = find(f1 >= 1 & f1 <= 10);
            
            % the points between foot of flow_sl1 and peak of flow_sl2
            pts = unique([pts1; pts2]);
            
            % calculate y_hat
            y_hat = y/(sum(sum(abs(y))));
            
            psi = atan(imag(y))./repmat(f1'*2*pi,[1 size(y,2)]);
%             aa = (y_hat(ind_f,pts)).*psi(ind_f,pts);        % eqn 7 in Bargiotas paper
            aa =    dot(y_hat(ind_f,pts),psi(ind_f,pts));  % eqn 7 in Bargiotas paper
            
             tempDelay = abs(sum(aa(:)))*1000*4 * scale;
            if (dist_total(slice-1) > 20 && tempDelay < 1) || ...
                    (dist_total(slice-1) > 40 && tempDelay < 2) || ...
                    (dist_total(slice-1) > 60 && tempDelay < 3) || ...
                    (dist_total(slice-1) > 80 && tempDelay < 4) || ...
                    (dist_total(slice-1) > 100 && tempDelay < 6)
                % force the delays to make sense along vessel
                ind = cat(1,ind,slice);     % to remove at the end for fitting
            end
            D(slice - 1) = tempDelay;
            
    end

                
    if plot_steps
        figure(400); clf;
        plot(flow_sl1,'r'); hold on;
        plot(pts1,flow_sl1(pts1),'r*')
        plot(flow_sl2,'b')
        plot(pts2,flow_sl2(pts2),'b*')
        legend('flow 1','flow 1 upslope', 'flow 2','flow 2 upslope')
        xlabel('time (ms)')
        title(['Point ' num2str(slice) ', delay = ' num2str(D(slice-1)) ' ms'])
        pause(0.1)
    end
    
end
% remove unused slices
D(ind-1) = [];
dist_total(ind-1) = [];

% now do a linear fit to calculate slope, save also rmse
[fitObject, ~] = polyfit(dist_total,D,1);
