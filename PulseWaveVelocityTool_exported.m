classdef PulseWaveVelocityTool_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PulseWaveVelocityToolUIFigure   matlab.ui.Figure
        TabGroup                        matlab.ui.container.TabGroup
        DataLoadingandPreprocessingTab  matlab.ui.container.Tab
        LoadDataPanel                   matlab.ui.container.Panel
        LoadReconstructedDataButton     matlab.ui.control.Button
        DataDirectoryEditFieldLabel     matlab.ui.control.Label
        DataDirectoryEditField          matlab.ui.control.EditField
        LoadSegmentationDicomsButton    matlab.ui.control.Button
        SegmentationDirectoryEditFieldLabel  matlab.ui.control.Label
        SegmentationDirectoryEditField  matlab.ui.control.EditField
        ScanInfoTable                   matlab.ui.control.Table
        CropPanel                       matlab.ui.container.Panel
        XrangeEditFieldLabel            matlab.ui.control.Label
        XrangeEditField                 matlab.ui.control.EditField
        toXEditFieldLabel               matlab.ui.control.Label
        toXEditField                    matlab.ui.control.EditField
        toYEditFieldLabel               matlab.ui.control.Label
        toYEditField                    matlab.ui.control.EditField
        YrangeEditFieldLabel            matlab.ui.control.Label
        YrangeEditField                 matlab.ui.control.EditField
        toZEditFieldLabel               matlab.ui.control.Label
        toZEditField                    matlab.ui.control.EditField
        ZrangeEditFieldLabel            matlab.ui.control.Label
        ZrangeEditField                 matlab.ui.control.EditField
        AxesX                           matlab.ui.control.UIAxes
        AxesY                           matlab.ui.control.UIAxes
        AxesZ                           matlab.ui.control.UIAxes
        CropButton                      matlab.ui.control.Button
        CropButton_2                    matlab.ui.control.Button
        CropButton_3                    matlab.ui.control.Button
        DVisualizationPanel             matlab.ui.container.Panel
        View3D                          matlab.ui.control.UIAxes
        RotateLeft                      matlab.ui.control.Button
        RotateRight                     matlab.ui.control.Button
        Rotate                          matlab.ui.control.Label
        RotateDown                      matlab.ui.control.Button
        RotateUp                        matlab.ui.control.Button
        ProcessingPanel                 matlab.ui.container.Panel
        UnwrapVelocity                  matlab.ui.control.Button
        CenterlineExtraction            matlab.ui.control.Button
        FlowandPulseWaveVelocityTab     matlab.ui.container.Tab
        SegmentationAndCenterline       matlab.ui.container.Panel
        View3D_2                        matlab.ui.control.UIAxes
        Reset3DviewButton               matlab.ui.control.Button
        BranchNumberTitle               matlab.ui.control.Label
        BranchNumberLabel               matlab.ui.control.Label
        BranchNumbers                   matlab.ui.control.EditField
        CheckcenterlinecalculateflowButton  matlab.ui.control.Button
        FlipcenterlinepointsButton      matlab.ui.control.Button
        PWVPointsTitle                  matlab.ui.control.Label
        PWVPointsLabel                  matlab.ui.control.Label
        PWVPoints                       matlab.ui.control.EditField
        PlotWaveformsButton             matlab.ui.control.Button
        WaveformsDisplay                matlab.ui.control.UIAxes
        CalculatePWV                    matlab.ui.control.Button
        PWVType                         matlab.ui.control.DropDown
        PWVDisplayTitle                 matlab.ui.control.Label
        PWVDisplay                      matlab.ui.control.EditField
        PWVCalcDisplay                  matlab.ui.control.UIAxes
        SavingTitle                     matlab.ui.control.Label
        SaveResultsCallback             matlab.ui.control.Button
        SaveName                        matlab.ui.control.DropDown
        timeResolvedSegCheckbox         matlab.ui.control.CheckBox
        ResetWorkSpace                  matlab.ui.container.Tab
        CleardataandrestartanalysisButton  matlab.ui.control.Button
    end

    
    properties (Access = private)
        directory;                  % the data directory
        segDirectory;               % the directory for dicoms from a pre-defined manual segmentation
        v;                          % the 5D velocity matrix (X x Y x Z x t x v)
        nframes;                    % the reconstrucated cardiac time frames
        res;                        % image dimensions in X, Y, and Z
        fov;                        % the acquired field of view, in cm
        pixdim;                     % resolution in X, Y, and Z, in mm
        ori;                        % orientation (1-axial, 2-sagittal, 3-coronal)
        timeres;                    % temporal resolution (per cardiac frame), in ms
        MAG;                        % the 4D magitude matrix (X x Y x Z x t)
        magWeightVel;               % the calculated magnitude weighted velocity
        angio;                      % a maximum intensity PCMRA
        vMean;                      % the mean velocity over time
        VENC;                       % velocity encoding, in mm/s
        segment;                    % the segmentation, updated throughout
        isSegmentationLoaded = 0;   % is the manual segmentation loaded
        isCropped = 0;              % have we performed any cropping?
        mask;                       % the mask
        isRawDataCropped;           % have we cropped the raw data yet?
        aorta_seg;                  % the specific aorta segmentation
        
        branchList;                 % the list of all unique branches following centerline extraction
        branchActual;               % the chosen branch for PWV measurements
        flowPerHeartCycle_vol;      % the resulting flow over the cardiac cycle in the aorta_seg voluem
        flowPulsatile_vol;          % the pulsatile flow waveforms in the aorta_seg volume 
        hpatch1;                    % the initial 3D patch for 3D vis
        hpatch2;                    % the segmentation 3D patch for 3D vis
        rotAngles;                  % the rotation angles used for viewing, can be changed by viewer
        
    end
    
    methods (Access = private)
        
        function View3DSegmentation(app)
            
            cla(app.View3D);
            ss = smooth3(app.segment);
            
            app.hpatch1 = patch(app.View3D, isosurface(ss,.5),'FaceColor','red','EdgeColor', 'none','FaceAlpha',0.35);
            reducepatch(app.hpatch1 ,0.6);
            if (app.isSegmentationLoaded)
                aa = smooth3(app.aorta_seg);
                hold(app.View3D,'on')
                app.hpatch2 = patch(app.View3D, isosurface(aa,.5),'FaceColor',[0.9 0.9 0.9],'EdgeColor', 'none','FaceAlpha',0.75);
                reducepatch(app.hpatch2,0.6);
            end
            
            % Make it all look good
            camlight(app.View3D);
            lighting(app.View3D,'gouraud');
            view(app.View3D, [0 0 -1]);
            daspect(app.View3D,[1 1 1])
            m_xstart = 1; m_ystart = 1; m_zstart = 1;
            m_xstop = size(app.segment,1); m_ystop = size(app.segment,2); m_zstop = size(app.segment,3);
            xlim(app.View3D,[m_ystart m_ystop]);
            ylim(app.View3D,[m_xstart m_xstop]);
            axis(app.View3D,'off');
            
            % if rotation angles are non-zero, rotate now
            if sum(app.rotAngles) > 0
                rotate(app.hpatch1,[1 0 0], app.rotAngles(1))
                rotate(app.hpatch1,[0 1 0], app.rotAngles(2))
                if (app.isSegmentationLoaded)
                    rotate(app.hpatch2,[1 0 0], app.rotAngles(1))
                    rotate(app.hpatch2,[0 1 0], app.rotAngles(2))
                end
            end
        end
        
        function reset3DSegmentationAndCenterline(app)
            % Initialize figure
            colorbar(app.View3D_2,'off')
            cla(app.View3D_2);
            
            if app.isSegmentationLoaded
                ss = smooth3(app.aorta_seg);
                hpatch = patch(app.View3D_2,isosurface(ss,0.5),'FaceAlpha',0.20);
                reducepatch(hpatch,0.6);
                set(hpatch,'FaceColor',[0.7 0.7 0.7],'EdgeColor', 'none','PickableParts','none');
            else
                ss = smooth3(app.segment);
                hpatch = patch(app.View3D_2,isosurface(ss,0.5),'FaceAlpha',0.20);
                reducepatch(hpatch,0.6);
                set(hpatch,'FaceColor',[0.7 0.7 0.7],'EdgeColor', 'none','PickableParts','none');
            end
            
            unqBranches = unique(app.branchList(:,4));
            c = lines(length(unqBranches));
            for b = 1:length(unqBranches)
                % extract coordinates for branch and plot
                currBranch = find(app.branchList(:,4) == b);
                hline(b) = line(app.View3D_2, ...
                    app.branchList(currBranch,2),app.branchList(currBranch,1),app.branchList(currBranch,3),...
                    'Color',c(b,:),'Marker','.','MarkerSize',12,'LineStyle','none');
            end
            
            % make it look good
            axis(app.View3D_2, 'vis3d')
            axis(app.View3D_2, 'off')
            colormap(app.View3D_2,'colorcube')
            camlight(app.View3D_2);
            lighting(app.View3D_2,'gouraud');
            view(app.View3D_2, [0 0 -1]);
            daspect(app.View3D_2,[1 1 1])
            % to update
            m_xstart = 1; m_ystart = 1; m_zstart = 1;
            m_xstop = app.res(1); m_ystop = app.res(2); m_zstop = app.res(3);
            xlim(app.View3D,[m_ystart m_ystop]);
            ylim(app.View3D,[m_xstart m_xstop]);
            zlim(app.View3D,[m_zstart m_zstop]);
            
            % use rotation angles from first screen to update the view
            % angle
            rotate(hpatch,[1 0 0], app.rotAngles(1))
            rotate(hpatch,[0 1 0], app.rotAngles(2))
            rotate(hline,[1 0 0], app.rotAngles(1))
            rotate(hline,[0 1 0], app.rotAngles(2))
            
            % Put the number labels on the CenterlinePlot
            numString_val = num2str(unqBranches);
            for i = 1:length(unqBranches)
                %find rows where textint = branchList(:,4) and take the mid point
                temp = app.branchList(app.branchList(:,4) == i,1:3);
                textLoc(i,1:3) = temp(round(size(temp,1)/2),:);
            end
            Ntxt = text(app.View3D_2,textLoc(unqBranches,2)+1,textLoc(unqBranches,1)+1,textLoc(unqBranches,3)+1,...
                numString_val,'Color','k','FontSize',14,'FontWeight', 'bold','Margin', 1,...
                'HitTest','off','PickableParts','none');
            
            % update view angle
            rotate(Ntxt,[1 0 0], app.rotAngles(1))
            rotate(Ntxt,[0 1 0], app.rotAngles(2))
            
            if size(app.branchList,4) > 1
                app.BranchNumbers.Value = ['1:' num2str(length(unqBranches))];
            else
                app.BranchNumbers.Value = '1';
            end
        end
        
        function viewAortaSegwFlows(app)
            % indices for flow plotting
            x = app.branchActual(:,1); y = app.branchActual(:,2); z = app.branchActual(:,3);
            index = sub2ind(size(app.aorta_seg),x,y,z);
            
            %reset figure
            cla(app.View3D_2);
            colorbar(app.View3D_2,'off');
            
            hpatch = patch(app.View3D_2,isosurface(smooth3(app.segment),0.5),'FaceAlpha',0.25);
            reducepatch(hpatch,0.6);
            set(hpatch,'FaceColor',[0.7 0.7 0.7],'EdgeColor', 'none','PickableParts','none');
            cdata = app.flowPerHeartCycle_vol(index);
            
            hSurface = surface(app.View3D_2,'XData',[y(:) y(:)],'YData',[x(:) x(:)],'ZData',[z(:) z(:)],...
                'CData',[cdata(:) cdata(:)],'FaceColor','none','EdgeColor','flat',...
                'Marker','.','MarkerSize',12);
            
            caxis(app.View3D_2,[min(cdata) max(cdata)]);
            colormap(app.View3D_2,jet)
            cbar = colorbar(app.View3D_2);
            caxis(app.View3D_2,[0 0.8*max(app.flowPerHeartCycle_vol(:))])
            set(get(cbar,'xlabel'),'string','Flow (mL/cycle)','fontsize',16,'Color','black');
            set(cbar,'FontSize',16,'color','black','Location','west');
            
            % make it look good
            axis(app.View3D_2, 'vis3d')
            axis(app.View3D_2, 'off')
            camlight(app.View3D_2);
            lighting(app.View3D_2,'gouraud');
            view(app.View3D_2, [0 0 -1]);
            daspect(app.View3D_2,[1 1 1])
            %to update
            m_xstart = 1; m_ystart = 1; m_zstart = 1;
            m_xstop = size(app.segment,1); m_ystop = size(app.segment,2); m_zstop = size(app.segment,3);
            xlim(app.View3D,[m_ystart m_ystop]);
            ylim(app.View3D,[m_xstart m_xstop]);
            zlim(app.View3D,[m_zstart m_zstop]);
            
            % update view angles
            rotate(hpatch,[1 0 0], app.rotAngles(1))
            rotate(hpatch,[0 1 0], app.rotAngles(2))
            rotate(hSurface,[1 0 0], app.rotAngles(1))
            rotate(hSurface,[0 1 0], app.rotAngles(2))
            
            % Put the number labels on the CenterlinePlot
            str = app.PWVPoints.Value;
            eval(['ptRange=[' str '];']);
            textint = ptRange(1:5:end);
            numString_val = num2str(textint);
            numString_val = strsplit(numString_val);
            
            c = winter(length(textint));
            for C = 1:length(textint)
                Ntxt(C) = text(app.View3D_2,app.branchActual(textint(C),2)-3,app.branchActual(textint(C),1)-2,app.branchActual(textint(C),3)+1,numString_val{C},...
                    'Color','k','HorizontalAlignment','right',...
                    'FontSize',14,'FontWeight','Bold','HitTest','off','PickableParts','none');
            end
            
            % update view angle
            rotate(Ntxt,[1 0 0], app.rotAngles(1))
            rotate(Ntxt,[0 1 0], app.rotAngles(2))
        end
        
        function plotFlowWaveforms(app)
            
            % grab waveforms
            x = app.branchActual(:,1); y = app.branchActual(:,2); z = app.branchActual(:,3);
            index = sub2ind(size(app.aorta_seg),x,y,z);
            waveforms = app.flowPulsatile_vol(index,:);
            
            viewAortaSegwFlows(app);
            
            % parse points
            str = app.PWVPoints.Value;
            eval(['ptRange=[' str '];']);
            waveforms = waveforms(ptRange,:);
            
            % plot
            card_time = [0:app.nframes-1]*app.timeres;
            c = winter(size(waveforms,1));
            if size(waveforms,1) > 5
                alpha = linspace(0.8,0.3,size(waveforms,1));
            else
                alpha = 0.8*ones(1,size(waveforms,1));
            end
            c = cat(2,c,alpha');
                        
            colorbar(app.WaveformsDisplay,'off')
            cla(app.WaveformsDisplay);
            hold(app.WaveformsDisplay,'on');
            for ii = 1:size(waveforms,1)
                plot(app.WaveformsDisplay,card_time,smooth(waveforms(ii,:)','sgolay'),'Color',c(ii,:),...
                    'LineWidth',2)
            end
            xlim(app.WaveformsDisplay,[0 max(card_time)])
            app.WaveformsDisplay.XLabel.String = 'Cardiac Time (ms)';
            app.WaveformsDisplay.YLabel.String = 'Flow (mL/s)';
            
            if numel(ptRange) > 1
                colormap(app.WaveformsDisplay,winter);
                cbar = colorbar(app.WaveformsDisplay);
                set(get(cbar,'title'),'string','Point number','fontsize',16,'Color','black');
                set(cbar,'FontName','Calibri','FontSize',10,'color','black');
            end
            
            % display a max of 5 points on cbar
            if size(waveforms,1) < 11
                cbar.Ticks = linspace(0, 1, size(waveforms,1)) ;
                cbar.TickLabels = num2cell(ptRange);
            else
                cbar.Ticks = linspace(0, 1, 5) ;
                cbar.TickLabels = num2cell(round(linspace(double(min(ptRange)),double(max(ptRange)),5)));
            end
            hold(app.WaveformsDisplay,'off');
        end
        
        function maskSz = cropImage(app,img)
            
        choice = 0;
        while choice == 0
            cropFig = figure(100);
            set(cropFig,'Units','normalized');
            set(cropFig,'Position',[0.0016 0.0481 0.4969 0.8454])
            set(cropFig,'Name','Draw rectangle to crop image')
            
            % View MIP
            imagesc(img);
            colormap('gray')
            axis equal off
            daspect([1 1 1]);
            
            h = drawrectangle(gca,'DrawingArea',[1 1 size(img,2)-1 size(img,1)-1]);
            maskSz = h.Position;
            maskSz(1:2) = floor(maskSz(1:2));
            maskSz(3:4) = ceil(maskSz(3:4));
            maskSz(maskSz<1) = 1;
            mask = zeros(size(img));
            mask(maskSz(2):(maskSz(2)+maskSz(4)), maskSz(1):(maskSz(1)+maskSz(3))) = 1;
            imgCropped = img.*mask;
            
            clf(cropFig);
            imagesc(imgCropped);
            colormap('gray')
            axis equal off
            daspect([1 1 1]);
            set(cropFig,'Name','Cropped image')
            choice = checkCrop;
            
        end
        % if cancel, reset the mask and img
        if choice == 2
            maskSz = [1 1 size(img,2)-1 size(img,1)-1];
        end
        % update cropped state
        if choice == 1
            app.isCropped = 1;
        end
            
            close(cropFig)
        end
        
        function updateMIPs(app, m_xstart, m_xstop, m_ystart, m_ystop, m_zstart, m_zstop)
            
            % update for ROI drawing
            app.XrangeEditField.Value = num2str(m_xstart);
            app.toXEditField.Value = num2str(m_xstop);
            
            app.YrangeEditField.Value = num2str(m_ystart);
            app.toYEditField.Value = num2str(m_ystop);
            
            app.ZrangeEditField.Value = num2str(m_zstart);
            app.toZEditField.Value = num2str(m_zstop);
            
             % Set mips
            cla(app.AxesX);
            imagesc(app.AxesX,reshape(max(app.angio,[],1),[app.res(2) app.res(3)]));
            set(app.AxesX,'XTickLabel','','YTickLabel','')
            colormap(app.AxesX,'gray')
            axis(app.AxesX,'equal')
            daspect(app.AxesX,[1 1 1]);
            
            cla(app.AxesY);
            imagesc(app.AxesY,reshape(max(app.angio,[],2),[app.res(1) app.res(3)]));
            set(app.AxesY,'XTickLabel','','YTickLabel','')
            colormap(app.AxesY,'gray')
            axis(app.AxesY,'equal')
            daspect(app.AxesY,[1 1 1]);
            
            cla(app.AxesZ);
            imagesc(app.AxesZ,reshape(max(app.angio,[],3),[app.res(1) app.res(2)]));
            set(app.AxesZ,'XTickLabel','','YTickLabel','')
            colormap(app.AxesZ,'gray')
            axis(app.AxesZ,'equal')
            daspect(app.AxesZ,[1 1 1]);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % add path of subfolder function
            fPath = strsplit(mfilename('fullpath'),'/');
            addpath(genpath(fullfile(fPath{1:end-1}, 'utilities')))
            drawnow;
%             app.PulseWaveVelocityToolUIFigure.WindowState = 'maximized';
        end

        % Button pushed function: LoadReconstructedDataButton
        function LoadReconstructedDataButtonPushed(app, event)
            % when load data is pushed, make sure to add the utilities
            % folder (and subfolders) so all functions are in the current
            % path
            [appDirectory,~,~] = fileparts(mfilename("fullpath"));
            addpath(genpath([appDirectory '/utilities']));
            
            clc;
            % from load data
            [app.directory, app.nframes, app.res, app.fov, app.pixdim, app.timeres, app.v, app.MAG, ...
                app.magWeightVel, app.angio, app.vMean, app.VENC, app.ori] = loadPARREC();
            
            % initialize the mask
            app.mask = ones(size(app.angio));
            
            % add to info table
            app.ScanInfoTable.Data = cat(2,cellstr([num2str(app.res(1)) ' x ' num2str(app.res(2)) ' x ' num2str(app.res(3))]),...
                cellstr([num2str(round(app.pixdim(1),1)) ' x ' num2str(round(app.pixdim(2),1)) ' x ' num2str(round(app.pixdim(3),1))]),...
                cellstr(num2str(app.timeres)),cellstr(num2str(app.nframes)));
            
            app.DataDirectoryEditField.Value = app.directory;
            m_xstart = 1; m_ystart = 1; m_zstart = 1;
            m_xstop = app.res(1); m_ystop = app.res(2); m_zstop = app.res(3);
            
            updateMIPs(app, m_xstart, m_xstop, m_ystart, m_ystop, m_zstart, m_zstop);
            
            clc;
            disp('View 3D Vasculature')
            
            normed_MIP = app.angio./max(app.angio(:));
            % fit a Gaussian to non-zero elements to determine threshold
            [muhat,sigmahat] = norm_fit(normed_MIP(:));
            
            app.segment = zeros(size(app.angio));
            app.segment(normed_MIP>muhat+2*sigmahat) = 1;
            
            app.segment = bwareaopen(app.segment,round(sum(app.segment(:)).*0.005),6); %The value at the end of the commnad in the minimum area of each segment to keep
            app.segment = imfill(app.segment,'holes'); % Fill in holes created by slow flow on the inside of vessels
            app.segment = single(app.segment);
            
            % initialize the mask
            app.mask = ones(size(app.angio));
            
            View3DSegmentation(app);
            
            % set intial parameters
            app.rotAngles = [0 0];
            app.isRawDataCropped = 0;
        end

        % Button pushed function: CleardataandrestartanalysisButton
        function CleardataandrestartanalysisButtonPushed(app, event)
            delete(app.PulseWaveVelocityToolUIFigure);  % close the app
            PulseWaveVelocityTool;                      % re-open it
        end

        % Button pushed function: LoadSegmentationDicomsButton
        function LoadSegmentationDicomsButtonPushed(app, event)
            
            % if clicked, let the user pick the directory containg the pre-segmented
            % dicoms, load them in and save them, also update 3D view
            app.segDirectory = uigetdir(app.directory, 'Select Segmentation Dicom Folder');
            app.SegmentationDirectoryEditField.Value = app.segDirectory;
            
            files = dir([app.segDirectory '/*.dcm']);
            % reset the aorta segmentation
            app.aorta_seg = zeros(size(app.angio));
            for ii = 1:numel(files)
                app.aorta_seg(:,:,ii) = dicomread([app.segDirectory '/' files(ii).name]);
            end
            app.isSegmentationLoaded = 1;
            % flip the slice dimension to match our other data
            app.aorta_seg = flip(app.aorta_seg,3);
            % normalize
            app.aorta_seg = single(app.aorta_seg/max(abs(app.aorta_seg(:))));
            app.aorta_seg(find(app.aorta_seg)) = 1;
            
            View3DSegmentation(app)
            
        end

        % Button pushed function: UnwrapVelocity
        function UnwrapVelocityButtonPushed(app, event)
            disp('Performing 4D velocity unwrapping...')
            % if raw data is not yet cropped, do it now!
            if ~app.isRawDataCropped
                [x, y, z] = ind2sub(size(app.mask),find(app.mask));
                lx = length(unique(x)); ly = length(unique(y)); lz = length(unique(z));
                maskIdx = find(app.mask);
                
                % crop velocity
                tempV = reshape(app.v,[prod(app.res),3,app.nframes]);
                tempV = tempV(maskIdx,:,:);
                tempV = reshape(tempV,lx,ly,lz,3,app.nframes);
                app.v = tempV;
                app.vMean = mean(app.v,5);
                clear tempV
                
                % crop MAG
                tempMAG = reshape(app.MAG,[prod(app.res),app.nframes]);
                tempMAG = tempMAG(maskIdx,:);
                tempMAG = reshape(tempMAG,lx,ly,lz,app.nframes);
                app.MAG = tempMAG;
                clear tempMAG;
                
                % update magWeightVel and angio
                [app.magWeightVel, app.angio] = calc_angio(app.MAG, app.v, app.VENC);
                
                % others to crop, segment and aorta_seg
                tempS = app.segment(:);
                tempS = tempS(maskIdx);
                tempS = reshape(tempS,lx,ly,lz);
                app.segment = tempS;
                clear tempS;
                if (app.isSegmentationLoaded)
                    tempS = app.aorta_seg(:);
                    tempS = tempS(maskIdx);
                    tempS = reshape(tempS,lx,ly,lz);
                    app.aorta_seg = tempS;
                    clear tempS;
                end
                app.isRawDataCropped = 1;
            end
            
            % first remove outliers (force everything to +/- VENC)
            V2 = app.v;
            V2(V2 < -app.VENC) = -app.VENC;
            V2(V2 > app.VENC) = app.VENC;
            % now scale V2 to +/- pi for unwrapping
            V2 = V2./app.VENC.*pi;
            
            % grab velocities
            phi_w_x = squeeze(V2(:,:,:,1,:));
            phi_w_y = squeeze(V2(:,:,:,2,:));
            phi_w_z = squeeze(V2(:,:,:,3,:));
            
            % perform unwrapping
            if (size(phi_w_x,4))==1
                phi_w_x_unwrapped = phi_w_x + 2*pi .* double(unwrap_3D(phi_w_x));
                phi_w_y_unwrapped = phi_w_y + 2*pi .* double(unwrap_3D(phi_w_y));
                phi_w_z_unwrapped = phi_w_z + 2*pi .* double(unwrap_3D(phi_w_z));
            else
                phi_w_x_unwrapped = phi_w_x + 2*pi .* double(unwrap_4D(phi_w_x));
                phi_w_y_unwrapped = phi_w_y + 2*pi .* double(unwrap_4D(phi_w_y));
                phi_w_z_unwrapped = phi_w_z + 2*pi .* double(unwrap_4D(phi_w_z));
            end
            
            % find the absolute maximum phi value across unwrapped data
            max_phi = max(abs(phi_w_x_unwrapped(:)));
            if max(abs(phi_w_y_unwrapped(:))) > max_phi
                max_phi = max(abs(phi_w_y_unwrapped(:)));
            end
            if max(abs(phi_w_z_unwrapped(:))) > max_phi
                max_phi = max(abs(phi_w_z_unwrapped(:)));
            end
            
            % inform of new maximum velocity
            new_venc = max_phi/pi*app.VENC;
            disp('Unwrapping finished')
            msgbox(sprintf('New maximum velocity = %4d cm/s\n',round(new_venc/10)),'Unwrapping complete');
            app.VENC = new_venc;
            
            % rescale images based on 'new' venc
            app.v(:,:,:,1,:) = phi_w_x_unwrapped./pi*app.VENC;
            app.v(:,:,:,2,:) = phi_w_y_unwrapped./pi*app.VENC;
            app.v(:,:,:,3,:) = phi_w_z_unwrapped./pi*app.VENC;
            app.vMean = mean(app.v,5);
            
        end

        % Button pushed function: CenterlineExtraction
        function CenterlineExtractionButtonPushed(app, event)
            
            % if raw data is not yet cropped, do it now!
            if ~app.isRawDataCropped
                [x, y, z] = ind2sub(size(app.mask),find(app.mask));
                lx = length(unique(x)); ly = length(unique(y)); lz = length(unique(z));
                maskIdx = find(app.mask);
                
                % crop velocity
                tempV = reshape(app.v,[prod(app.res),3,app.nframes]);
                tempV = tempV(maskIdx,:,:);
                tempV = reshape(tempV,lx,ly,lz,3,app.nframes);
                app.v = tempV;
                app.vMean = mean(app.v,5);
                clear tempV
                
                % crop MAG
                tempMAG = reshape(app.MAG,[prod(app.res),app.nframes]);
                tempMAG = tempMAG(maskIdx,:);
                tempMAG = reshape(tempMAG,lx,ly,lz,app.nframes);
                app.MAG = tempMAG;
                clear tempMAG;
                
                % update magWeightVel and angio
                [app.magWeightVel, app.angio] = calc_angio(app.MAG, app.v, app.VENC);
                
                % others to crop, segment and aorta_seg
                tempS = app.segment(:);
                tempS = tempS(maskIdx);
                tempS = reshape(tempS,lx,ly,lz);
                app.segment = tempS;
                clear tempS;
                if (app.isSegmentationLoaded)
                    tempS = app.aorta_seg(:);
                    tempS = tempS(maskIdx);
                    tempS = reshape(tempS,lx,ly,lz);
                    app.aorta_seg = tempS;
                    clear tempS;
                end
                app.isRawDataCropped = 1;
            end
            
            % these are hard-coded for now
            sortingCriteria = 3;
            spurLength = 3;
            
            se = strel('sphere',1);
            if app.isSegmentationLoaded
                ss = imerode(app.aorta_seg,se);
            else
                ss = imerode(app.segment,se);
            end
            [~,~, app.branchList, ~] = feature_extraction( ...
                sortingCriteria, spurLength, app.vMean, ss);
            
            reset3DSegmentationAndCenterline(app);
            app.TabGroup.SelectedTab = app.FlowandPulseWaveVelocityTab;
        end

        % Button pushed function: CheckcenterlinecalculateflowButton
        function CheckcenterlinecalculateflowButtonPushed(app, event)
            % Grab the branches from user input, then perform aorta segmentation, check
            % the points/segmentation is correct, and calculate flow waveforms
            
            % parse points
            str = app.BranchNumbers.Value;
            if contains(str,':')
                nums = textscan(str,'%d','Delimiter',':');
                if length(nums{1}) == 3
                    ptRange = (nums{1,1}(1)):(nums{1,1}(2)):nums{1,1}(3);
                else
                    ptRange = (nums{1,1}(1)):(nums{1,1}(2));
                end
            elseif contains(str,',')
                nums = textscan(str,'%d','Delimiter',',');
                ptRange = nums{1};
            elseif isnumeric(str2double(str))
                ptRange = str2double(str);
            else
                errordlg('Cannot parse centerline points for PWV measurement')
                return;
            end
            
            idx = [];
            % extract branches
            for b = 1:numel(ptRange)
                idx = cat(1,idx,find(app.branchList(:,4)==ptRange(b)));
            end
            app.branchActual = app.branchList(idx,1:3);
            
            reset3DSegmentationAndCenterline(app);
            hline2 = line(app.View3D_2,app.branchActual(:,2),app.branchActual(:,1),app.branchActual(:,3),...
                'Color','g','Marker','*','MarkerSize',12,'LineStyle','none');
            
            % update view angle
            rotate(hline2,[1 0 0], app.rotAngles(1))
            rotate(hline2,[0 1 0], app.rotAngles(2))
            
            choice = choosedialog;
%             checkCount = checkCount+1;
            
            if choice
                
                clc;
                % now we've found the centerline:
                % 1. calculate aorta segmentation (if not already available
                % 2. perform non-rigid registration to get time-resolved aortic
                % segmentation
                % 3. calculate flow
                
                % first smooth the centerline for better flow measures
                windowWidth = 5;    % the smoothing window
                polynomialOrder = 1;
                xsg=sgolayfilt(app.branchActual(:,1),polynomialOrder, windowWidth);
                ysg=sgolayfilt(app.branchActual(:,2),polynomialOrder, windowWidth);
                zsg=sgolayfilt(app.branchActual(:,3),polynomialOrder, windowWidth);
                app.branchActual = round([xsg,ysg,zsg]);
                
                % calculate aorta segmentation, if not already available
                if app.isSegmentationLoaded == 0   % create a new aorta_seg
                    x = app.branchActual(:,1); y = app.branchActual(:,2); z = app.branchActual(:,3);
                    index = sub2ind(size(app.segment),x,y,z);
                    g = zeros(size(app.segment));
                    g(index) = 1;
                    
                    se = strel('sphere',4);
                    gg = imdilate(g,se);
                    
                    app.aorta_seg = smooth3(gg);
                end
                
                % perform 3D nonrigid registration of segmentation and apply
                % displacement field to aorta_seg
                
                % create a mask based on the extent of aorta_seg, this helps in the
                % registration step (masking out erroneous fat signal and foldover),
                % keep up to 5 voxels on each side
                [x,y,z] = ind2sub(size(app.aorta_seg),find(app.aorta_seg));
                mask = zeros(size(app.aorta_seg));
                x1 = min(x)-5; if x1<1; x1=1; end
                x2 = max(x)+5; if x2>size(app.aorta_seg,1); x2=size(app.aorta_seg,1); end
                y1 = min(y)-5; if y1<1; y1=1; end
                y2 = max(y)+5; if y2>size(app.aorta_seg,2); y2=size(app.aorta_seg,2); end
                z1 = min(z)-5; if z1<1; z1=1; end
                z2 = max(z)+5; if z2>size(app.aorta_seg,3); z2=size(app.aorta_seg,3); end
                mask(x1:x2,y1:y2,z1:z2) = 1;
                
                % progress bar
                h = waitbar(0, sprintf('Registering aorta segmentation time frames...'));
                aortaSeg_timeResolved = zeros([size(app.angio) app.nframes]);
                fixed = mask.*app.angio; fixed = fixed/max(fixed(:));
                for j = 1:app.nframes
                    
                    if app.timeResolvedSegCheckbox.Value == 1
                        moving = app.magWeightVel(:,:,:,j);
                        moving = mask.*moving; moving = moving/max(moving(:));
                        moving = imhistmatch(moving,fixed);
                        % calculate structual similarity, if low, skip registration and use
                        % default segmentation
                        ssimval = ssim(moving,fixed);
                        if j > app.nframes/2 && j < 3/4*app.nframes
                            smoothFactor = 3.5;
                        elseif j >= 3/4*app.nframes
                            smoothFactor = 4.5;
                        else
                            smoothFactor = 2;
                        end
                        [MOVINGREG.DisplacementField,MOVINGREG.RegisteredImage] = imregdemons(moving,fixed,[100 50 25 10],'AccumulatedFieldSmoothing',smoothFactor,'PyramidLevels',4,'DisplayWaitBar',false);
                        % apply transforms directly to aortaSeg
                        aortaSeg_timeResolved(:,:,:,j) = imwarp(app.aorta_seg,-MOVINGREG.DisplacementField);
                        waitbar (j/app.nframes, h)
                    else
                        
                        aortaSeg_timeResolved(:,:,:,j) = app.aorta_seg;
                    end
                end
                close(h);
                
                % Calculate flow over whole aorta
                displayWaitBar = true;
                [app.flowPerHeartCycle_vol, app.flowPulsatile_vol, segment1, area_val] = ...
                    params_timeResolved(app.branchActual, app.angio, app.v, app.nframes, app.pixdim, aortaSeg_timeResolved, app.ori, displayWaitBar);
                
                app.branchActual = flipud(app.branchActual);
                
                app.PWVPoints.Value = ['1: ' num2str(length(app.branchActual))];
                app.PWVPointsLabel.Text = ['PWV Points (1:' num2str(length(app.branchActual)) ')'];
                
                % view the flows at each centerline point, and plot the waveforms
                viewAortaSegwFlows(app);
                plotFlowWaveforms(app);
            else
                msgbox('Change selected branch numbers')
                reset3DSegmentationAndCenterline(app);
                return;
            end
        end

        % Button pushed function: Reset3DviewButton
        function Reset3DviewButtonPushed(app, event)
            reset3DSegmentationAndCenterline(app);
        end

        % Button pushed function: PlotWaveformsButton
        function PlotWaveformsButtonPushed(app, event)
            plotFlowWaveforms(app);
        end

        % Button pushed function: CalculatePWV
        function CalculatePWVButtonPushed(app, event)
            % grab waveforms
            x = app.branchActual(:,1); y = app.branchActual(:,2); z = app.branchActual(:,3);
            index = sub2ind(size(app.segment),x,y,z);
            waveforms = app.flowPulsatile_vol(index,:);
            
            % parse points
            str = app.PWVPoints.Value;
            eval(['ptRange=[' str '];']);
            waveforms = waveforms(ptRange,:);
            
            % grab PWV calc type: 1 is cross correlation, 2 is TTF, 3
            % is Wavelet, can be updated with update to calc_pwv
            switch app.PWVType.Value
                case 'Cross-correlation'
                    PWVcalctype = 1;
                case 'Wavelet'
                    PWVcalctype = 3;
                case 'Time-to-foot'
                    PWVcalctype = 2;
            end
            
            if numel(ptRange) < 3
                errordlg('Need at least 3 points for cross-correlation PWV calculation')
                return;
            end           
            
            % pass data into calc_pwv
            branch = app.branchActual(ptRange,:);
            vox = mean(app.pixdim);
            for i=2:size(branch,1)
                dist_vec(i-1) = norm(branch(i,:)-branch(i-1,:))*vox;
            end
            % total distance along centerline
            if PWVcalctype == 2
                dist_total = [0 cumsum(dist_vec)];
            else
                dist_total = cumsum(dist_vec);
            end
            
            % the scale to interpolate by, default (1) is 1 ms
            scale = 1;
            % calculate PWV using the delay times using xcorrelation
            [D,fitObject, dist_total] = calc_pwv(waveforms,dist_total,app.timeres,PWVcalctype,scale);
            
            % the PWV, 1/slope of fit
            PWV = 1/fitObject(1);
            app.PWVDisplay.Value = num2str(round(PWV,2));
            y1 = polyval(fitObject,dist_total,'r');
            
%             % how good is the fit? calculate RMSE
%             RMSE = sqrt(mean((D-y1).^2));
%             set(hfull.RMSE_display,'string',num2str(round(RMSE,2)))
            
            % plot and display slope, PWV, RMSE
            cla(app.PWVCalcDisplay)
            scatter(app.PWVCalcDisplay,dist_total,D,'.k','SizeData',75);
            hold(app.PWVCalcDisplay,'on');
            plot(app.PWVCalcDisplay,dist_total,y1,'b','LineWidth',2);
            legend(app.PWVCalcDisplay,'delays','linear fit','Location','Northwest')
            app.PWVCalcDisplay.XLabel.String = 'distance (mm)';
            switch PWVcalctype
                case 1
                    str = 'cross-corr delay (ms)';
                    app.PWVCalcDisplay.YLim = [0 max(D)+1];
                case 2
                     str = 'time-to-foot (ms)';
                case 3
                    str = 'wavelet delay (ms)';
                    app.PWVCalcDisplay.YLim = [0 max(D)+1];
            end
            app.PWVCalcDisplay.YLabel.String = str;
            
            viewAortaSegwFlows(app);
            plotFlowWaveforms(app);
        end

        % Button pushed function: SaveResultsCallback
        function SaveResultsCallbackButtonPushed(app, event)
                        
            savePrefix = app.SaveName.Value;
            saveFolder = fullfile(app.directory, 'PWV_results'); mkdir(saveFolder);
            saveName =  fullfile(saveFolder, savePrefix);
            PWV = table(str2double(app.PWVDisplay.Value),...
                string(app.PWVPoints.Value));
            PWV.Properties.VariableNames = {'PWV','Save_Points'};
            writetable(PWV,[saveName '.xlsx']);
            
            % grab and save image
            robot = java.awt.Robot();
            temp = app.PulseWaveVelocityToolUIFigure.Position; % returns position as [left bottom width height]
            allMonPos = get(0,'MonitorPositions');
            curMon = find(temp(1)<(allMonPos(:,1)+allMonPos(:,3)),1,'first');
            curMonHeight = allMonPos(curMon,4)+1;
            pos = [temp(1) curMonHeight-(temp(2)+temp(4)) temp(3)-1 temp(4)]; % [left top width height].... UL X, UL Y, width, height
            rect = java.awt.Rectangle(pos(1),pos(2),pos(3),pos(4));
            cap = robot.createScreenCapture(rect);
            % Convert to an RGB image
            rgb = typecast(cap.getRGB(0,0,cap.getWidth,cap.getHeight,[],0,cap.getWidth),'uint8');
            imgData = zeros(cap.getHeight,cap.getWidth,3,'uint8');
            imgData(:,:,1) = reshape(rgb(3:4:end),cap.getWidth,[])';
            imgData(:,:,2) = reshape(rgb(2:4:end),cap.getWidth,[])';
            imgData(:,:,3) = reshape(rgb(1:4:end),cap.getWidth,[])';
            imwrite(imgData, [saveName '.tiff']);
            
            % save the processing information in a results struct
            results = [];
            results.directory = app.directory;
            results.segDirectory = app.segDirectory;
            results.branchActual = app.branchActual;
            results.segment = app.segment;
%             results.area_val = area_val;
            results.flowPulsatile = app.flowPulsatile_vol;
            results.flowPerHeartCycle_vol = app.flowPerHeartCycle_vol;
%             results.waveforms = waveforms;
            save([saveName '_results.mat'],'results')
            
            % inform of the saving
            msgbox(['results saved to ' saveName '.xlsx'], 'Saving complete')
        end

        % Button pushed function: CropButton
        function CropButtonPushed(app, event)
            img = app.AxesX.Children(length(app.AxesX.Children),1).CData;
            maskSz = cropImage(app,img);
            m_xstart = str2double(app.XrangeEditField.Value);
            m_xstop = str2double(app.toXEditField.Value);
            m_ystart = maskSz(2);m_ystop = maskSz(2)+maskSz(4);
            m_zstart = maskSz(1);m_zstop = maskSz(1)+maskSz(3);
            
            tempMask = zeros(size(img));tempMask(m_ystart:m_ystop,m_zstart:m_zstop) = 1;
            app.mask = app.mask.*repmat(permute(tempMask,[3 1 2]),[size(app.mask,1) 1 1]);
            
            % update angio and segmentation
            app.angio = app.angio.*app.mask;
            app.segment = app.segment.*app.mask;
            
            updateMIPs(app, m_xstart, m_xstop, m_ystart, m_ystop, m_zstart, m_zstop);
            View3DSegmentation(app);
            
        end

        % Button pushed function: CropButton_2
        function CropButton_2Pushed(app, event)
            img = app.AxesY.Children(length(app.AxesY.Children),1).CData;
            maskSz = cropImage(app,img);
            
            m_xstart = maskSz(2);m_xstop = maskSz(2)+maskSz(4);
            m_ystart = str2double(app.YrangeEditField.Value);
            m_ystop = str2double(app.toYEditField.Value);
            m_zstart = maskSz(1);m_zstop = maskSz(1)+maskSz(3);

            tempMask = zeros(size(img));tempMask(m_xstart:m_xstop,m_zstart:m_zstop) = 1;
            app.mask = app.mask.*repmat(permute(tempMask,[1 3 2]),[1 size(app.mask,2) 1]);
            
            % update angio and segmentation
            app.angio = app.angio.*app.mask;
            app.segment = app.segment.*app.mask;
            
            updateMIPs(app, m_xstart, m_xstop, m_ystart, m_ystop, m_zstart, m_zstop);
            View3DSegmentation(app);
            
        end

        % Button pushed function: CropButton_3
        function CropButton_3Pushed(app, event)
            img = app.AxesZ.Children(length(app.AxesZ.Children),1).CData;
            maskSz = cropImage(app,img);
            m_xstart = maskSz(2);m_xstop = maskSz(2)+maskSz(4);
            m_ystart = maskSz(1);m_ystop = maskSz(1)+maskSz(3);
            m_zstart = str2double(app.ZrangeEditField.Value);
            m_zstop = str2double(app.toZEditField.Value);

            tempMask = zeros(size(img));tempMask(m_xstart:m_xstop,m_ystart:m_ystop) = 1;
            app.mask = app.mask.*repmat(tempMask,[1 1 size(app.mask,3)]);
            
            % update angio and segmentation
            app.angio = app.angio.*app.mask;
            app.segment = app.segment.*app.mask;
            
            updateMIPs(app, m_xstart, m_xstop, m_ystart, m_ystop, m_zstart, m_zstop);
            View3DSegmentation(app);
        end

        % Button pushed function: RotateLeft
        function RotateLeftButtonPushed(app, event)
            rotate(app.hpatch1,[0 1 0],10)
            if (app.isSegmentationLoaded)
                rotate(app.hpatch2,[0 1 0],10)
            end
            % update rotate angles
            app.rotAngles = [app.rotAngles(1) app.rotAngles(2) + 10];
        end

        % Button pushed function: RotateRight
        function RotateRightButtonPushed(app, event)
            rotate(app.hpatch1,[0 1 0],-10)
            if (app.isSegmentationLoaded)
                rotate(app.hpatch2,[0 1 0],-10)
            end
            % update rotate angles
            app.rotAngles = [app.rotAngles(1) app.rotAngles(2) - 10];
        end

        % Button pushed function: RotateDown
        function RotateDownButtonPushed(app, event)
            rotate(app.hpatch1,[1 0 0],10)
            if (app.isSegmentationLoaded)
                rotate(app.hpatch2,[1 0 0],10)
            end
            % update rotate angles
            app.rotAngles = [app.rotAngles(1)+10 app.rotAngles(2)];
        end

        % Button pushed function: RotateUp
        function RotateUpButtonPushed(app, event)
            rotate(app.hpatch1,[1 0 0],-10)
            if (app.isSegmentationLoaded)
                rotate(app.hpatch2,[1 0 0],-10)
            end
            % update rotate angles
            app.rotAngles = [app.rotAngles(1)-10 app.rotAngles(2)];
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PulseWaveVelocityToolUIFigure and hide until all components are created
            app.PulseWaveVelocityToolUIFigure = uifigure('Visible', 'off');
            app.PulseWaveVelocityToolUIFigure.Color = [1 1 1];
            app.PulseWaveVelocityToolUIFigure.Colormap = [0.2431 0.149 0.6588;0.251 0.1647 0.7059;0.2588 0.1804 0.7529;0.2627 0.1961 0.7961;0.2706 0.2157 0.8353;0.2745 0.2353 0.8706;0.2784 0.2549 0.898;0.2784 0.2784 0.9216;0.2824 0.302 0.9412;0.2824 0.3216 0.9569;0.2784 0.3451 0.9725;0.2745 0.3686 0.9843;0.2706 0.3882 0.9922;0.2588 0.4118 0.9961;0.2431 0.4353 1;0.2196 0.4588 0.9961;0.1961 0.4863 0.9882;0.1843 0.5059 0.9804;0.1804 0.5294 0.9686;0.1765 0.549 0.9529;0.1686 0.5686 0.9373;0.1529 0.5922 0.9216;0.1451 0.6078 0.9098;0.1373 0.6275 0.898;0.1255 0.6471 0.8902;0.1098 0.6627 0.8745;0.0941 0.6784 0.8588;0.0706 0.6941 0.8392;0.0314 0.7098 0.8157;0.0039 0.7216 0.7922;0.0078 0.7294 0.7647;0.0431 0.7412 0.7412;0.098 0.749 0.7137;0.1412 0.7569 0.6824;0.1725 0.7686 0.6549;0.1922 0.7765 0.6235;0.2157 0.7843 0.5922;0.2471 0.7922 0.5569;0.2902 0.7961 0.5176;0.3412 0.8 0.4784;0.3922 0.8039 0.4353;0.4471 0.8039 0.3922;0.5059 0.8 0.349;0.5608 0.7961 0.3059;0.6157 0.7882 0.2627;0.6706 0.7804 0.2235;0.7255 0.7686 0.1922;0.7725 0.7608 0.1647;0.8196 0.749 0.1529;0.8627 0.7412 0.1608;0.902 0.7333 0.1765;0.9412 0.7294 0.2118;0.9725 0.7294 0.2392;0.9961 0.7451 0.2353;0.9961 0.7647 0.2196;0.9961 0.7882 0.2039;0.9882 0.8118 0.1882;0.9804 0.8392 0.1765;0.9686 0.8627 0.1647;0.9608 0.8902 0.1529;0.9608 0.9137 0.1412;0.9647 0.9373 0.1255;0.9686 0.9608 0.1059;0.9765 0.9843 0.0824];
            app.PulseWaveVelocityToolUIFigure.Position = [357 92 1234 760];
            app.PulseWaveVelocityToolUIFigure.Name = 'Pulse Wave Velocity Tool';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.PulseWaveVelocityToolUIFigure);
            app.TabGroup.TabLocation = 'bottom';
            app.TabGroup.Position = [1 1 1234 760];

            % Create DataLoadingandPreprocessingTab
            app.DataLoadingandPreprocessingTab = uitab(app.TabGroup);
            app.DataLoadingandPreprocessingTab.Title = 'Data Loading and Preprocessing';
            app.DataLoadingandPreprocessingTab.BackgroundColor = [1 1 1];

            % Create LoadDataPanel
            app.LoadDataPanel = uipanel(app.DataLoadingandPreprocessingTab);
            app.LoadDataPanel.BorderType = 'none';
            app.LoadDataPanel.TitlePosition = 'centertop';
            app.LoadDataPanel.Title = 'Load Data';
            app.LoadDataPanel.BackgroundColor = [1 1 1];
            app.LoadDataPanel.FontName = 'ZapfDingbats';
            app.LoadDataPanel.FontWeight = 'bold';
            app.LoadDataPanel.FontSize = 16;
            app.LoadDataPanel.Position = [1 496 617 240];

            % Create LoadReconstructedDataButton
            app.LoadReconstructedDataButton = uibutton(app.LoadDataPanel, 'push');
            app.LoadReconstructedDataButton.ButtonPushedFcn = createCallbackFcn(app, @LoadReconstructedDataButtonPushed, true);
            app.LoadReconstructedDataButton.FontName = 'ZapfDingbats';
            app.LoadReconstructedDataButton.FontSize = 16;
            app.LoadReconstructedDataButton.Position = [210 180 198 28];
            app.LoadReconstructedDataButton.Text = 'Load Reconstructed Data';

            % Create DataDirectoryEditFieldLabel
            app.DataDirectoryEditFieldLabel = uilabel(app.LoadDataPanel);
            app.DataDirectoryEditFieldLabel.HorizontalAlignment = 'right';
            app.DataDirectoryEditFieldLabel.FontName = 'ZapfDingbats';
            app.DataDirectoryEditFieldLabel.FontSize = 14;
            app.DataDirectoryEditFieldLabel.Position = [11 147 95 22];
            app.DataDirectoryEditFieldLabel.Text = 'Data Directory';

            % Create DataDirectoryEditField
            app.DataDirectoryEditField = uieditfield(app.LoadDataPanel, 'text');
            app.DataDirectoryEditField.Editable = 'off';
            app.DataDirectoryEditField.FontName = 'ZapfDingbats';
            app.DataDirectoryEditField.FontSize = 9;
            app.DataDirectoryEditField.Position = [121 147 486 22];

            % Create LoadSegmentationDicomsButton
            app.LoadSegmentationDicomsButton = uibutton(app.LoadDataPanel, 'push');
            app.LoadSegmentationDicomsButton.ButtonPushedFcn = createCallbackFcn(app, @LoadSegmentationDicomsButtonPushed, true);
            app.LoadSegmentationDicomsButton.FontName = 'ZapfDingbats';
            app.LoadSegmentationDicomsButton.FontSize = 16;
            app.LoadSegmentationDicomsButton.Position = [203 108 213 28];
            app.LoadSegmentationDicomsButton.Text = 'Load Segmentation Dicoms';

            % Create SegmentationDirectoryEditFieldLabel
            app.SegmentationDirectoryEditFieldLabel = uilabel(app.LoadDataPanel);
            app.SegmentationDirectoryEditFieldLabel.HorizontalAlignment = 'right';
            app.SegmentationDirectoryEditFieldLabel.FontName = 'ZapfDingbats';
            app.SegmentationDirectoryEditFieldLabel.FontSize = 14;
            app.SegmentationDirectoryEditFieldLabel.Position = [1 75 152 22];
            app.SegmentationDirectoryEditFieldLabel.Text = 'Segmentation Directory';

            % Create SegmentationDirectoryEditField
            app.SegmentationDirectoryEditField = uieditfield(app.LoadDataPanel, 'text');
            app.SegmentationDirectoryEditField.Editable = 'off';
            app.SegmentationDirectoryEditField.FontSize = 9;
            app.SegmentationDirectoryEditField.Position = [168 75 439 22];

            % Create ScanInfoTable
            app.ScanInfoTable = uitable(app.LoadDataPanel);
            app.ScanInfoTable.ColumnName = {'matrix'; 'resolution (mm)'; 'time resolution (ms)'; 'cardiac frames'};
            app.ScanInfoTable.RowName = {};
            app.ScanInfoTable.ColumnSortable = [false false false false];
            app.ScanInfoTable.ColumnEditable = [false false false false];
            app.ScanInfoTable.RowStriping = 'off';
            app.ScanInfoTable.FontName = 'ZapfDingbats';
            app.ScanInfoTable.Position = [25 15 570 45];

            % Create CropPanel
            app.CropPanel = uipanel(app.DataLoadingandPreprocessingTab);
            app.CropPanel.BorderType = 'none';
            app.CropPanel.TitlePosition = 'centertop';
            app.CropPanel.Title = 'Crop';
            app.CropPanel.BackgroundColor = [1 1 1];
            app.CropPanel.FontName = 'ZapfDingbats';
            app.CropPanel.FontWeight = 'bold';
            app.CropPanel.FontSize = 16;
            app.CropPanel.Position = [1 128 617 365];

            % Create XrangeEditFieldLabel
            app.XrangeEditFieldLabel = uilabel(app.CropPanel);
            app.XrangeEditFieldLabel.HorizontalAlignment = 'right';
            app.XrangeEditFieldLabel.FontName = 'ZapfDingbats';
            app.XrangeEditFieldLabel.FontSize = 14;
            app.XrangeEditFieldLabel.Position = [45 61 55 22];
            app.XrangeEditFieldLabel.Text = 'X-range';

            % Create XrangeEditField
            app.XrangeEditField = uieditfield(app.CropPanel, 'text');
            app.XrangeEditField.Editable = 'off';
            app.XrangeEditField.HorizontalAlignment = 'right';
            app.XrangeEditField.Position = [105 61 30 22];
            app.XrangeEditField.Value = '1';

            % Create toXEditFieldLabel
            app.toXEditFieldLabel = uilabel(app.CropPanel);
            app.toXEditFieldLabel.HorizontalAlignment = 'right';
            app.toXEditFieldLabel.FontName = 'ZapfDingbats';
            app.toXEditFieldLabel.FontSize = 14;
            app.toXEditFieldLabel.Position = [133 61 20 22];
            app.toXEditFieldLabel.Text = 'to';

            % Create toXEditField
            app.toXEditField = uieditfield(app.CropPanel, 'text');
            app.toXEditField.HorizontalAlignment = 'right';
            app.toXEditField.Position = [157 61 32 22];
            app.toXEditField.Value = 'res';

            % Create toYEditFieldLabel
            app.toYEditFieldLabel = uilabel(app.CropPanel);
            app.toYEditFieldLabel.HorizontalAlignment = 'right';
            app.toYEditFieldLabel.FontName = 'ZapfDingbats';
            app.toYEditFieldLabel.FontSize = 14;
            app.toYEditFieldLabel.Position = [132 30 20 22];
            app.toYEditFieldLabel.Text = 'to';

            % Create toYEditField
            app.toYEditField = uieditfield(app.CropPanel, 'text');
            app.toYEditField.HorizontalAlignment = 'right';
            app.toYEditField.Position = [156 30 32 22];
            app.toYEditField.Value = 'res';

            % Create YrangeEditFieldLabel
            app.YrangeEditFieldLabel = uilabel(app.CropPanel);
            app.YrangeEditFieldLabel.HorizontalAlignment = 'right';
            app.YrangeEditFieldLabel.FontName = 'ZapfDingbats';
            app.YrangeEditFieldLabel.FontSize = 14;
            app.YrangeEditFieldLabel.Position = [45 30 54 22];
            app.YrangeEditFieldLabel.Text = 'Y-range';

            % Create YrangeEditField
            app.YrangeEditField = uieditfield(app.CropPanel, 'text');
            app.YrangeEditField.HorizontalAlignment = 'right';
            app.YrangeEditField.FontName = 'ZapfDingbats';
            app.YrangeEditField.Position = [104 30 30 22];
            app.YrangeEditField.Value = '1';

            % Create toZEditFieldLabel
            app.toZEditFieldLabel = uilabel(app.CropPanel);
            app.toZEditFieldLabel.HorizontalAlignment = 'right';
            app.toZEditFieldLabel.FontName = 'ZapfDingbats';
            app.toZEditFieldLabel.FontSize = 14;
            app.toZEditFieldLabel.Position = [132 0 20 22];
            app.toZEditFieldLabel.Text = 'to';

            % Create toZEditField
            app.toZEditField = uieditfield(app.CropPanel, 'text');
            app.toZEditField.HorizontalAlignment = 'right';
            app.toZEditField.Position = [156 0 32 22];
            app.toZEditField.Value = 'res';

            % Create ZrangeEditFieldLabel
            app.ZrangeEditFieldLabel = uilabel(app.CropPanel);
            app.ZrangeEditFieldLabel.HorizontalAlignment = 'right';
            app.ZrangeEditFieldLabel.FontName = 'ZapfDingbats';
            app.ZrangeEditFieldLabel.FontSize = 14;
            app.ZrangeEditFieldLabel.Position = [45 0 54 22];
            app.ZrangeEditFieldLabel.Text = 'Z-range';

            % Create ZrangeEditField
            app.ZrangeEditField = uieditfield(app.CropPanel, 'text');
            app.ZrangeEditField.HorizontalAlignment = 'right';
            app.ZrangeEditField.FontName = 'ZapfDingbats';
            app.ZrangeEditField.Position = [104 0 30 22];
            app.ZrangeEditField.Value = '1';

            % Create AxesX
            app.AxesX = uiaxes(app.CropPanel);
            title(app.AxesX, '')
            xlabel(app.AxesX, '')
            ylabel(app.AxesX, '')
            app.AxesX.XColor = 'none';
            app.AxesX.XTick = 0;
            app.AxesX.YColor = 'none';
            app.AxesX.YTick = 0;
            app.AxesX.BackgroundColor = [1 1 1];
            app.AxesX.Position = [11 91 200 200];

            % Create AxesY
            app.AxesY = uiaxes(app.CropPanel);
            title(app.AxesY, '')
            xlabel(app.AxesY, '')
            ylabel(app.AxesY, '')
            app.AxesY.FontUnits = 'normalized';
            app.AxesY.FontSize = 0.0691244237613988;
            app.AxesY.XColor = 'none';
            app.AxesY.XTick = [];
            app.AxesY.YColor = 'none';
            app.AxesY.YTick = [];
            app.AxesY.BackgroundColor = [1 1 1];
            app.AxesY.Position = [212 91 200 200];

            % Create AxesZ
            app.AxesZ = uiaxes(app.CropPanel);
            title(app.AxesZ, '')
            xlabel(app.AxesZ, '')
            ylabel(app.AxesZ, '')
            app.AxesZ.XColor = 'none';
            app.AxesZ.XTick = [];
            app.AxesZ.YColor = 'none';
            app.AxesZ.YTick = [];
            app.AxesZ.BackgroundColor = [1 1 1];
            app.AxesZ.Position = [412 91 200 200];

            % Create CropButton
            app.CropButton = uibutton(app.CropPanel, 'push');
            app.CropButton.ButtonPushedFcn = createCallbackFcn(app, @CropButtonPushed, true);
            app.CropButton.IconAlignment = 'center';
            app.CropButton.FontName = 'ZapfDingbats';
            app.CropButton.FontSize = 16;
            app.CropButton.Position = [65 299 93 28];
            app.CropButton.Text = 'Crop';

            % Create CropButton_2
            app.CropButton_2 = uibutton(app.CropPanel, 'push');
            app.CropButton_2.ButtonPushedFcn = createCallbackFcn(app, @CropButton_2Pushed, true);
            app.CropButton_2.IconAlignment = 'center';
            app.CropButton_2.FontName = 'ZapfDingbats';
            app.CropButton_2.FontSize = 16;
            app.CropButton_2.Position = [266 299 93 28];
            app.CropButton_2.Text = 'Crop';

            % Create CropButton_3
            app.CropButton_3 = uibutton(app.CropPanel, 'push');
            app.CropButton_3.ButtonPushedFcn = createCallbackFcn(app, @CropButton_3Pushed, true);
            app.CropButton_3.IconAlignment = 'center';
            app.CropButton_3.FontName = 'ZapfDingbats';
            app.CropButton_3.FontSize = 16;
            app.CropButton_3.Position = [466 299 93 28];
            app.CropButton_3.Text = 'Crop';

            % Create DVisualizationPanel
            app.DVisualizationPanel = uipanel(app.DataLoadingandPreprocessingTab);
            app.DVisualizationPanel.BorderType = 'none';
            app.DVisualizationPanel.TitlePosition = 'centertop';
            app.DVisualizationPanel.Title = '3D Visualization';
            app.DVisualizationPanel.BackgroundColor = [1 1 1];
            app.DVisualizationPanel.FontName = 'ZapfDingbats';
            app.DVisualizationPanel.FontWeight = 'bold';
            app.DVisualizationPanel.FontSize = 16;
            app.DVisualizationPanel.Position = [617 6 615 730];

            % Create View3D
            app.View3D = uiaxes(app.DVisualizationPanel);
            title(app.View3D, '')
            xlabel(app.View3D, '')
            ylabel(app.View3D, '')
            app.View3D.XColor = 'none';
            app.View3D.XTick = [];
            app.View3D.YColor = 'none';
            app.View3D.YTick = [];
            app.View3D.BackgroundColor = [1 1 1];
            app.View3D.Position = [8 9 598 691];

            % Create RotateLeft
            app.RotateLeft = uibutton(app.DVisualizationPanel, 'push');
            app.RotateLeft.ButtonPushedFcn = createCallbackFcn(app, @RotateLeftButtonPushed, true);
            app.RotateLeft.IconAlignment = 'center';
            app.RotateLeft.VerticalAlignment = 'bottom';
            app.RotateLeft.BackgroundColor = [1 1 1];
            app.RotateLeft.FontName = 'ZapfDingbats';
            app.RotateLeft.FontSize = 24;
            app.RotateLeft.FontWeight = 'bold';
            app.RotateLeft.Position = [9 33 28 28];
            app.RotateLeft.Text = '<';

            % Create RotateRight
            app.RotateRight = uibutton(app.DVisualizationPanel, 'push');
            app.RotateRight.ButtonPushedFcn = createCallbackFcn(app, @RotateRightButtonPushed, true);
            app.RotateRight.IconAlignment = 'center';
            app.RotateRight.VerticalAlignment = 'bottom';
            app.RotateRight.BackgroundColor = [1 1 1];
            app.RotateRight.FontName = 'ZapfDingbats';
            app.RotateRight.FontSize = 24;
            app.RotateRight.FontWeight = 'bold';
            app.RotateRight.Position = [88 33 28 28];
            app.RotateRight.Text = '>';

            % Create Rotate
            app.Rotate = uilabel(app.DVisualizationPanel);
            app.Rotate.HorizontalAlignment = 'center';
            app.Rotate.FontName = 'ZapfDingbats';
            app.Rotate.FontSize = 16;
            app.Rotate.Position = [36 36 53 22];
            app.Rotate.Text = 'Rotate';

            % Create RotateDown
            app.RotateDown = uibutton(app.DVisualizationPanel, 'push');
            app.RotateDown.ButtonPushedFcn = createCallbackFcn(app, @RotateDownButtonPushed, true);
            app.RotateDown.IconAlignment = 'center';
            app.RotateDown.VerticalAlignment = 'bottom';
            app.RotateDown.BackgroundColor = [1 1 1];
            app.RotateDown.FontName = 'ZapfDingbats';
            app.RotateDown.FontSize = 24;
            app.RotateDown.FontWeight = 'bold';
            app.RotateDown.Position = [49 8 28 28];
            app.RotateDown.Text = 'ÿ';

            % Create RotateUp
            app.RotateUp = uibutton(app.DVisualizationPanel, 'push');
            app.RotateUp.ButtonPushedFcn = createCallbackFcn(app, @RotateUpButtonPushed, true);
            app.RotateUp.IconAlignment = 'center';
            app.RotateUp.VerticalAlignment = 'bottom';
            app.RotateUp.BackgroundColor = [1 1 1];
            app.RotateUp.FontName = 'ZapfDingbats';
            app.RotateUp.FontSize = 24;
            app.RotateUp.FontWeight = 'bold';
            app.RotateUp.Position = [49 59 28 28];
            app.RotateUp.Text = 'ÿ';

            % Create ProcessingPanel
            app.ProcessingPanel = uipanel(app.DataLoadingandPreprocessingTab);
            app.ProcessingPanel.BorderType = 'none';
            app.ProcessingPanel.TitlePosition = 'centertop';
            app.ProcessingPanel.Title = 'Processing';
            app.ProcessingPanel.BackgroundColor = [1 1 1];
            app.ProcessingPanel.FontName = 'ZapfDingbats';
            app.ProcessingPanel.FontWeight = 'bold';
            app.ProcessingPanel.FontSize = 16;
            app.ProcessingPanel.Position = [1 1 617 124];

            % Create UnwrapVelocity
            app.UnwrapVelocity = uibutton(app.ProcessingPanel, 'push');
            app.UnwrapVelocity.ButtonPushedFcn = createCallbackFcn(app, @UnwrapVelocityButtonPushed, true);
            app.UnwrapVelocity.IconAlignment = 'center';
            app.UnwrapVelocity.FontName = 'ZapfDingbats';
            app.UnwrapVelocity.FontSize = 16;
            app.UnwrapVelocity.Position = [219 59 187 28];
            app.UnwrapVelocity.Text = 'Unwrap Velocity';

            % Create CenterlineExtraction
            app.CenterlineExtraction = uibutton(app.ProcessingPanel, 'push');
            app.CenterlineExtraction.ButtonPushedFcn = createCallbackFcn(app, @CenterlineExtractionButtonPushed, true);
            app.CenterlineExtraction.IconAlignment = 'center';
            app.CenterlineExtraction.FontName = 'ZapfDingbats';
            app.CenterlineExtraction.FontSize = 16;
            app.CenterlineExtraction.Position = [219 13 187 28];
            app.CenterlineExtraction.Text = 'Centerline Extraction';

            % Create FlowandPulseWaveVelocityTab
            app.FlowandPulseWaveVelocityTab = uitab(app.TabGroup);
            app.FlowandPulseWaveVelocityTab.Title = 'Flow and Pulse Wave Velocity';
            app.FlowandPulseWaveVelocityTab.BackgroundColor = [1 1 1];

            % Create SegmentationAndCenterline
            app.SegmentationAndCenterline = uipanel(app.FlowandPulseWaveVelocityTab);
            app.SegmentationAndCenterline.BorderType = 'none';
            app.SegmentationAndCenterline.TitlePosition = 'centertop';
            app.SegmentationAndCenterline.Title = '3D View';
            app.SegmentationAndCenterline.BackgroundColor = [1 1 1];
            app.SegmentationAndCenterline.FontName = 'ZapfDingbats';
            app.SegmentationAndCenterline.FontWeight = 'bold';
            app.SegmentationAndCenterline.FontSize = 16;
            app.SegmentationAndCenterline.Position = [1 18 450 702];

            % Create View3D_2
            app.View3D_2 = uiaxes(app.SegmentationAndCenterline);
            title(app.View3D_2, '')
            xlabel(app.View3D_2, '')
            ylabel(app.View3D_2, '')
            app.View3D_2.XColor = 'none';
            app.View3D_2.XTick = [];
            app.View3D_2.YColor = 'none';
            app.View3D_2.YTick = [];
            app.View3D_2.BackgroundColor = [1 1 1];
            app.View3D_2.Position = [1 0 459 669];

            % Create Reset3DviewButton
            app.Reset3DviewButton = uibutton(app.SegmentationAndCenterline, 'push');
            app.Reset3DviewButton.ButtonPushedFcn = createCallbackFcn(app, @Reset3DviewButtonPushed, true);
            app.Reset3DviewButton.FontName = 'ZapfDingbats';
            app.Reset3DviewButton.FontSize = 14;
            app.Reset3DviewButton.Position = [343 4 108 29];
            app.Reset3DviewButton.Text = 'Reset 3D view';

            % Create BranchNumberTitle
            app.BranchNumberTitle = uilabel(app.FlowandPulseWaveVelocityTab);
            app.BranchNumberTitle.HorizontalAlignment = 'right';
            app.BranchNumberTitle.FontName = 'ZapfDingbats';
            app.BranchNumberTitle.FontSize = 18;
            app.BranchNumberTitle.FontWeight = 'bold';
            app.BranchNumberTitle.Position = [548 696 198 22];
            app.BranchNumberTitle.Text = 'Set branches for aorta';

            % Create BranchNumberLabel
            app.BranchNumberLabel = uilabel(app.FlowandPulseWaveVelocityTab);
            app.BranchNumberLabel.HorizontalAlignment = 'right';
            app.BranchNumberLabel.FontName = 'ZapfDingbats';
            app.BranchNumberLabel.FontSize = 18;
            app.BranchNumberLabel.Position = [518 665 160 22];
            app.BranchNumberLabel.Text = 'Branch number(s): ';

            % Create BranchNumbers
            app.BranchNumbers = uieditfield(app.FlowandPulseWaveVelocityTab, 'text');
            app.BranchNumbers.FontName = 'ZapfDingbats';
            app.BranchNumbers.FontSize = 16;
            app.BranchNumbers.Position = [677 664 70 23];

            % Create CheckcenterlinecalculateflowButton
            app.CheckcenterlinecalculateflowButton = uibutton(app.FlowandPulseWaveVelocityTab, 'push');
            app.CheckcenterlinecalculateflowButton.ButtonPushedFcn = createCallbackFcn(app, @CheckcenterlinecalculateflowButtonPushed, true);
            app.CheckcenterlinecalculateflowButton.FontName = 'ZapfDingbats';
            app.CheckcenterlinecalculateflowButton.FontSize = 18;
            app.CheckcenterlinecalculateflowButton.Position = [511 622 271 29];
            app.CheckcenterlinecalculateflowButton.Text = 'Check centerline, calculate flow';

            % Create FlipcenterlinepointsButton
            app.FlipcenterlinepointsButton = uibutton(app.FlowandPulseWaveVelocityTab, 'push');
            app.FlipcenterlinepointsButton.FontName = 'ZapfDingbats';
            app.FlipcenterlinepointsButton.FontSize = 18;
            app.FlipcenterlinepointsButton.Visible = 'off';
            app.FlipcenterlinepointsButton.Position = [933 679 213 29];
            app.FlipcenterlinepointsButton.Text = 'Flip centerline points';

            % Create PWVPointsTitle
            app.PWVPointsTitle = uilabel(app.FlowandPulseWaveVelocityTab);
            app.PWVPointsTitle.HorizontalAlignment = 'right';
            app.PWVPointsTitle.FontName = 'ZapfDingbats';
            app.PWVPointsTitle.FontSize = 18;
            app.PWVPointsTitle.FontWeight = 'bold';
            app.PWVPointsTitle.Position = [866 556 264 22];
            app.PWVPointsTitle.Text = 'Set PWV measurement points';

            % Create PWVPointsLabel
            app.PWVPointsLabel = uilabel(app.FlowandPulseWaveVelocityTab);
            app.PWVPointsLabel.HorizontalAlignment = 'right';
            app.PWVPointsLabel.FontName = 'ZapfDingbats';
            app.PWVPointsLabel.FontSize = 18;
            app.PWVPointsLabel.Position = [863 522 160 22];
            app.PWVPointsLabel.Text = 'PWV points: ';

            % Create PWVPoints
            app.PWVPoints = uieditfield(app.FlowandPulseWaveVelocityTab, 'text');
            app.PWVPoints.FontName = 'ZapfDingbats';
            app.PWVPoints.FontSize = 16;
            app.PWVPoints.Tooltip = {'Select the point labels in the centerline that will be used for PWV calculation. Minimum of 3 points needed.'};
            app.PWVPoints.Position = [1027 522 106 23];

            % Create PlotWaveformsButton
            app.PlotWaveformsButton = uibutton(app.FlowandPulseWaveVelocityTab, 'push');
            app.PlotWaveformsButton.ButtonPushedFcn = createCallbackFcn(app, @PlotWaveformsButtonPushed, true);
            app.PlotWaveformsButton.FontName = 'ZapfDingbats';
            app.PlotWaveformsButton.FontSize = 18;
            app.PlotWaveformsButton.Position = [539 522 216 29];
            app.PlotWaveformsButton.Text = 'Examine flow waveforms';

            % Create WaveformsDisplay
            app.WaveformsDisplay = uiaxes(app.FlowandPulseWaveVelocityTab);
            title(app.WaveformsDisplay, '')
            xlabel(app.WaveformsDisplay, 'Cardiac time (ms)')
            ylabel(app.WaveformsDisplay, 'Flow (mL/s)')
            app.WaveformsDisplay.FontName = 'ZapfDingbats';
            app.WaveformsDisplay.FontSize = 14;
            app.WaveformsDisplay.BackgroundColor = [1 1 1];
            app.WaveformsDisplay.Position = [470 269 711 232];

            % Create CalculatePWV
            app.CalculatePWV = uibutton(app.FlowandPulseWaveVelocityTab, 'push');
            app.CalculatePWV.ButtonPushedFcn = createCallbackFcn(app, @CalculatePWVButtonPushed, true);
            app.CalculatePWV.FontSize = 18;
            app.CalculatePWV.Position = [690 216 137 29];
            app.CalculatePWV.Text = 'Calculate PWV';

            % Create PWVType
            app.PWVType = uidropdown(app.FlowandPulseWaveVelocityTab);
            app.PWVType.Items = {'Cross-correlation', 'Wavelet', 'Time-to-foot'};
            app.PWVType.FontName = 'ZapfDingbats';
            app.PWVType.FontSize = 14;
            app.PWVType.Position = [540 219 140 22];
            app.PWVType.Value = 'Cross-correlation';

            % Create PWVDisplayTitle
            app.PWVDisplayTitle = uilabel(app.FlowandPulseWaveVelocityTab);
            app.PWVDisplayTitle.HorizontalAlignment = 'right';
            app.PWVDisplayTitle.FontName = 'ZapfDingbats';
            app.PWVDisplayTitle.FontSize = 18;
            app.PWVDisplayTitle.FontWeight = 'bold';
            app.PWVDisplayTitle.Position = [944 197 192 22];
            app.PWVDisplayTitle.Text = 'Calculated PWV (m/s)';

            % Create PWVDisplay
            app.PWVDisplay = uieditfield(app.FlowandPulseWaveVelocityTab, 'text');
            app.PWVDisplay.Editable = 'off';
            app.PWVDisplay.FontName = 'ZapfDingbats';
            app.PWVDisplay.FontSize = 18;
            app.PWVDisplay.Tooltip = {'Select the point labels in the centerline that will be used for PWV calculation. Minimum of 2 points needed.'};
            app.PWVDisplay.Position = [1001 163 80 26];

            % Create PWVCalcDisplay
            app.PWVCalcDisplay = uiaxes(app.FlowandPulseWaveVelocityTab);
            title(app.PWVCalcDisplay, '')
            xlabel(app.PWVCalcDisplay, 'delay (ms)')
            ylabel(app.PWVCalcDisplay, 'distance (mm)')
            app.PWVCalcDisplay.FontName = 'ZapfDingbats';
            app.PWVCalcDisplay.FontSize = 14;
            app.PWVCalcDisplay.BackgroundColor = [1 1 1];
            app.PWVCalcDisplay.Position = [491 25 385 190];

            % Create SavingTitle
            app.SavingTitle = uilabel(app.FlowandPulseWaveVelocityTab);
            app.SavingTitle.HorizontalAlignment = 'right';
            app.SavingTitle.FontName = 'ZapfDingbats';
            app.SavingTitle.FontSize = 18;
            app.SavingTitle.FontWeight = 'bold';
            app.SavingTitle.Position = [1008 97 65 22];
            app.SavingTitle.Text = 'Saving';

            % Create SaveResultsCallback
            app.SaveResultsCallback = uibutton(app.FlowandPulseWaveVelocityTab, 'push');
            app.SaveResultsCallback.ButtonPushedFcn = createCallbackFcn(app, @SaveResultsCallbackButtonPushed, true);
            app.SaveResultsCallback.FontSize = 18;
            app.SaveResultsCallback.Position = [991 25 100 29];
            app.SaveResultsCallback.Text = 'Save';

            % Create SaveName
            app.SaveName = uidropdown(app.FlowandPulseWaveVelocityTab);
            app.SaveName.Items = {'Global PWV', ' AAo PWV', ' DAo PWV'};
            app.SaveName.FontName = 'ZapfDingbats';
            app.SaveName.FontSize = 14;
            app.SaveName.Position = [971 66 140 22];
            app.SaveName.Value = 'Global PWV';

            % Create timeResolvedSegCheckbox
            app.timeResolvedSegCheckbox = uicheckbox(app.FlowandPulseWaveVelocityTab);
            app.timeResolvedSegCheckbox.Tooltip = {'If checked'; ' non-rigid registration is performed for each cardiac time frame to get time-resolved segmentations. Takes time!'};
            app.timeResolvedSegCheckbox.Text = 'Time-resolved segmentation';
            app.timeResolvedSegCheckbox.FontSize = 14;
            app.timeResolvedSegCheckbox.Position = [806 628 199 22];

            % Create ResetWorkSpace
            app.ResetWorkSpace = uitab(app.TabGroup);
            app.ResetWorkSpace.Title = 'Reset Workspace';

            % Create CleardataandrestartanalysisButton
            app.CleardataandrestartanalysisButton = uibutton(app.ResetWorkSpace, 'push');
            app.CleardataandrestartanalysisButton.ButtonPushedFcn = createCallbackFcn(app, @CleardataandrestartanalysisButtonPushed, true);
            app.CleardataandrestartanalysisButton.Position = [101 550 178 22];
            app.CleardataandrestartanalysisButton.Text = 'Clear data and restart analysis';

            % Show the figure after all components are created
            app.PulseWaveVelocityToolUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = PulseWaveVelocityTool_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PulseWaveVelocityToolUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PulseWaveVelocityToolUIFigure)
        end
    end
end