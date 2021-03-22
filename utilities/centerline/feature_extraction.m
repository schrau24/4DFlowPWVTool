function [CL,branchMat, branchList, branchTextList] = feature_extraction( ...
    sortingCriteria, spurLength, vMean, segment)


% Thinning operation
Y = Skeleton3D(segment);
disp('Finished thinning operation')

% VASCULAR TREE CONSTUCTION
% specify sortingCriteria as either
% = 2 to get all branches connected to each other sorting (few branches)
% = 3 to get branch by branch sorting (many branches)
[CL, branchMat, branchList, branchTextList, junctionMat, junctionList] = ...
    centerline(Y, vMean, spurLength, sortingCriteria);
disp('Feature Extraction finished')

end