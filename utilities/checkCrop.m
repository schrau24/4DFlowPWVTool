function choice = checkCrop
answer = questdlg('Is this correct?', ...
	'Check cropped image', ...
	'Yes, continue','Cancel Cropping', 'No','Yes, continue');
% Handle response
switch answer
    case 'Yes, continue'
        choice = 1;
    case 'No'
        choice = 0;
    case 'Cancel Cropping'
        choice = 2;
end