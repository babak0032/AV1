% Demo to burn lines the user draws into an image.  
% Multiple lines may be drawn and they all get burned in.

%----- Initializing steps -----
% Clean up
clc;
clear all;
close all;
workspace; % Display the workspace panel.
fontSize = 20;

hasIPT = license('test', 'image_toolbox');
if ~hasIPT
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
end

% Display images to prepare for the demo.
monochromeImage = imread('pout.tif');
subplot(2, 2, 1);
imshow(monochromeImage);
title('Original Image', 'FontSize', fontSize);
subplot(2, 2, 2);
imshow(monochromeImage);
title('DRAW LINE HERE!!!', 'FontSize', fontSize);
subplot(2, 2, 4);
imshow(monochromeImage);
title('Original Image with lines burned into image', 'FontSize', fontSize);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
set(gcf,'name','Image Analysis Demo','numbertitle','off') 

%----- Burn line into image -----
burnedImage = imread('pout.tif');
% Create a binary image for all the lines we will draw.
cumulativeBinaryImage = false(size(burnedImage));
subplot(2, 2, 3);
imshow(monochromeImage);
title('Binary Image', 'FontSize', fontSize);
% Create line mask, h, as an ROI object over the second image in the bottom row.
axis on;
again = true;
lineCount = 0;
while again && lineCount < 20
	promptMessage = sprintf('Draw line #%d in the upper right image,\nor Quit?', lineCount + 1);
	titleBarCaption = 'Continue?';
	button = questdlg(promptMessage, titleBarCaption, 'Draw', 'Quit', 'Draw');
	if strcmpi(button, 'Quit')
		break;
	end
	lineCount = lineCount + 1;
	subplot(2, 2, 2);
	hLine = imline(gca);
	caption = sprintf('DRAW HERE.  Original Image with %d lines in overlay.', lineCount);
	title(caption, 'FontSize', fontSize);

	% Create a binary image ("mask") from the ROI object.
	singleLineBinaryImage = hLine.createMask();
	% OR it in to the "all lines" binary image mask we're building up.
	cumulativeBinaryImage = cumulativeBinaryImage | singleLineBinaryImage;
	% Display the lines mask.
	subplot(2, 2, 3);
	imshow(cumulativeBinaryImage);
	caption = sprintf('Binary mask of the %d lines', lineCount);
	title(caption, 'FontSize', fontSize);
	
	% Burn line into image by setting it to 255 wherever the mask is true.
	burnedImage(cumulativeBinaryImage) = 255;
	% Display the image with the "burned in" line.
	subplot(2, 2, 4);
	cla;
	imshow(burnedImage);
	caption = sprintf('New image with %d lines burned into image', lineCount);
	title(caption, 'FontSize', fontSize);
end

