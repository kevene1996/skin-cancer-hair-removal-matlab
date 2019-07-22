% ImageMask function detects and computes the edges(or hair particles of 
% a given image 'Image'
% An colored RGB image is formed by 3 matrixes each contain their RGB
% colors
% REQUIREMENTS FOR PROPER REMOVAL:
%  - Image should be clear with visible hair for a better removal
%  - hair particles should be black 
%  - Images with a lot of hair might get distorted
%  - hair particles inside tumors might not get removed if tumor is same
%  color as hair
function Imagemask = ImageMask(Image) 
RImage = Image(:,:,1); % get the Red matrix pixel values from image (1 is the 3D index number representing the red channel)
GImage = Image(:,:,2); % get the Green matrix pixel values from image (2 is the 3D index number representing the green channel)
BImage = Image(:,:,3); % get the Blue matrix pixel values from image(3 is the 3D index number representing the blue channel)

% Compute threshold for each color channel using different edge detection methods
% Each method is used depending on the hair characteristics in the image.
method = 'sobel'; % method can be 'sobel' or 'canny'
[~,threshold1] = edge(RImage,method); 
[~,threshold2] = edge(GImage,method);
[~,threshold3] = edge(BImage,method);

% Apply the edge detector to every RGB matrix using specified method with
% respective threshold, edge() function ignores all edges that are not
% stronger than the threshold
IedgeR = edge(RImage,method,threshold1);
IedgeG = edge(GImage,method,threshold2);
IedgeB = edge(BImage,method,threshold3);
 
% create morphological structure used to dilate the RGB images, the
% structure is a line that ranges between 0 to 90 degrees and has a certain thickness
% used to increase the size of the detected edge.
thickness = 5; % change thickness to improve image mask, varies between 1 and 10.
se90 = strel('line',thickness,90);
se0 = strel('line',thickness,0);

% image dilation is applied to every RGB color channel in order to make the edges 
% thicker
ImagedilateR = imdilate(IedgeR,[se90 se0]);
ImagedilateG = imdilate(IedgeG,[se90 se0]);
ImagedilateB = imdilate(IedgeB,[se90 se0]);
Imagedilate = ImagedilateR&ImagedilateG&ImagedilateB; % convert the 3 images into 1 image    
Imagecomp = imcomplement(Imagedilate); % compliment image so white pixels of 'hair' become black
Imagemask  = double(mat2gray(im2double(Imagecomp)) == 1); % convert binary image to grayscale double
% first binary image is converted to double, then to grayscale and to double again

% makes sure that mask and image are of the same size, repmat() replicates
% array of image acording to the specifies size ( 3 in this case)z
if size(Imagemask,3)==1 && size(Image,3)>1 
    Imagemask = repmat(Imagemask,[1,1,size(Image,3)]); 
end
end
