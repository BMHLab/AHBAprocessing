%CMAP Create a custom colormap from RGB values
% 
% FUNCTION: cmap_out = cmap(in_arg,in_size,in_cutd,in_cutl)
%
% Create custom colormaps from any number of colors
% Inputs: 1) in_arg = RGB triplets: [A B C] with values between 0-1
%         2) in_size (optional): Length of colorbar (i.e. 10 = 10 colors)
%         3) in_cutd/l (optional): Cutoff values D,L (%): Remove D% darkest and L% lightest colors
%
% If only color is given, default cutoff = 10%
% If only color and length are given, X% = Y%
%
% If CMAP_OUT is given, the function returns a matrix.
% If CMAP_OUT is not given, colormap is set as specified by input arguments.
%
% NOTE: Please download the excellent RGB.m code from Kristjn Jnasson, University of Iceland.
% Colors may then be specified as 'steelblue', 'forestgreen' etc. in stead of RGB triplets.
% See: http://www.mathworks.com/matlabcentral/fileexchange/24497-rgb-triple-of-color-name-version-2
%
% Any of the following expressions give the same result:
%   cmap('green');
%   a = cmap('green'); colormap(a)
%   colormap(cmap('green'));
%
% Syntax: 1) optional_variable = cmap([1 0 0]); (Red tones)
%         2) optional_variable = cmap('red'); (Type rgb chart to see available color names)
%         3) optional_variable = cmap([1 0 0],12); (Red tones, 12 colors from black to white)
%         4) optional_variable = cmap('red',12,15) (Red tones, 12 colors, with the 15% darkest & lightest colors removed)
%         5) optional_variable = cmap('red',12,30,0) (Red tones, 12 colors, with only the 30% darkest colors removed)
%
% % Example 1: Red colors 
% figure; pcolor([0:20;0:20]); cmap('red'); colorbar('horiz'); title('Ex 1: Red, using rgb.m function')
% figure; pcolor([0:20;0:20]'); cmap([0 1 0]); colorbar('horiz'); title('Ex 1: Green, using RGB values [0 1 0]')
%   
% % Example 2: 10 blue colors
% figure; pcolor([0:20;0:20]); colormap(cmap('steelblue',10));
% shading flat; colorbar('horiz'); title('Ex 2: Steel blue, 10 colors')
%
% % Example 3: dark blue - dark red colorbar, 20 colors, 5% darkest values and 20% lightest values cut
% figure; pcolor([0:20;0:20]);
% colormap([cmap('steelblue',10,5,20);flipud(cmap('firebrick',10,5,20))]); colorbar('horiz');
% title('Ex 3: Blue -> Red, 20 colors, 5% darkest and 20% lightest colors are cut')
%
% % Example 4: 2/3 green + 1/3 yellow, individually cut
% p = rand(20,20)*100;
% figure; pcolor(p);
% colormap([cmap('seagreen',ceil(length(unique(p))-1)/3*2,20,1);flipud(cmap('yellow',floor(length(unique(p))-1)/3,20,8))]);
% colorbar('horiz'); title('Ex 4: 2/3 green, 1/3 yellow colors, darkest/lightest colors cut individually')
%
% % Example 5: Mix a new spectrum with different colors and tones:
% p = rand(20,20)*100; figure; pcolor(p)
% colmap1 = cmap('steelblue',10,10,15);
% colmap2 = cmap('orange',6,40,20);
% colmap3 = rgb('orangered');
% colmap4 = cmap('red',8,15,45);
% colormap([colmap1;flipud(colmap2);colmap3;flipud(colmap4)]);
% colorbar('horiz'); title('Ex 5: Spectrum')
%
% ABOUT CMAP
% This program is public domain and may be distributed freely.
% Author: Erik Kvaleberg, Norwegian Naval Training Establishment (ekvaleberg@gmail.com)
% July 2013
%---------------


function [cmap_out] = make_cmap(in_arg,in_size,in_cutd,in_cutl)

 if ischar(in_arg)
  in_arg = lower(in_arg);
  cvec = rgb(in_arg);
 elseif isnumeric(in_arg)
  if sum(size(in_arg))==4 && min(size(in_arg))==1
   cvec = in_arg;
  else
   fprintf(1,'\n\n\n Error: Invalid color \n\n');
   return
  end
 end

if nargin == 1
 cinc = 0.01;
else
 cinc = 0.0001;
end

if nargin == 3
 if in_cutd <= 0; in_cutd = 1; end
else if nargin == 4
 if in_cutd <= 0; in_cutd = 1; end
 if in_cutl <= 0; in_cutl = 1; end
end; end

map = NaN(length(-max(cvec)-1:cinc:max(cvec)+1),3);
for cii = 1:length(-max(cvec)-1:cinc:max(cvec)+1)
 map(cii,:) = (cvec-1-max(cvec)) + cinc*cii;
end
clear cii

map(find(map<0)) = 0; map(find(map>1)) = 1;

for cii = 1:length(map(:,1))
 if (map(cii,1)==1) & (map(cii,2)==1) & (map(cii,3)==1) %#ok<AND2>
  endind(cii) = cii;
 else
  endind(cii) = NaN;
 end
end
clear cii

for cii = 1:length(map(:,1))
if (map(cii,1)==0) & (map(cii,2)==0) & (map(cii,3)==0) %#ok<AND2>
  startind(cii) = cii;
 else
  startind(cii) = NaN;
 end
end
clear cii

map = map(max(startind)+1:min(endind)-1,:);

if nargout == 0
 if nargin == 2
  % Cut the 10% darkest and 10% lightest colors if in_size,in_cutxy is not given
  map = map(ceil(length(map(:,1))*(10/100)):end - floor(length(map(:,1))*(10/100)),:);
  % Divide into Z colors
  colormap(map(floor(length(map(:,1))/in_size):floor(length(map(:,1))/in_size):length(map(:,1)),:));
 else if nargin == 3
  % Cut the X% darkest and X% lightest colors if in_size is given
  map = map(ceil(length(map(:,1))*(in_cutd/100)):end - floor(length(map(:,1))*(in_cutd/100)),:);
  % Divide into Z colors
  colormap(map(length(map(:,1))/in_size:floor(length(map(:,1))/in_size):length(map(:,1)),:));
 else if nargin == 4
  % Cut the X% darkest and Y% lightest colors if in_size is given
  map = map(ceil(length(map(:,1))*(in_cutd/100)):end - floor(length(map(:,1))*(in_cutl/100)),:);
  % Divide into Z colors
  colormap(map(length(map(:,1))/in_size:floor(length(map(:,1))/in_size):length(map(:,1)),:));
 else
  % Cut the 10% darkest and 10% lightest colors if in_size,in_cutxy is not given
  colormap(map(floor(length(map(:,1))*(10/100)):end - floor(length(map(:,1))*(10/100)),:));
 end
 end
 end

else

if nargin == 2
 % Cut the 10% darkest and 10% lightest colors if in_size,in_cutxy is not given
 map = map(ceil(length(map(:,1))*(10/100)):end - floor(length(map(:,1))*(10/100)),:);
 % Divide into Z colors
 cmap_out = map(floor(length(map(:,1))/in_size):floor(length(map(:,1))/in_size):length(map(:,1)),:);
else if nargin == 3
 % Cut the X% darkest and X% lightest colors if in_size is given
 map = map(ceil(length(map(:,1))*(in_cutd/100)):end - floor(length(map(:,1))*(in_cutd/100)),:);
 % Divide into Z colors
 cmap_out = map(length(map(:,1))/in_size:floor(length(map(:,1))/in_size):length(map(:,1)),:);
else if nargin == 4
 % Cut the X% darkest and Y% lightest colors if in_size is given
 map = map(ceil(length(map(:,1))*(in_cutd/100)):end - floor(length(map(:,1))*(in_cutl/100)),:);
 % Divide into Z colors
 cmap_out = map(floor(length(map(:,1))/in_size):floor(length(map(:,1))/in_size):length(map(:,1)),:);
else
 % Cut the 10% darkest and 10% lightest colors if in_size,in_cutxy is not given
 cmap_out = map(floor(length(map(:,1))*(10/100)):end - floor(length(map(:,1))*(10/100)),:);
end
end
end

end


