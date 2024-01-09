function s = jupiterPlot

% Generates a plot of Jupiter and returns it's graphic handler
%
% Usage
% s = jupiterPlot
%
% Output arguments:
% ---------------------------------------------------------------------
% s             [1x1]   Jupiter - obj                [surface/-]

% Earth settings
[x1,y1,z1] = sphere(50);
mult = 69911;                     % Jupiter radius
s = surface(x1*mult,y1*mult,z1*mult);

load topo 
s.FaceColor = 'texturemap';
s.CData = topo;
s.EdgeColor = 'none';          % remove edges
s.FaceLighting = 'gouraud';    % lighting for specific surfaces
s.SpecularStrength = 0.4;      % reflected light
axis equal

% image_file = 'https://svs.gsfc.nasa.gov/vis/a000000/a003600/a003615/flat_earth03.jpg';
image_file = "flat_jupiter.jpg";
cdata = imread(image_file);
cdata=cdata([end:-1:1],:,:);
set(s, 'FaceColor', 'texturemap', 'CData', cdata, 'EdgeColor', 'none');

hold on
axis equal, grid on

end

% Earth’s flat image used for 3d representations and animations was downloaded from Nasa website
% (https://svs.gsfc.nasa.gov/vis/a010000/a012000/a012021/Hubble_Jupiter_maps_compared_H264_1920x1080_print.jpg) as there's no restrictions under 
% "Still Images, Audio Recordings, Video, and Related Computer Files for Non-Commercial Use" 
% (https://www.nasa.gov/multimedia/guidelines/index.html) for our use. 
% We do not own this image: NASA/Goddard Space Flight Center Scientific Visualization Studio 
% The Blue Marble Next Generation data is courtesy of Reto Stockli (NASA/GSFC) 
% and NASA's Earth Observatory.