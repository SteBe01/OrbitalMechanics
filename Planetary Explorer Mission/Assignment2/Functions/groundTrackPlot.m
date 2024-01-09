function [] = groundTrackPlot(lon, lat, EarthPlot_name)

% plot for ground tracks
%
% usage:
% groundTrackPlot(lon, lat, EarthPlot_name)
%
% Input arguments:
% ----------------------------------------------------------------
% lon               [3x1]       longitude vector                [deg]
% lat               [3x1]       latitude vector                 [deg]
% EarthPlot_name    [string]    name of the Earth Plot (img)    [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% N/A
% 
% CONTRIBUITORS:
% Pier Francesco A. Bachini
% Stefano Belleti
% Chiara Giardini
% Carolina Gómez Sánchez

earth_img = imread(EarthPlot_name);

lon_range = [-180, 180];
lat_range = [-90, 90];

figure()
imshow(earth_img, 'XData', lon_range, 'YData', lat_range);
hold on;

len = length(lon);
first = 1;
split = 0;

for i=1:len-1
    if split == 1
        split = 0;
    end

    if abs(lon(i+1))+abs(lon(i)) > 300
        if lon(i+1)*lon(i) < 0
            split = 1;
        end
    end
    if abs(lat(i+1))+abs(lat(i)) > 150
        if lat(i+1)*lat(i) < 0
            split = 1;
        end
    end

    if split
        last = i;
        plot(lon(first:last),lat(first:last), LineWidth=1.5, Color="green", HandleVisibility="off")
        first = i + 1;
    end
end
plot(lon(first:end),lat(first:end), LineWidth=1.5, Color="green", HandleVisibility="off")

plot(lon(1),lat(1), "or",Color="magenta", LineWidth=2)
plot(lon(end),lat(end), "square", Color="yellow", LineWidth=2)

axis equal
axis([lon_range, lat_range])
axis on

xlabel("Latitude [deg]")
ylabel("Longitude [deg]")
legend("Start", "End", Location="best")

hold off

end
