function [] = groundTrackPlot2(lon, lat, col, line)

% plot for ground tracks (second layer, to be used after groundTrackPlot)
%
% usage:
% [] = groundTrackPlot2(lon, lat, col, line)
%
% Input arguments:
% ----------------------------------------------------------------
% lon               [Nx1]       longitude vector                [deg]
% lat               [Nx1]       latitude vector                 [deg]
% col               [string]    color for the plot              [-]
% line              [1x1]       lineWidth                       [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% N/A
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

lon_range = [-180, 180];
lat_range = [-90, 90];

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
        plot(lon(first:last),lat(first:last), LineWidth=line, Color=col, HandleVisibility="off")
        first = i + 1;
    end
end
plot(lon(first:end),lat(first:end), LineWidth=line, Color=col, HandleVisibility="off")

plot(lon(end),lat(end), "square", Color="black", LineWidth=2)

axis equal
axis([lon_range, lat_range])
axis on

xlabel("Latitude [deg]")
ylabel("Longitude [deg]")

leg=legend("Start", "End 1st Plot", 'End 2nd Plot', Location="best");
fontsize(leg, 10, "points")

hold off

end

