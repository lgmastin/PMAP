function [Z,P,T,RH] = MetReader(filename)
% reads meteorological data files from either of the following sources:
%   1) University of Wyoming - https://weather.uwyo.edu/upperair/sounding.html
%   2) NOAA - https://www.ncdc.noaa.gov/cdo-web/

fileID = fopen(filename);

% loop through lines of meteorological until the recognized header is found
found_header = false;
while ~found_header
    str = fgetl(fileID);
    if numel(str)>=10 && strcmp(str(1:10),'----------')
        format = 'Wyoming';
        found_header = true;
    elseif numel(str)>=5 && strcmp(str(1:5),'PRESS')
        format = 'NOAA';
        found_header = true;
    elseif str==-1
        error('Met file not recognized: Check if header is formatted propertly')
    end
end    
    
if strcmp(format,'Wyoming')
    C = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f %f','HeaderLines',3);
    Z = C{2}; % elevation, m
    P = C{1}*100; % pressure, Pa
    T = C{3}+273.15; % temperature, K
    RH = C{5}; % relative humidity, %
elseif strcmp(format,'NOAA')
    for i = 1:4
        str = fgetl(fileID);
    end
    C1 = textscan(str,'%f %fE %f %f %f %f');
    C2 = textscan(fileID,'%f %f %f %f %f %f');
    Z = [C1{2}; C2{2}]; % elevation, m
    P = [C1{1}; C2{1}]*100; % pressure, Pa
    T = [C1{3}; C2{3}]+273.15; % temperature, K
    TD = [C1{4}; C2{4}]+273.15; % dew point temperature, K
    for i = 1:numel(TD)
        RH(i) = 100*get_sat_partial_pressure(TD(i))/get_sat_partial_pressure(T(i)); % relative humidity, %
    end
end
fclose(fileID);

end