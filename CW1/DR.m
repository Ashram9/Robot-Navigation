% Car Dead Reckoning All-Epoch
addpath('Coursework 1 Data/');
addpath('Coursework 1 Software for Students to Use/');
% Define_Constants
Define_Constants

% load .csv data file
data_DR = readmatrix('Dead_reckoning.csv');
times = data_DR(:, 1); % Time in seconds
speed = mean(data_DR(:, 2:5), 2); % Forward average speed in m/s
heading_deg = data_DR(:, 7); % Heading in degrees
heading_rad = heading_deg * deg_to_rad; % Heading in radians

data_GNSS = readmatrix('results_GNSS_KF_with_OD.csv');
initial_latitude_deg = data_GNSS(1, 2); % Initial latitude in degrees
initial_longitude_deg = data_GNSS(1, 3); % Initial longitude in degrees


% (1) Resolve the average velocity in North and East direction
% store average velocities
avg_velocities = zeros(length(times), 2);
% store damped instantaneous velocities
damped_inst_velocities = zeros(length(times), 2);

% Initialisation
avg_velocities(1, :) = [cos(heading_rad(1)), sin(heading_rad(1))] * speed(1);
damped_inst_velocities(1, :) = avg_velocities(1, :);

for k = 2:length(times)
    heading_cur = heading_rad(k);
    heading_prev = heading_rad(k-1);

    % Resolving Vector
    M_NE = 1/2 * [cos(heading_cur) + cos(heading_prev), sin(heading_cur) + sin(heading_prev)];
    avg_velocities(k, :) = M_NE * speed(k);

    damped_inst_velocities(k, :) = 1.7 * avg_velocities(k, :) - 0.7 * damped_inst_velocities(k-1, :);
end


% (2) Calculate the latitude and longitude at epoch k
% Define variables to store latitude and longitude in radians for calculating purpose
latitudes_rad = zeros(length(times), 1);
longitudes_rad = zeros(length(times), 1);

% Initialisation
latitude_rad = initial_latitude_deg * deg_to_rad; % Convert initial latitude to radians
longitude_rad = initial_longitude_deg * deg_to_rad; % Convert initial longitude to radians
latitudes_rad(1) = latitude_rad;
longitudes_rad(1) = longitude_rad;

for k = 2:length(times)
    delta_time = times(k) - times(k-1); % Time difference between epochs k and k-1
    [RN, RE] = Radii_of_curvature(latitude_rad); % Compute meridian and transverse radii of curvature
    height = data_GNSS(k, 4);

    velocity_N = avg_velocities(k, 1);
    velocity_E = avg_velocities(k, 2);

    % Calculate latitude change using the formula
    delta_latitude = velocity_N * delta_time / (RN + height);
    % Update latitude
    latitude_rad = latitude_rad + delta_latitude;
    % Calculate longitude change using the formula
    delta_longitude = velocity_E * delta_time / ((RE + height) * cos(latitude_rad));
    % Update longitude
    longitude_rad = longitude_rad + delta_longitude;
    
    % Store the new latitude and longitude
    latitudes_rad(k) = latitude_rad;
    longitudes_rad(k) = longitude_rad;
end

% Convert the final latitude and longitude to degrees
latitudes_deg = latitudes_rad * rad_to_deg;
longitudes_deg = longitudes_rad * rad_to_deg;

%%
% Define the filename for the CSV file
filename = 'results_DR.csv';
% Open the CSV file for writing
fileID = fopen(filename, 'w');
% Check if the file was opened successfully
if fileID == -1
    error('Failed to open the file for writing');
end

% Define the header as a cell array of strings
header = {'Time(s)', 'Latitude(deg)', 'Longitude(deg)', 'Damped Inst North(m/s)', 'Damped Inst East(m/s)'};
% Write the header to the CSV file
fprintf(fileID, '%s,', header{1:end-1}); % Write all header elements except the last one
fprintf(fileID, '%s\n', header{end});     % Write the last header element and a newline

% Write the data to the CSV file
for k = 1:length(times)
    time = times(k);
    data_DR = [time, latitudes_deg(k), longitudes_deg(k), damped_inst_velocities(k, :)];
    
    % Write the data to the CSV file
    fprintf(fileID, '%f,', data_DR(1:end-1)); % Write all data elements except the last one
    fprintf(fileID, '%f\n', data_DR(end));     % Write the last data element and a newline
end

% Close the CSV file
fclose(fileID);
disp(['Results saved to "', filename, '"']);