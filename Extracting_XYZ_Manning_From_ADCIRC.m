% This script extract xyz and manning coeffients from fort.13 & 14
%                               By Amin Kiaghadi
%                          Last Update: March 30, 2020



tic
clear variables;
clc;
close all;

%% SETTING THE UPPER AND LOWER BOUNDS FOR X and Y

% please input the coordinates of the box you want to obtain data
Xupper=-93.25;
Yupper=30.5;
Xlower=-94.25;
Ylower=29.8;


%% R E A D   A D C I R C   fort.14


fort_14_file = 'IO_Files/fort_S2G.14';
fort_13_file = 'IO_Files/newfile.txt';

fort_14_id = fopen( fort_14_file, 'r');
First_Line = fgetl ( fort_14_id );           % first line- no need
Second_Line = fscanf(fort_14_id, '%d %d ', 2);  % second line- summary
ne      = Second_Line(1);     %number of elements
nn      = Second_Line(2);     %number of nodes
%
nodes_n_x_y_z = fscanf(fort_14_id, '%lf', [4,nn]);
n  = nodes_n_x_y_z(1,:);     % node number
x  = nodes_n_x_y_z(2,:);     % node x
y  = nodes_n_x_y_z(3,:);     % node y
z = nodes_n_x_y_z(4,:);      % node z
%
fclose(fort_14_id);

%% R E A D   A D C I R C   fort.13

fort_13_id = fopen( fort_13_file, 'r');

%finding the location of manning coeffients in fort.13

text = fileread(fort_13_file);
TextAsCells = regexp(text, '\n', 'split');
mask = ~cellfun(@isempty, strfind(TextAsCells,...
    'mannings_n_at_sea_floor'));
Mann_loc_beg=find(mask,2);      % starting line of mannings in fort.13

mask1 = ~cellfun(@isempty, strfind(TextAsCells, ...
    'surface_submergence_state'));
Mann_loc_end=find(mask1,2);      % end line of manning defualt in fort.13

%getting rid of the lines up to manning deafult value
for i=1:Mann_loc_end(1)-2
    dummy=fgetl (fort_13_id);
end

%getting manning defult value

Manning_diff = fscanf(fort_13_id, '%f ', 1);  % Manning defualt value

% reseting the location of cursor 
fseek(fort_13_id,0,'bof')


%getting rid of the lines up to number of manning values
for i=1:Mann_loc_beg(2)
    dummy=fgetl (fort_13_id);
end



%getting number and values of manning coefficients
Manning_count = fscanf(fort_13_id, '%d ', 1);  % number of coefficients
Manning_info = fscanf(fort_13_id, '%lf', [2,Manning_count]); % manning values

fclose(fort_13_id);

Manning_node= Manning_info(1,:);          % node numbers with manning
Manning_values= Manning_info(2,:);        % manning coefficients

% Associating manning coefficients to each ADCIRC node
manning=ones(nn,1)*Manning_diff;         % setting all to defult
manning(Manning_node)=Manning_values;    % assing available coefficients

%% BUILDING A MATRIX WITH ALL NODE NUMBER, XYZ, and MANNING 

info=[n', x', y', z', manning];


%% TRUNCATING THE BUILT MATRIX EITH THE CUSTOMIZED CORRDINATION BOX

Truncated_info_order=find(info(:,2)> Xlower & info(:,2)< Xupper &...
    info(:,3)> Ylower & info(:,3)< Yupper);
Truncated_info=info(Truncated_info_order,:);


dlmwrite('ADCIRC_TRUNCATED.csv', Truncated_info, 'delimiter', ',', 'precision', 10);


disp('Done ! Hopefully successfully !');
toc
load handel
sound(y,Fs)