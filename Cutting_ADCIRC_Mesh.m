% This script can be used to custom cut an ADCIRC mesh
% and generate a new for.14 
%                               By Amin Kiaghadi
%                          Last Update: May 02, 2020

tic
clear all;
clc;
close all;


%% Clipping option

% Please note Flag=1 means using Shapefile and for that your MATLAB should
% have Mapping Tool installed
% Defualt value for Flag is 0 to use coordinate box

Flag=1;

if Flag==1

    addpath('C:\Kiaghadi\UT\ADCIRC_Tools\IO_Files\Shapefiles');
    Shape=shaperead('Lower_Neches_Area.shp');
end

% SETTING THE UPPER AND LOWER BOUNDS FOR X and Y

% please input the coordinates of the box you want to obtain data
Xupper=-93.25;
Yupper=30.5;
Xlower=-94.25;
Ylower=29.8;


%% R E A D   A D C I R C   fort.14


fort_14_file = 'Latest_For14_Amin.grd';

% Reading nodes and elements
fort_14_id = fopen( fort_14_file, 'r');
First_Line = fgetl ( fort_14_id );           % first line- no need
Second_Line = fscanf(fort_14_id, '%d %d ', 2);  % second line- summary
ne      = Second_Line(1);     %number of elements
nn      = Second_Line(2);     %number of nodes
%
nodes_n_x_y_z = fscanf(fort_14_id, '%lf', [4,nn]);
% n  = nodes_n_x_y_z(1,:);     % node number
% x  = nodes_n_x_y_z(2,:);     % node x
% y  = nodes_n_x_y_z(3,:);     % node y
% z = nodes_n_x_y_z(4,:);      % node z
Nodes=nodes_n_x_y_z';
cells_14 = fscanf(fort_14_id, '%ld', [5, ne])';
cells_14 = cells_14(:, 1:5);
%
fclose(fort_14_id);



% Reading boundary conditions

fort_14_id = fopen( fort_14_file, 'r');
txt = fileread(fort_14_file);
TextAsCells = regexp(txt, '\n', 'split');
mask = ~cellfun(@isempty, strfind(TextAsCells,...
    '= Number of node'));
All_Bound_loc_beg=find(mask);      % starting line of all boundaries


% Reading open boundary nodes

OpenB_loc_beg=All_Bound_loc_beg(1);      % starting line of Open Boundary

for i=1:OpenB_loc_beg-1
    dummy=fgetl (fort_14_id);
end


% fseek(fort_14_id,605965795,'bof')
Num_OpenB=fscanf(fort_14_id, '%d\n', 1)';
dummy=fgetl (fort_14_id);
OpenB_Nodes_Original = fscanf(fort_14_id, '%d\n', [1 Num_OpenB])';
% cursor_End_Open_Bound=ftell(fort_14_id);
Line_number=ne+nn+2+3+Num_OpenB;

%  fseek(fort_14_id,605967741,'bof')


% Reading land boundaries
Boundary=cell(size(All_Bound_loc_beg,2),1);
Boundary{1}=OpenB_Nodes_Original;
%  i=2end fseek(fort_14_id,605999379,'bof')
num_extra_lines=All_Bound_loc_beg(2)-Line_number-1;

   for j=1:num_extra_lines
         dummy=fgetl (fort_14_id); % moving extra line
   end

   
Num_Bound=zeros(size(All_Bound_loc_beg,2),1);
Category=zeros(size(All_Bound_loc_beg,2),1);
Num_Bound(1)=Num_OpenB;   
for i=2:size(All_Bound_loc_beg,2)
    

   num_cat_land=fscanf(fort_14_id, '%d', 2)';
   Category(i)=num_cat_land(2);
   Num_Bound(i)=num_cat_land(1);
   
    switch Category(i)
      case 24
          dummy=fgetl (fort_14_id);
          Boundary{i} = fscanf(fort_14_id, '%lf', [5 Num_Bound(i)])';
        case 23
          dummy=fgetl (fort_14_id);
          Boundary{i} = fscanf(fort_14_id, '%lf', [3 Num_Bound(i)])';
      otherwise
          dummy=fgetl (fort_14_id);
          Boundary{i} = fscanf(fort_14_id, '%d\n', [1 Num_Bound(i)])';
    end
   Line_number=Line_number+Num_Bound(i)+num_extra_lines+2;
end

fclose(fort_14_id);

%% TRUNCATING THE BUILT MATRIX 
if Flag==0
    Truncated_Nodes_order=find(Nodes(:,2)> Xlower & Nodes(:,2)< Xupper &...
                          Nodes(:,3)> Ylower & Nodes(:,3)< Yupper);
    Truncated_Nodes=Nodes(Truncated_Nodes_order,:);
    
else
    Xupper=max(Shape.X);
    Yupper=max(Shape.Y);
    Xlower=min(Shape.X);
    Ylower=min(Shape.Y);
   
    Truncated_Nodes_order=find(Nodes(:,2)> Xlower & Nodes(:,2)< Xupper &...
                          Nodes(:,3)> Ylower & Nodes(:,3)< Yupper);
    Truncated_Nodes=Nodes(Truncated_Nodes_order,:);
 
    for i=1:size(Truncated_Nodes,1)
            in(i)=inpolygon(Truncated_Nodes(i,2),Truncated_Nodes(i,3)...
                  ,Shape.X,Shape.Y);
    end
    Loc_tr_temp=find (in)';
    Truncated_Nodes=Truncated_Nodes(Loc_tr_temp,:);  
end
    
    

%% Switching node numbers in the connectivity matrix

% New Node Numbers
New_Node_Number=1:size(Truncated_Nodes,1);
Numbers=[Truncated_Nodes(:,1),New_Node_Number'];

% Finding truncated nodes in the connectivity matrix
[val1,pos1]=intersect(cells_14(:,3),Truncated_Nodes(:,1));
indices1 = find(ismember(cells_14(:,3),val1,'rows'));
[val2,pos2]=intersect(cells_14(:,4),Truncated_Nodes(:,1));
indices2 = find(ismember(cells_14(:,4),val2,'rows'));
[val3,pos3]=intersect(cells_14(:,5),Truncated_Nodes(:,1));
indices3 = find(ismember(cells_14(:,5),val3,'rows'));
Position=[indices1;indices2;indices3];
Position=unique(Position);     % Cell # for of cells with at least one node within the box
Cells_14_Temp=cells_14(Position,1:5);   % Getting the connectivity matrix info

% Switching Node Numbers
results1=zeros(size(Cells_14_Temp,1),1);
results2=zeros(size(Cells_14_Temp,1),1);
results3=zeros(size(Cells_14_Temp,1),1);
[tf,rowWithElement] = ismember(Cells_14_Temp(:,3),Numbers(:,1));
results1(tf) = Numbers(rowWithElement(tf),2);
[tf,rowWithElement] = ismember(Cells_14_Temp(:,4),Numbers(:,1));
results2(tf) = Numbers(rowWithElement(tf),2);
[tf,rowWithElement] = ismember(Cells_14_Temp(:,5),Numbers(:,1));
results3(tf) = Numbers(rowWithElement(tf),2);
Results=[results1,results2,results3];

% Removing elements that are partially within the new domain
for i=1:3
Results_with_Zero = find(~Results(:,i));
Results(Results_with_Zero,:)=[];
end


%% Switching node numbers in  boundaries

% total_weir_number=0;

Boundary_updated = Boundary;
Boundary_count=zeros(size(Boundary,1),1);
for i=1:size(Boundary,1)
   
    switch Category(i)
      case 24
        Boundary_first=Boundary{i}(:,1:1);
        Boundary_second=Boundary{i}(:,1:2);
        [tf,rowWithElement] = ismember(Boundary_first(:,1),Numbers(:,1));
        Boundary_res_first=zeros(size(tf));
        Boundary_res_first(tf) = Numbers(rowWithElement(tf),2);
        [tf,rowWithElement] = ismember(Boundary_second(:,2),Numbers(:,1));
        Boundary_res_second=zeros(size(tf));
        Boundary_res_second(tf) = Numbers(rowWithElement(tf),2);
           if sum(tf)~=0
              Boundary_count(i)=1;
           end
        Boundary_updated{i}(:,1)=Boundary_res_first;
        Boundary_updated{i}(:,2)=Boundary_res_second;
        clear Boundary_res_first
        clear Boundary_res_second
        
      otherwise
        [tf,rowWithElement] = ismember(Boundary{i}(:,1:1),Numbers(:,1));
        Temp=zeros(size(tf));
        Temp(tf) = Numbers(rowWithElement(tf),2); 
           if sum(tf)~=0
              Boundary_count(i)=1;
           end
        Boundary_updated{i}=Temp;
        clear Temp
    end
end
       
Boundary_location=find(Boundary_count);        
            
Boundary_Final= cell(size(Boundary_location,1),1);

Final_number=zeros(size(Boundary_location,1),1);

for i=1:size(Boundary_location,1)

    Boundary_temp=Boundary_updated{Boundary_location(i)};
    Boundary_temp=Boundary_temp(all(Boundary_temp,2),:);
    Boundary_Final{i}=Boundary_temp;
    Final_number(i)=size(Boundary_temp,1);

end

Final_Category=Category (Boundary_location)';


total_Boundary_number=0;


for i=1:size(Final_Category,2)
    switch Final_Category(i)
        case 24
             temp=2*Final_number(i);
        case 0
             temp=0;
        otherwise
             temp=Final_number(i);
    end
    total_Boundary_number=total_Boundary_number+temp;
end

%% Preparing final matrices for writing


New_Nodes=[New_Node_Number',Truncated_Nodes(:,2:4)];

New_Cell_Number=1:size(Results,1);

Cell_Constant= ones(size(Results,1),1)+2;

New_Cells=[New_Cell_Number',Cell_Constant,Results];

Dimensions=[size(New_Cells,1), size(New_Nodes,1)];

  

%% Writing

fid = fopen('C:\Kiaghadi\UT\ADCIRC_Tools\MATLAB_Neches_fort.14', 'wt');
fprintf(fid, '%s\n','Modified fort.14');

% Nodes

fprintf(fid, '%d %d\n', Dimensions);
for i=1:size(New_Nodes,1)
fprintf(fid, '%9d %16.10f %16.10f %16.10f\n', New_Nodes(i,:));
end


% Connectivity

for i=1:size(New_Cells,1)
fprintf(fid, '%9d %5d % d % d % d\n', New_Cells(i,:));
end


% Boundaries

%Open Boundary 
if Boundary_location(1)==1
    fprintf(fid, '%s %s\n','1',' = Number of open boundaries');
    fprintf(fid, '%d %s\n',size(Boundary_updated{1},1)...
            ,' = Total number of open boundary nodes');
    fprintf(fid, '%d %s\n',size(Boundary_updated{1},1)...
            ,' = Number of nodes for open boundary 1');
     for j=1:size(Boundary_Final{1},1)
           fprintf(fid, '%d\n', Boundary_updated{1}(j));
     end   
else
    fprintf(fid, ' %s\n','0 = Number of open boundaries');
    fprintf(fid, ' %s\n','0 = Total number of open boundary nodes');  
 
end

fprintf(fid, '%d %s\n',size(Boundary_Final,1),...
            ' = Number of land boundaries');
fprintf(fid, '%d %s\n',total_Boundary_number,...
            ' = Total number of land boundary nodes');        
    
    
for i=1:size(Final_Category,2)
    i
    switch Final_Category(i)
        case 0
            continue
        case 24
             fprintf(fid, '%d',Final_number(i));
             fprintf(fid, '%s %d %s\n',...
              ' 24 = Number of node pairs for weir (land boundary ',i,')');
             for j=1:Final_number(i)
                fprintf(fid, '%d %d %f %d %d\n', Boundary_Final{i}(j,:));
             end
          
        case 23
             fprintf(fid, '%d',Final_number(i));
             fprintf(fid, '%s %d\n'...
              ,' 23 = Number of nodes for land boundary ',i);
             for j=1:Final_number(i)
                fprintf(fid, '%d %d  %d\n', Boundary_Final{i}(j,:));
             end            
          
        otherwise
             fprintf(fid, '%d',Final_number(i));
             fprintf(fid, '%s %d\n'...
              ,' 20 = Number of nodes for land boundary ',i);
             for j=1:Final_number(i)
                fprintf(fid, '%d\n', Boundary_Final{i}(j));   
             end
    end    
end
    
   
fclose(fid);

disp('Done ! Hopefully successfully !');
toc
load handel
sound(y,Fs)