
clear all; close all;

main = '/Users/gino/Dropbox (Brain Conquistadores)/OPERATING_ROOM/SPEECH_STUDY/NO_ARREST/14.-35414460/%s';


inDIR = sprintf(main,'J_C_data');
inNoN = sprintf(main,'NoN');
outCI = sprintf(main,'CI');
inDTI = sprintf(main,'DTI_tracks');

NR_name_sort_module = 'NoN_nodes_mod.txt';
outlinks = 'matrix_outlinks.txt';
% dti_name = 'dti_table.txt';

th1 = 0.;
th2 = 0;

cd(inNoN); % move to module file folder
voxel = importdata(NR_name_sort_module); % file which contains the nodes, sorted by module


size_mod = zeros(max(voxel(:,6))+1,2); % vector to store: label cluster / numb of nodes in cluster


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    COUNT ELEMENTS IN EACH MODULE   %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% create an array which contains: label of cluster / numb of nodes in that
% cluster, i.e. size

mod = 1; % module 1, first module
for i = 1:length(voxel) % from 1 to max number of nodes in NoN
    if voxel(i,6) == mod
        size_mod(mod,1) = mod;
        size_mod(mod,2) = size_mod(mod,2) + 1;
    else
        mod = mod + 1;
        size_mod(mod,1) = mod;
        size_mod(mod,2) = size_mod(mod,2) + 1; 
    end
end

mat_links = importdata(outlinks); % import outlinks matrix 
N = size(mat_links,1); % size matrix
A = zeros(N,N); 

% A IS THE FMRI MATRIX NORMALIZED BY SIZE OF THE MODULES and Thresholded
% in A the thresholding is done on the already normalized matrix
for i = 1:N
    for j = (i+1):N 
        
        link = mat_links(i,j)/(size_mod(i,2)+size_mod(j,2))/2;
        if link > th1
            A(i,j) = link;
            A(j,i) = link;
        end
        
    end
end


% in math_th the thresholding is done in the non normalized matrix
mat_th = mat_links;
mat_th(mat_links < th2) = 0;

myLabel = cell(length(A));

%myLabel = {'MFG(L)','BA(L)','Motor(L)','Tongue(L)','PreMotor(L)','SupM(L)','SMA','BA(R)','MFG(R)','Tongue(R)','SupM(R) 1','SupM(R) 2','IOA(R)'}';
myLabel = {'1 - BA(L)','2 - WA(L)','3 - PreM(L) 1','4 - PreM(L) 2','5 - MFG(L) 1','6 - MFG(L) 2','7 - IOA(L)','8 - SMA','9 - BA(R)','10 - MFG(R)'};


% cd(inDTI);
% dti = importdata(dti_name);
% 
% dti2 = zeros(size(dti,1),size(dti,1));
% dti3 = dti2; 
% 
% % DTI simmetrized but not normalized by the modules size
% for i = 1:N
%     for j = (i+1):N         
%             dti2(i,j) = (dti(i,j)+dti(j,i))/2; 
%             dti2(j,i) = dti2(i,j); 
%     end
% end
% 
% 
% % DTI simmetrized and normalized by the modules size
% for i = 1:N
%     for j = (i+1):N         
%             dti3(i,j) = dti2(i,j)/(size_mod(i,2)+size_mod(j,2)); 
%             dti3(j,i) = dti3(i,j); 
%     end
% end

% DTI simmetrized but not normalized by the modules size
% FIG3 = figure(3)
% circularGraph(dti2,'Label',myLabel);

% DTI simmetrized and normalized by the modules size


FIG1 = figure(1)
circularGraph(A,'Label',myLabel);  % FMRI matrix, normalized by the size of the modules and thresholded
saveas(FIG1,'circular_fMRI_1.png');

FIG2 = figure(2)
circularGraph(mat_th,'Label',myLabel); % FMRI matrix of links, normalized by max value and thresholded
saveas(FIG2,'circular_fMRI_2.png');


% FIG3 = figure(4)
% circularGraph(dti3,'Label',myLabel); % dti matrix, normalized by the size of modules
% saveas(FIG3,'circular_DTI.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           GraphViz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% fMRI 1 %%%%%%%%

fileID = fopen('graph_fMRI_1.txt','w');
fprintf(fileID,'graph {\n');

for i = 1:N
    for j = (i+1):N
        if A(i,j)~=0
            line = sprintf('"%s" -- "%s" [penwidth = %0.4f]\n',myLabel{i},myLabel{j},A(i,j));
            fprintf(fileID,line);
        end
    end
    if A(i,:) == 0
        line = sprintf('"%s"\n',myLabel{i});
        fprintf(fileID,line);
    end
end
fprintf(fileID,'}\n');
fclose(fileID);



%%% fMRI 2 %%%%%%%%

fileID = fopen('graph_fMRI_2.txt','w');
fprintf(fileID,'graph {\n');

mat_th(:) = mat_th(:)./max(mat_th(:));

for i = 1:N
    for j = (i+1):N
        if mat_th(i,j)~=0
            line = sprintf('"%s" -- "%s" [penwidth = %0.4f]\n',myLabel{i},myLabel{j},mat_th(i,j));
            fprintf(fileID,line);
        end
    end
    if mat_th(i,:) == 0
        line = sprintf('"%s"\n',myLabel{i});
        fprintf(fileID,line);
    end
end
fprintf(fileID,'}\n');
fclose(fileID);

% %%% DTI %%%%%%%%


% fileID = fopen('graph_DTI.txt','w');
% fprintf(fileID,'graph {\n');
% 
% for i = 1:N
%     for j = (i+1):N
%         if dti3(i,j)~=0
%             line = sprintf('"%s" -- "%s" [penwidth = %0.4f]\n',myLabel{i},myLabel{j},dti3(i,j));
%             fprintf(fileID,line);
%         end
%     end
%     if dti3(i,:) == 0
%         line = sprintf('"%s"\n',myLabel{i});
%         fprintf(fileID,line);
%     end
% end
% fprintf(fileID,'}\n');
% fclose(fileID);





