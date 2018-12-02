
% This code computes teh Jij by thresholding the Cij with a lambda
% threshold. Prints the Cij matrix and its figure for each lamdba value

clear all; close all;

main = '/Users/gino/Dropbox (Brain Conquistadores)/OPERATING_ROOM/SPEECH_STUDY/NO_ARREST/23.-35414460/%s';

inDIR = sprintf(main,'Brain_Maps');
outDIR_data = sprintf(main,'J_C_data');
outDIR_Fig = sprintf(main,'J_C_Figures');
inDIR_NoN = sprintf(main,'NoN');

%time_series_file = 'ts_th_0.5473.txt';
time_series_file = 'time_series.txt';


lambdaList = 0.6:0.01:1;
Delta = 0.02;


NR_name_sort_module = 'NoN_nodes_mod.txt';

cd(inDIR); % move to module file folder
voxel = importdata(NR_name_sort_module); % file which contains the nodes, sorted by module: NR / x / y / z / corr / lab_module

size_mod = zeros(max(voxel(:,6))+1,2); % vector to store: label cluster / numb of nodes in cluster


% create an array which contains: label of cluster / numb of nodes in that
% cluster, i.e. size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               MODULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module = 1; % module 1, first module
for i = 1:length(voxel) % from 1 to max number of nodes in NoN
    if voxel(i,6) == module
        size_mod(module,1) = module;
        size_mod(module,2) = size_mod(module,2) + 1;
    else
        module = module + 1;
        size_mod(module,1) = module;
        size_mod(module,2) = size_mod(module,2) + 1;
    end
end

M = sum(size_mod(:,2)); % size of the final matrix J

% Move to input directory
cd(inDIR);

% import data of the correlation matrix from time series
data = importdata(time_series_file);
ts_T = transpose(data);

Cij_orig = corr(ts_T); % correlation matrix
%  dlmwrite('Cij.txt',Cij_orig,'delimiter','\t');

% Cij_orig = importdata('Cij.txt'); %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           MATRIX COLORMAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minCor = min(Cij_orig(:));
maxCor = max(Cij_orig(:));

delta = (maxCor+abs(minCor))/256;
Zero_min = round(abs(minCor)/delta); % zero from the bottom on the colorbar
Zero_max = round(maxCor/delta); % zero from the top of the colorbar
step = 1; % space to be left for white color (times two)
up = Zero_min - step;
down = Zero_max - step;
if minCor < 0
    
    % redColorMap = [linspace(0.3, 1,up); zeros(1,up); zeros(1,up)]'; % set red color map
    % greenColorMap = [ zeros(1,down); linspace(1, 0.3, down); zeros(1,down)]';  %set green color map
    redColorMap =  autumn(up);
    greenColorMap = flipud(winter(down)); % above zero color map
    cmap = [redColorMap;white(2*step); greenColorMap]; % set full color map -> red / white / green
else
    % greenColorMap = [ zeros(1,256); linspace(1, 0., 256); zeros(1,256)]';  %set green color map
    greenColorMap = flipud(winter(down)); % above zero color map
    cmap =  greenColorMap;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FIGURE CORRELATION MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the correlation Matrix
blocks = figure;
fig_C = imagesc(Cij_orig);
% title('Cij Matrix - original', 'FontSize', 15); % set title

xticks(0:200:length(Cij_orig));
yticks(0:100:length(Cij_orig));


ax=gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
ax.GridAlpha=0.8;
ax.MinorGridAlpha = 1;
%grid on
%grid minor
%set(gca,'xtick')

%lambdaList = transpose(importdata('lambda_C_values.txt'));

%%%%%%%%  BLOCKS   %%%%%%%%
hold on
orig = 0;
edge = size_mod(1,2);

for i = 1:(length(size_mod)-1)
    
    h(i) = rectangle('position',[orig orig edge edge]);
    set(h(i),'EdgeColor',[0 0 0],'linewidth',2);
    orig = orig + size_mod(i,2);
    edge = size_mod(i+1,2);
    
end

hold on
orig = 0;
edge = size_mod(1,2);

for i = 1:(length(size_mod)-1)
    
    h(i) = rectangle('position', [0 orig edge size(Cij_orig,1)]);
    set(h(i),'EdgeColor',[0 0 0],'linewidth',1);
    orig = orig + size_mod(i,2);
    edge = edge + size_mod(i+1,2);
    
end

hold on
orig = 0;
edge = size_mod(1,2);

for i = 1:(length(size_mod)-1)
    
    h(i) = rectangle('position', [orig 0 size(Cij_orig,1) edge]);
    set(h(i),'EdgeColor',[0 0 0],'linewidth',1);
    orig = orig + size_mod(i,2);
    edge = edge + size_mod(i+1,2);
    
end


colormap(cmap);
colormap(jet);
colorbar;
set(gca, 'FontSize', 25,'Fontname','Garamond')
cd(inDIR_NoN);

saveas(fig_C,'Cij_original.png');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLD CORRELATION MATRIX %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

percolation = [];
for lambda = lambdaList
    
    Cij = Cij_orig; % correlation matrix
    
    % threshold Cij
    for i=1:length(Cij)
        for j=1:length(Cij)
            if abs(Cij(i,j)) < lambda
                Cij(i,j) = 0;
            end
        end
    end
    
    % compute size of the GC at this lambda and store the values: lambda, GC size
    GC = size(largestcomponent(Cij),2)/size(Cij_orig,1);
    percolation = [percolation; lambda, GC];   % lambda, GC size
    
    
    % Move to output directory and write matrix Cij
    cd(outDIR_data)
    name_out_C = sprintf('C_%0.4f.txt',lambda);
    dlmwrite(name_out_C,Cij,'delimiter','\t');
    
    
    % find the value zero on the color map scale (colorbar)
    % for the Connectivity Matrix

    if mod(lambda,Delta) == 0
        
        minCor = min(Cij(:))
        maxCor = max(Cij(:));
        delta = (maxCor+abs(minCor))/256;
        if minCor < 0
            delta = (maxCor+abs(minCor))/256;
            if delta == 0
                delta = 1;
            end
            
            step = round(2*lambda/delta); % space to be left for white color (times two)
            
            up = abs(round((minCor+lambda)/delta));
            down = abs(round((maxCor-lambda)/delta));
            
            % redColorMap = [linspace(0.3, 1,up); zeros(1,up); zeros(1,up)]'; % set red color map
            % greenColorMap = [ zeros(1,down); linspace(1, 0.3, down); zeros(1,down)]';  %set green color map
            redColorMap =  autumn(up);
            greenColorMap = flipud(winter(down)); % above zero color map
%             greenColorMap = jet; % above zero color map
            cmap = [redColorMap; white(step); greenColorMap]; % set full color map -> red / white / green
        else
            step = round(lambda/delta); % space to be left for white color (times two)
            down = abs(round((maxCor-lambda)*256));
            % greenColorMap = [ zeros(1,down); linspace(1, 0.3, down); zeros(1,down)]';  %set green color map
            greenColorMap = flipud(winter(down)); % above zero color map
            cmap = [white(step); greenColorMap]; % set full color map -> red / white / green
            
        end
        
        %         Move to output directory for figures
        cd(outDIR_Fig)
        
        figure(1);
        figC = imagesc(Cij)
        namefigC = sprintf('C_%0.4f.jpg',lambda)
        name_title = sprintf('Cij - lambda = %0.4f',lambda);
        title(name_title, 'FontSize', 15); % set title
%         xticks(0:100:length(Cij));
        yticks(0:100:length(Cij));
        ax=gca;
        ax.XAxis.FontSize = 12;
        ax.YAxis.FontSize = 12;
        ax.GridAlpha=1;
        ax.MinorGridAlpha = 1;
        %grid on;
        %grid minor
        
        %%%%%%%  BLOCKS   %%%%%%%%
        hold on
        orig = 0;
        edge = size_mod(1,2);
        
        for i = 1:(length(size_mod)-1)
            
            h(i) = rectangle('position',[orig orig edge edge]);
            set(h(i),'EdgeColor',[0 0 0],'linewidth',2);
            orig = orig + size_mod(i,2);
            edge = size_mod(i+1,2);
            
        end
        
        hold on
        orig = 0;
        edge = size_mod(1,2);
        
        for i = 1:(length(size_mod)-1)
            
            h(i) = rectangle('position', [0 orig edge size(Cij_orig,1)]);
            set(h(i),'EdgeColor',[0 0 0],'linewidth',1.25);
            orig = orig + size_mod(i,2);
            edge = edge + size_mod(i+1,2);
            
        end
        
        hold on
        orig = 0;
        edge = size_mod(1,2);
        
        for i = 1:(length(size_mod)-1)
            
            h(i) = rectangle('position', [orig 0 size(Cij_orig,1) edge]);
            set(h(i),'EdgeColor',[0 0 0],'linewidth',1.25);
            orig = orig + size_mod(i,2);
            edge = edge + size_mod(i+1,2);
            
        end
        
        colormap(cmap);  % set the colorscheme
        colorbar;
        set(gca, 'FontSize', 25,'Fontname','Garamond')
        saveas(figC,namefigC);
        
        
        
    end
    
    
end

cd(inDIR_NoN);
dlmwrite('percolation.txt',percolation,'delimiter','\t');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               BLOCKS MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% C_blocks = zeros(M,M);
%
%
% imax = size_mod(1,2) % size of first cluster
% jmin = 1;
% jmax = size_mod(1,2) + size_mod(2,2)
%
% for module = 1:(max(voxel(:,6))-1) % from 1 to the max numb of cluster
%
%     for i = imax % within a given block - rows
%
%         for j  = jmin:jmax      % within a given block - columns
%             C_blocks(i,j) = 1;
%             C_blocks(j,i) = C_blocks(i,j);
%         end
%     end
%     jmin = imax + 1;    % shift to next block
%     imax = imax + size_mod(module + 1,2);  % shift to next block
%     size_mod(module+2,2);
%     jmax = jmax + size_mod(module+2,2);
% end





