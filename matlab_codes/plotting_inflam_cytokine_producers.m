A = 'output0000000';
A2 = 'output000000';
A3 = 'output00000';
B = '.xml';

tcount = 24;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing1 = MCDS.discrete_cells.metadata.type(locs_proinflam);

tcount = 24*3;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing3 = MCDS.discrete_cells.metadata.type(locs_proinflam);


tcount = 24*5;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing5 = MCDS.discrete_cells.metadata.type(locs_proinflam);

tcount = 24*7;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing7 = MCDS.discrete_cells.metadata.type(locs_proinflam);

tcount = 24*9;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing9 = MCDS.discrete_cells.metadata.type(locs_proinflam);

tcount = 24*11;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing11 = MCDS.discrete_cells.metadata.type(locs_proinflam);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing9 = MCDS.discrete_cells.metadata.type(locs_proinflam);

tcount = 24*13;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing13 = MCDS.discrete_cells.metadata.type(locs_proinflam);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing9 = MCDS.discrete_cells.metadata.type(locs_proinflam);

tcount = 24*15;
if tcount<11
        K = [A num2str(tcount-1,'%d') B];
    elseif tcount<101
        K = [A2 num2str(tcount-1,'%d') B];
    else
        K = [A3 num2str(tcount-1,'%d') B];
    end
MCDS = read_MultiCellDS_xml(K);
    
% find cells producing pro-inflammatory cytokine
locs_proinflam = find(MCDS.discrete_cells.phenotype.transport_processes.variable(3).export_rate>0);
types_producing15 = MCDS.discrete_cells.metadata.type(locs_proinflam);

%plot bar chart of the type of cells producing pro-inflammatory cytokine
figure
subplot(3,3,1)
histogram(types_producing1)
title('DAY 1')
subplot(3,3,2)
histogram(types_producing3)
title('DAY 3')
subplot(3,3,3)
histogram(types_producing5)
title('DAY 5')
subplot(3,3,4)
histogram(types_producing7)
title('DAY 7')
subplot(3,3,5)
histogram(types_producing9)
title('DAY 9')
subplot(3,3,6)
histogram(types_producing11)
title('DAY 11')
subplot(3,3,7)
histogram(types_producing13)
title('DAY 13')
subplot(3,3,8)
histogram(types_producing15)
title('DAY 15')

