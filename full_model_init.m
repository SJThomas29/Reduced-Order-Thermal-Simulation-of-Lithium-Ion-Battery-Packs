%% Full model configuration
% See parameter justification for sources and explanation for parameters

% Initialise with user inputs if not already done
disp('Checking for previous init...');
clear; % Remove after sens
if ~exist('init', 'var')
    clear;
    s = false; % If this is true, inputs below will be used
    inputLayout = ["4" "3" "9" "3" "5" "3" "1" "y" "z" "y"];% - default
    inputTherm = ["25" "none" "+y" "1" "0" "1" "1" "0" "0"];% - default
    inputBoundSim = ["nan" "nan" "nan" "nan" "nan" "nan" "7" "1" "0.1" "10"]; % - default
    disp('Model not initialised, collecting variables from user...')
    [conf, confOk] = GetInputs(s,inputLayout,inputTherm,inputBoundSim);
    if ~confOk; clear; return; end
else
    disp('Model already initialised, proceeding to simulation...');
    return;
end

%% Global configuration
% All parameters and constants for the model are defined here
% Any parameters used in the model direction will be stored in a structure,
% this helps to de-clutter the workspace

% Load parameters for the cell model, materials dictionary and
% radiator/refrigerant tune
dataPath = [pwd '\Data\'];
load([dataPath 'Look_up_tables_32Ah_V2.mat']);
matDict = readtable([dataPath 'MaterialsDictionary.xlsx']);
load([dataPath 'rad_20C_tune.mat']);
load([dataPath 'ref_20C_tune.mat']);

% Vehicle and battery pack model based on a combination of Audi etron 50 and 55
% Vehicle model
vehM = 2565; % Vehicle mass [kg] - e-tron 50
vehCd = 0.28; % Drag coefficient [#] - e-tron 50 
vehAf = 3.15; % Frontal area [m^2] - e-tron 50 
% Cell params
cellThick = 16.5e-3; % Cell thickness [m] - etron 55
cellWidth = 100e-3; % Cell width [m] - etron 55
cellLength = 330e-3; % Cell length [m] - etron 55
cellEnergyEtron = 220; % Cell energy [Wh] - etron 55
elecCellCap = 60; % Cell capacity [Ah] - etron 55
% Pack dimensions 
cellWidthWallThickness = 4e-3; % Wall thickness at the width of the module [m] - etron 55 (est.)
cellLengthWallThickness = 2*cellWidthWallThickness; % Wall thickness at the length of the module [m] - etron 55 (est.)
crashStructureThickness = 6.5e-3; % Crash structure thickness inbetween modules [m] - etron 55 (est.) - higher to account for extra casing
packEnergy = 71000; % Total pack energy [Wh] - e-tron 50
coolThickness = 10e-3; % Thickness of coolant channel [m]

% Scale Team Giant cell model for the e-tron cells
% Reduce resistance linearly with capacity and do the opposite to capacitance
cellCapGiant = 32; % Cell capacity [Ah] - Team Giant
cellEnergyGiant = 118.65; % Cell energy [Wh] - Team Giant
elecR0 = R1_matrix*cellCapGiant/elecCellCap; % Scaled R0 [Ohm]
elecR1 = R2_matrix*cellCapGiant/elecCellCap; % Scaled R1 [Ohm]
elecC1 = C1_matrix*elecCellCap/cellCapGiant; % Scaled C1 [Ohm]
elecSOCBreak = SOC_breakpoints; % State of charge breakpoints [#]
elecTempBreak = Temp_breakpoints; % Temperature breakpoints [degC]
elecOCV = Em_matrix; % Open circuit voltage [V]

% Calculate power multiplier based on model and etron capacity
cellsS = conf.modules.x*conf.modules.y*conf.cells.s; % Series cells in pack [#]
cellsP = conf.cells.p; % Parallel cells in pack [#]
cellsTot = cellsS*cellsP; % Total cells in pack [#]
modelEnergy = cellsS*cellsP*cellEnergyEtron;%cellEnergyEtron; % Model energy [Wh]
vehPowerMult = modelEnergy/packEnergy; % Vehicle model power multiplier [#]

%% Electrical configuration 
% Inputs to the electrical modeul built in this script
% loadDenomInds - indices for the loading calculator denominator
% loadTerm12Inds - indices for terms 1 and 2 where p is 2
% loadTerm1Inds - indices for term 1 where p is more than 2
% loadTerm2Inds - indices for term 1 where p is more than 2
% loadTerm3Inds - indices for term 1 where p is more than 2

% Assign initial state of charge with some variation
variation = 0; % Variation amount for SOC and temp [#]
socAve = 0.95; % Average SOC value [#]
elecInitSOC = (1 - variation)*socAve + 2*variation*socAve*rand(cellsS,cellsP);  % Initial state of charge matrix (SxP) [#]

% Initalise matrices of indices for current distribution calculations
% Done using indices as it speeds up operations at each step
if cellsP == 1
    % No matrices required, same current in all cells
    % Define empty parameters to avoid matlab function errors
    elecLoadDenomInds = 0;
    elecLoadTerm12Inds = 0;
    elecLoadTerm1Inds = 0;
    elecLoadTerm2Inds = 0;
    elecLoadTerm3Inds = 0;
elseif cellsP == 2
    % Defining indices for the loading determination
    for i = 1:cellsS
        % Denominator
        % Indices of the two resistances in each parallel branch
        elecLoadDenomInds(:,:,i,1) = [i i + cellsS];
        elecLoadDenomInds(:,:,i,2) = [i i + cellsS];
        % Terms 1 and 2 equal in this configuration
        % Index of the other cell in the parallel branch
        elecLoadTerm12Inds(:,:,i,1) = i + cellsS;
        elecLoadTerm12Inds(:,:,i,2) = i;
    end
    % Define empty parameters to avoid matlab function errors
    elecLoadTerm1Inds = 0;
    elecLoadTerm2Inds = 0;
    elecLoadTerm3Inds = 0;
else
    % Defining indices for the loading determination
    for i = 1:cellsS
        branchStartInd = i; % Start of parallel branch
        branchEndInd = i + (cellsP - 1)*cellsS; % End of parallel branch
        branchInds = branchStartInd:cellsS:branchEndInd; % All cells in branch
        for j = 1:cellsP
            % Denominator
            % Indices for the sum of resistance combinations in the parallel
            % branch of length parallel count - 1 
            elecLoadDenomInds(:,:,i,j) = nchoosek(branchInds,cellsP - 1);
            % Term 1
            % Indices for the sum of resistance combinations in the parallel
            % branch not including the focused cell of length parallel count - 2 
            loadTerm1IndsFull = nchoosek(branchInds,cellsP - 2); % All combinations
            elecLoadTerm1Inds(:,:,i,j) = loadTerm1IndsFull(~any(loadTerm1IndsFull == branchInds(j),2),:); % Removing combinations including focused
            % Term 2
            % Indices above alongside the index for the missing voltage not
            % including the focused cell. Only used in P > 2.
            % e.g. for 4 parallel focusing on cell 1, [R2 R3] becomes [V4 R2 R3]
            if cellsP > 2 
                term1Temp = elecLoadTerm1Inds(:,:,i,j);
                term2Temp = [zeros(height(term1Temp),1) term1Temp];
                vNeeded = branchInds; vNeeded(j) = []; % All missing voltages
                for k = 1:length(vNeeded)
                    missing = logical(floor(mean(~ismember(term1Temp,vNeeded(k)),2))); % Finds where the missing voltage should go
                    term2Temp(missing,1) = vNeeded(k); % Assignes voltage to correct place
                end
                elecLoadTerm2Inds(:,:,i,j) = term2Temp;
            end
            % Term 3
            % Indices for all resistances in the parallel branch that do not
            % include the focused cell
            inds = branchInds; inds(j) = []; % Inds on branch not including focused
            elecLoadTerm3Inds(:,:,i,j) = inds;                 
        end
    end
    % Define empty parameters to avoid matlab function errors
    elecLoadTerm12Inds = 0;
end

%% Thermal configuration
% Inputs to thermal model built in this script
% Where n1,n2 and n3 are the number of segments in the x,y,z direction respectively
% x - x boundary positions [n1 + 1, 1]
% y - y boundary positions [n2 + 1, 1]
% z - z boundary positions [n3 + 1, 1]
% C - heat capacity of each volume [n1 x n2 x n3]
% Kx - thermal conductivity at boundary in the x direction [n1 + 1, n2, n3]
% Ky - thermal conductivity at boundary in the y direction [n1, n2 + 1, n3]
% Kz - thermal conductivity at boundary in the z direction [n1, n2, n3 + 1]
% mdotDistX - flow distribution in the x direction [n1, n2, n3]
% mdotDistY - flow distribution in the y direction [n1, n2, n3]
% Inputs to be defined at realtime based on control - mdotZ could be added in future work using the same methods
% mdotX - mass flow rate for each volume in the x direction [n1, n2, n3] - function of time (uses mdotDistX)
% mdotY - mass flow rate for each volume in the y direction [n1, n2, n3] - function of time (uses mdotDistY)
% Q - heat generation of each volume [n1, n2, n3] - function of time

% Define ambient temp
tempAmb = conf.bound.amb; % Ambient temperature [degC]

% Define x,y,z vectors using configuration
% First define a single cell
cell = struct; % Base structure
cell.l = [0; (cellLength*(1:conf.regions.l)/conf.regions.l)']; % Length
cell.w = [0; (cellWidth*(1:conf.regions.w)/conf.regions.w)']; % Width
cell.t = [0; (cellThick*(1:conf.regions.t)/conf.regions.t)']; % Thickness
% Uncomment to plot cell
% PlotGrid(cell.l, cell.w, cell.t, 'Cell', '3d'); % view can be, 3d, xy, xz, yz
% axis equal;

% Define the cell stack and rotate based on configuration
cellTemp = cell; % Temporary variable for stacking
for i = 2:conf.cells.tot
    prev = cellTemp.t(end); % End of last cell
    cellTemp.t = [cellTemp.t; prev + cell.t(2:end)]; % Add next cell
end
cellStack.(conf.orient.l) = cellTemp.l; % _CellStack is equal to lCell where _ is the length dimension
cellStack.(conf.orient.w) = cellTemp.w; % _CellStack is equal to wCell where _ is the width dimension
cellStack.(conf.orient.t) = cellTemp.t; % _CellStack is equal to tCell where _ is the thickness dimension
% Uncomment to plot cell stack
% plot_grid(cellStack.x, cellStack.y, cellStack.z, 'Cell stack', '3d'); % view can be, 3d, xy, xz, yz
% axis equal;

% Define full pack using cell stacks, walls and cooling - to model structure
% Combines walls and structural support of neighboroughing modules to
% reduce the overall number of volumes
thermX = 0; % Initialise empty variable
for i = 1:conf.modules.x % Stack modules in the x direction
    % Check which walls are required
    if conf.orient.l == 'x' % End of module lengthways - thicker walls
       wallThickness = cellLengthWallThickness;
    else
       wallThickness = cellWidthWallThickness;
    end
    % Add module - wall/crash - cells - wall/crash
    if strcmp(conf.crash, 'y') % Crash structure in use 
        thermX = [thermX; thermX(end) + wallThickness]; % first wall
        thermX = [thermX; thermX(end) + cellStack.x(2:end)]; % cell stack
        thermX = [thermX; thermX(end) + wallThickness]; % last wall
        if i < conf.modules.x % Only add crash if not the last
            thermX = [thermX; thermX(end) + crashStructureThickness]; % crash structure
        end
    else
        thermX = [thermX; thermX(end) + wallThickness]; % first wall
        thermX = [thermX; thermX(end) + cellStack.x(2:end)]; % cell stack
        thermX = [thermX; thermX(end) + wallThickness]; % last wall
    end
end
thermY = 0; % Initialise empty variable
for i = 1:conf.modules.y % Stack modules in the y direction
    % Check which walls are required
    if conf.orient.l == 'y' % End of module lengthways - thicker walls
       wallThickness = cellLengthWallThickness;
    else
       wallThickness = cellWidthWallThickness;
    end
    % Add module - wall/crash - cells - wall/crash
    if strcmp(conf.crash, 'y') % Crash structure in use 
        thermY = [thermY; thermY(end) + wallThickness]; % first wall
        thermY = [thermY; thermY(end) + cellStack.y(2:end)]; % cell stack
        thermY = [thermY; thermY(end) + wallThickness]; % last wall
        if i < conf.modules.y % Only add crash if not the last
            thermY = [thermY; thermY(end) + crashStructureThickness]; % crash structure
        end
    else
        thermY = [thermY; thermY(end) + wallThickness]; % first wall
        thermY = [thermY; thermY(end) + cellStack.y(2:end)]; % cell stack
        thermY = [thermY; thermY(end) + wallThickness]; % last wall
    end
end
thermZ = 0; % Simple as there is only one module in z direction
% Check which walls are required
if conf.orient.l == 'z' % End of module lengthways - thicker walls and gap
   wallThickness = cellLengthWallThickness;
else
   wallThickness = cellWidthWallThickness;
end
% Add main body of module - wall - cells - wall
thermZ = [thermZ; thermZ(end) + wallThickness]; % first wall
thermZ = [thermZ; thermZ(end) + cellStack.z(2:end)]; % cell stack
thermZ = [thermZ; thermZ(end) + wallThickness]; % last wall
% Add coolant if needed
topCool = abs(conf.cooling.t.x) + abs(conf.cooling.t.y) == 1;
if topCool % Top cooling in use
    thermZ = [thermZ; thermZ(end) + coolThickness]; % Add end coolant
end
bottomCool =  abs(conf.cooling.b.x) + abs(conf.cooling.b.y) == 1;
if bottomCool % Top cooling in use
    thermZ = [thermZ(1); thermZ(1) + coolThickness; thermZ(2:end) + coolThickness]; % Add start coolant
end

% Define pack dimensions - number of volumes
thermDims(1) = length(thermX) - 1;
thermDims(2) = length(thermY) - 1;
thermDims(3) = length(thermZ) - 1;
thermDimsStruct.x = thermDims(1); % Temp struct to be deleted after init
thermDimsStruct.y = thermDims(2);
thermDimsStruct.z = thermDims(3);
% Uncomment to plot pack
% plot_grid(thermX, thermY, thermZ, 'Pack', '3d'); % view can be, 3d, xy, xz, yz
% axis equal;

% Define Kx, Ky, Kz and Cp matrices using configuration
% First step, define material matrix, this will be [n1 x n2 x n3] containing the
% ID of the volume material. Also define another matrix of the same
% dimensions containing a 1 for in plane conductivity and 2 for cross plane
% Define direction of heat transfer for cells
dir.(conf.orient.l) = 1; % length direction
dir.(conf.orient.w) = 1; % width direction
dir.(conf.orient.t) = 2; % thickness direction
% Define volumes that contain aluminium
if conf.orient.l == 'z' % Different approach for z modules, no crash
    al.(conf.orient.l) = (1:conf.modules.(conf.orient.l) + 1) + ((1:conf.modules.(conf.orient.l) + 1) - 1)*conf.regions.l;
else
    if strcmp(conf.crash, 'y') % crash structure in use
        % mid vals - crash structure layer, low and high are modules 
        vals = (2:2:2*(conf.modules.(conf.orient.l) - 1)) + ((2:conf.modules.(conf.orient.l)) - 1)*conf.regions.l + ((2:conf.modules.(conf.orient.l)) - 1);
        vals = [vals, vals - 1, vals + 1]; % add high and low values
    else
        % high vals - second layer at each boundary, lowVals is the first
        vals = (2:conf.modules.(conf.orient.l)) + ((2:conf.modules.(conf.orient.l)) - 1)*conf.regions.l + ((2:conf.modules.(conf.orient.l)) - 1);
        vals = [vals, vals - 1]; % add low values
    end
    vals = [vals 1 thermDimsStruct.(conf.orient.l)]; % add start and end values
    al.(conf.orient.l) = vals;
    al.(conf.orient.l) = sort(al.(conf.orient.l));
end
if conf.orient.w == 'z' % Different approach for z modules, no crash
    al.(conf.orient.w) = (1:conf.modules.(conf.orient.w) + 1) + ((1:conf.modules.(conf.orient.w) + 1) - 1)*conf.regions.w;
else
    if strcmp(conf.crash, 'y') % crash structure in use
        % mid vals - crash structure layer, low and high are modules 
        vals = (2:2:2*(conf.modules.(conf.orient.w) - 1)) + ((2:conf.modules.(conf.orient.w)) - 1)*conf.regions.w + ((2:conf.modules.(conf.orient.w)) - 1);
        vals = [vals, vals - 1, vals + 1]; % add high and low values
    else
        % high vals - second layer at each boundary, lowVals is the first
        vals = (2:conf.modules.(conf.orient.w)) + ((2:conf.modules.(conf.orient.w)) - 1)*conf.regions.w + ((2:conf.modules.(conf.orient.w)) - 1);
        vals = [vals, vals - 1]; % add low values
    end
    vals = [vals 1 thermDimsStruct.(conf.orient.w)]; % add start and end values
    al.(conf.orient.w) = vals;
    al.(conf.orient.w) = sort(al.(conf.orient.w));
end
if conf.orient.t == 'z' % Different approach for z modules, no crash
    al.(conf.orient.t) = (1:conf.modules.(conf.orient.t) + 1) + ((1:conf.modules.(conf.orient.t) + 1) - 1)*conf.regions.t*conf.cells.tot;
else
    if strcmp(conf.crash, 'y') % crash structure in use
        % mid vals - crash structure layer, low and high are modules 
        vals = (2:2:2*(conf.modules.(conf.orient.t) - 1)) + ((2:conf.modules.(conf.orient.t)) - 1)*conf.regions.t*conf.cells.tot + ((2:conf.modules.(conf.orient.t)) - 1);
        vals = [vals, vals - 1, vals + 1]; % add high and low values
    else
        % high vals - second layer at each boundary, lowVals is the first
        vals = (2:conf.modules.(conf.orient.t)) + ((2:conf.modules.(conf.orient.t)) - 1)*conf.regions.t*conf.cells.tot + ((2:conf.modules.(conf.orient.t)) - 1);
        vals = [vals, vals - 1]; % add low values
    end
    vals = [vals 1 thermDimsStruct.(conf.orient.t)]; % add start and end values
    al.(conf.orient.t) = vals;
    al.(conf.orient.t) = sort(al.(conf.orient.t));
end
% Add values to z aluminium if there is bottom cooling
if bottomCool; al.z = al.z + 1; end
% Define material indices to model structure
thermMatID = ones(thermDims); % all to default (cell)
thermMatID(al.x,:,:) = 2; % x values to aluminium
thermMatID(:,al.y,:) = 2; % y values to aluminium
thermMatID(:,:,al.z) = 2; % z values to aluminium
if bottomCool; thermMatID(:,:,1) = 3; end % Bottom coolant
if topCool; thermMatID(:,:,end) = 3; end % Top coolant
% Second step, define Cp, Kx, Ky and Kz matrices using the material dictionary
% Apply volume based specific heat capacity based on material indices
thermCp = matDict.Mass_Cp(thermMatID);
% Apply material thermal conductivites in each volume in the x,y,z directions,
% uses the cell direction for all as the aluminium is homogenous
KxMat = reshape(matDict(thermMatID,4 + dir.x).Variables,size(thermMatID));
KyMat = reshape(matDict(thermMatID,4 + dir.y).Variables,size(thermMatID));
KzMat = reshape(matDict(thermMatID,4 + dir.z).Variables,size(thermMatID));
% Calculate thermal conductivities at each boundary based on the material
% properties and volume dimensions
% Converting from thermal conductivity in the form W/mK to W/K for the
% simulation using the following equation:
% K[W/K] = K[W/mK]*A[m^2]/L[m] - L (length), k (material conductivity), A (area)
% Get length of volumes in x,y,z directions and map to pack dimensions
Lx = repmat(thermX(2:end) - thermX(1:end - 1),[1 thermDims(2) thermDims(3)]);
Ly = repmat((thermY(2:end) - thermY(1:end - 1))',[thermDims(1) 1 thermDims(3)]);
Lz = []; Lz(1,1,:) = thermZ(2:end) - thermZ(1:end - 1);
Lz = repmat(Lz,[thermDims(1) thermDims(2) 1]);
% Get cross sectional area of volumes in x,y,z directions
% Need to duplicate first layer in direction to match dimensions of boundaries
Ax = Ly.*Lz;
Ay = Lx.*Lz;
Az = Lx.*Ly;
% Calculate thermal conductivities [W/K] of volumes in x,y,z directions
Kx = (KxMat.*Ax)./Lx;
Ky = (KyMat.*Ay)./Ly;
Kz = (KzMat.*Az)./Lz;

% Apply boundary conductivities based on the config
% Define default values using averages from adjacent volumes
tempK.x = convn(Kx,[0.5; 0.5]); 
tempK.y = convn(Ky,[0.5 0.5]);
tempK.z = convn(Kz,0.5*ones(1,1,2));
% Define values that arent aluminium for boundary overwrites
cellVals.x = 1:thermDims(1); cellVals.x(al.x) = [];
cellVals.y = 1:thermDims(2); cellVals.y(al.y) = [];
cellVals.z = 1:thermDims(3); cellVals.z(al.z) = [];
if topCool; cellVals.z(end) = []; end % Dont adjust top if cooling
if bottomCool; cellVals.z(1) = []; end % Dont adjust botom if cooling
% Redefine areas with same dimensions as boundaries
AxBound = Ax; AxBound(end+1,:,:) = AxBound(end,:,:);
AyBound = Ay; AyBound(:,end+1,:) = AyBound(:,end,:);
AzBound = Az; AzBound(:,:,end+1) = AzBound(:,:,end);
% Calculate all boundary conductivities as if they used the pouch
pouchCond = 0.269; % Pouch thermal cond [W/mK]
pouchThick = 0.254e-3; % Pouch thickness [m]
pouchK.x = (pouchCond*AxBound)/pouchThick;
pouchK.y = (pouchCond*AyBound)/pouchThick;
pouchK.z = (pouchCond*AzBound)/pouchThick;
% Calculate all boundary conductivities as if they were pure radiation
radUnitCond = 5.4268; % Linear radiation thermal conductivity with a unit area
radK.x = radUnitCond*AxBound;
radK.y = radUnitCond*AyBound;
radK.z = radUnitCond*AzBound;
% Calculate all boundary conductivities as if they used thermal gel 
thermGelThickness = 1.5e-3; % Thickness of thermal gel connection to coolant [m]
thermGelCond = 3.5; % Thermal gel conductivity [W/mk]
gelK.x = (thermGelCond*AxBound)/thermGelThickness; 
gelK.y = (thermGelCond*AyBound)/thermGelThickness; 
gelK.z = (thermGelCond*AzBound)/thermGelThickness; 
% % Cell - cell % %
% Minimum - radiation, Maximum - pouch conductivity
% Get volumes with cell boundaries - only in cell thickness direction
cellBound.(conf.orient.l) = []; % No boundaries in length
cellBound.(conf.orient.w) = []; % No boundaries in width
cellBound.(conf.orient.t) = []; % Initialise thickness boundary
for i = 1:length(al.(conf.orient.t)) - 1
    vals = (al.(conf.orient.t)(i):conf.regions.t:al.(conf.orient.t)(i + 1) - 1);
    vals = vals(2:end - 1) + 1; % Trim start and end values
    cellBound.(conf.orient.t) = [cellBound.(conf.orient.t) vals];
end
% Apply cell-cell conductivities, half the cond as thickness is doubled
tempK.x(cellBound.x,cellVals.y,cellVals.z) = radK.x(cellBound.x,cellVals.y,cellVals.z) + conf.cond.cellCell*(pouchK.x(cellBound.x,cellVals.y,cellVals.z)/2 - radK.x(cellBound.x,cellVals.y,cellVals.z));
tempK.y(cellVals.x,cellBound.y,cellVals.z) = radK.y(cellVals.x,cellBound.y,cellVals.z) + conf.cond.cellCell*(pouchK.y(cellVals.x,cellBound.y,cellVals.z)/2 - radK.y(cellVals.x,cellBound.y,cellVals.z));
tempK.z(cellVals.x,cellVals.y,cellBound.z) = radK.z(cellVals.x,cellVals.y,cellBound.z) + conf.cond.cellCell*(pouchK.z(cellVals.x,cellVals.y,cellBound.z)/2 - radK.z(cellVals.x,cellVals.y,cellBound.z));
% Define module boundaries for use in boundary conductivities below
if conf.orient.l == 'z' % Different approach for z modules, no crash
    modBound.(conf.orient.l) = unique([al.(conf.orient.l), al.(conf.orient.l) + 1]);
    modBound.(conf.orient.l)(1) = []; modBound.(conf.orient.l)(end) = [];
else
    modBound.(conf.orient.l) = unique([al.(conf.orient.l), al.(conf.orient.l) + 1]);
    modBound.(conf.orient.l)(1) = []; modBound.(conf.orient.l)(end) = [];
    if strcmp(conf.crash, 'y')
        vals = (2:2:2*(conf.modules.(conf.orient.l) - 1)) + ((2:conf.modules.(conf.orient.l)) - 1)*conf.regions.l + ((2:conf.modules.(conf.orient.l)) - 1);
        vals = [vals, vals + 1]; % add high and low values
    else
        vals = (2:conf.modules.(conf.orient.l)) + ((2:conf.modules.(conf.orient.l)) - 1)*conf.regions.l + ((2:conf.modules.(conf.orient.l)) - 1);
    end
    modBound.(conf.orient.l)(ismember(modBound.(conf.orient.l),vals)) = []; % Remove middle values, al to al
end
if conf.orient.w == 'z' % Different approach for z modules, no crash
    modBound.(conf.orient.w) = unique([al.(conf.orient.w), al.(conf.orient.w) + 1]);
    modBound.(conf.orient.w)(1) = []; modBound.(conf.orient.w)(end) = [];
else
    modBound.(conf.orient.w) = unique([al.(conf.orient.w), al.(conf.orient.w) + 1]);
    modBound.(conf.orient.w)(1) = []; modBound.(conf.orient.w)(end) = [];
    if strcmp(conf.crash, 'y')
        vals = (2:2:2*(conf.modules.(conf.orient.w) - 1)) + ((2:conf.modules.(conf.orient.w)) - 1)*conf.regions.w + ((2:conf.modules.(conf.orient.w)) - 1);
        vals = [vals, vals + 1]; % add high and low values
    else
        vals = (2:conf.modules.(conf.orient.w)) + ((2:conf.modules.(conf.orient.w)) - 1)*conf.regions.w + ((2:conf.modules.(conf.orient.w)) - 1);
    end
    modBound.(conf.orient.w)(ismember(modBound.(conf.orient.w),vals)) = []; % Remove middle values, al to al
end
if conf.orient.t == 'z' % Different approach for z modules, no crash
    modBound.(conf.orient.t) = unique([al.(conf.orient.t), al.(conf.orient.t) + 1]);
    modBound.(conf.orient.t)(1) = []; modBound.(conf.orient.t)(end) = [];
else
    modBound.(conf.orient.t) = unique([al.(conf.orient.t), al.(conf.orient.t) + 1]);
    modBound.(conf.orient.t)(1) = []; modBound.(conf.orient.t)(end) = [];
    if strcmp(conf.crash, 'y')
        vals = (2:2:2*(conf.modules.(conf.orient.t) - 1)) + ((2:conf.modules.(conf.orient.t)) - 1)*conf.regions.t*conf.cells.tot + ((2:conf.modules.(conf.orient.t)) - 1);
        vals = [vals, vals + 1]; % add high and low values
    else
        vals = (2:conf.modules.(conf.orient.t)) + ((2:conf.modules.(conf.orient.t)) - 1)*conf.regions.t*conf.cells.tot + ((2:conf.modules.(conf.orient.t)) - 1);
    end
    modBound.(conf.orient.t)(ismember(modBound.(conf.orient.t),vals)) = []; % Remove middle values, al to al
end
% % Cell - module % %
% Minimum - radiation, Maximum - pouch conductivity
% Length direction
cellValsTemp.(conf.orient.l) = modBound.(conf.orient.l);
cellValsTemp.(conf.orient.w) = cellVals.(conf.orient.w);
cellValsTemp.(conf.orient.t) = cellVals.(conf.orient.t);
tempK.(conf.orient.l)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z) =  radK.(conf.orient.l)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z) + conf.cond.cellLenMod*(pouchK.(conf.orient.l)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z)/2 - radK.(conf.orient.l)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z));
% Width direction
cellValsTemp.(conf.orient.l) = cellVals.(conf.orient.l);
cellValsTemp.(conf.orient.w) = modBound.(conf.orient.w);
cellValsTemp.(conf.orient.t) = cellVals.(conf.orient.t);
tempK.(conf.orient.w)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z) =  radK.(conf.orient.w)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z) + conf.cond.cellWidMod*(pouchK.(conf.orient.w)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z)/2 - radK.(conf.orient.w)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z));
% Thickness direction
cellValsTemp.(conf.orient.l) = cellVals.(conf.orient.l);
cellValsTemp.(conf.orient.w) = cellVals.(conf.orient.w);
cellValsTemp.(conf.orient.t) = modBound.(conf.orient.t);
tempK.(conf.orient.t)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z) =  radK.(conf.orient.t)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z) + conf.cond.cellThickMod*(pouchK.(conf.orient.t)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z)/2 - radK.(conf.orient.t)(cellValsTemp.x,cellValsTemp.y,cellValsTemp.z));
% % Module - crash structure/module % %
% Minimum - radiation, Maximum - thermal gel
% Crash structure/module boundaries span full pack, apart from coolant
valsX = 1:thermDims(1); valsY = 1:thermDims(2); valsZ = 1:thermDims(3);
boundsX = modBound.x(2:2:end-1); % Remove start and end
boundsY = modBound.y(2:2:end-1); % Remove start and end
if topCool; valsZ(end) = []; end % remove top if cooled
if bottomCool; valsZ(1) = []; end % remove bottom if cooled
% Boundaries vary based on the use of crash
if strcmp(conf.crash, 'y')
    boundsX = [boundsX + 1, boundsX + 2]; boundsX = sort(boundsX);
    boundsY = [boundsY + 1, boundsY + 2]; boundsY = sort(boundsY);
    valsX(ismember(valsX,boundsX(1:2:end))) = [];
    valsY(ismember(valsY,boundsY(1:2:end))) = [];
else
    boundsX = boundsX + 1;
    boundsY = boundsY + 1;
end
% Apply boundary conductivities
tempK.x(boundsX,valsY,valsZ) = radK.x(boundsX,valsY,valsZ) + conf.cond.modCrash*(gelK.x(boundsX,valsY,valsZ) - radK.x(boundsX,valsY,valsZ));
tempK.y(valsX,boundsY,valsZ) = radK.y(valsX,boundsY,valsZ) + conf.cond.modCrash*(gelK.y(valsX,boundsY,valsZ) - radK.y(valsX,boundsY,valsZ));
% % Module - coolant % %
% Locked as thermal gel
if topCool; tempK.z(1:thermDims(1),1:thermDims(2),end - 1) = gelK.z(1:thermDims(1),1:thermDims(2),end - 1); end % Top overwrite
if bottomCool; tempK.z(1:thermDims(1),1:thermDims(2),2) = gelK.z(1:thermDims(1),1:thermDims(2),2); end % Bottom overwrite

% Define terms for model from temp structure
thermKx = tempK.x;
thermKy = tempK.y;
thermKz = tempK.z;

% Apply insulation boundary conditions
if isnan(conf.bound.xa) % x start insulation
    thermKx(1,:,:) = 0;
    thermBoundXA = 0;
else
    thermBoundXA = conf.bound.amb;
end
if isnan(conf.bound.xb) % x end insulation
    thermKx(end,:,:) = 0;
    thermBoundXB = 0;
else
    thermBoundXB = conf.bound.amb;
end
if isnan(conf.bound.ya) % y start insulation
    thermKy(:,1,:) = 0;
    thermBoundYA = 0;
else
    thermBoundYA = conf.bound.amb;
end
if isnan(conf.bound.yb) % y end insulation
    thermKy(:,end,:) = 0;
    thermBoundYB = 0;
else
    thermBoundYB = conf.bound.amb;
end
if isnan(conf.bound.za) % z start insulation
    thermKz(:,:,1) = 0;
    thermBoundZA = 0;
else
    thermBoundZA = conf.bound.amb;
end
if isnan(conf.bound.zb) % z end insulation
    thermKz(:,:,end) = 0;
    thermBoundZB = 0;
else
    thermBoundZB = conf.bound.amb;
end

% Use lengths and areas from above to create mass matrix
% Calculate volumes
V = Lx.*Ly.*Lz;
% Get densities from material dictionary
rho = matDict.Density(thermMatID);
% Calculate masses from volume and density
thermM = rho.*V;

% Define mdotDistX and mdotDistY based on the areas and flow directions
% Assuming there is one pump for all flow, so flow is halved if top and
% bottom cooling is being used
if topCool && bottomCool
    flowMult = 0.5; 
else
    flowMult = 1;
end
mdotDistX = zeros(thermDims); mdotDistY = mdotDistX; % Init empty matrices
if topCool % top cooling
    if conf.cooling.t.x ~= 0 % x direction flow
        areas = Ax(1,:,end);
        areasDist = areas/sum(areas); % distribution of areas
        mdotDistX(:,:,end) = conf.cooling.t.x*repmat(areasDist,thermDims(1),1); % Accounts for direction of flow
    else % y direction flow
        areas = Ay(:,1,end); % areas of each 'tube'
        areasDist = areas/sum(areas); % distribution of areas
        mdotDistY(:,:,end) = conf.cooling.t.y*repmat(areasDist,1,thermDims(2)); % Accounts for direction of flow
    end
end
if bottomCool % top cooling
    if conf.cooling.b.x ~= 0 % x direction flow
        areas = Ax(1,:,1);
        areasDist = areas/sum(areas); % distribution of areas
        mdotDistX(:,:,1) = conf.cooling.b.x*repmat(areasDist,thermDims(1),1); % Accounts for direction of flow
    else % y direction flow
        areas = Ay(:,1,1); % areas of each 'tube'
        areasDist = areas/sum(areas); % distribution of areas
        mdotDistY(:,:,1) = conf.cooling.b.y*repmat(areasDist,1,thermDims(2)); % Accounts for direction of flow
    end
end
% Apply flow rate multiplier to split coolant if using top and bottom
mdotDistX = mdotDistX*flowMult; 
mdotDistY = mdotDistY*flowMult; 

% Mapping the cell positions to the heat generation for the map and vice versa
% Order of the modules does not matter as they all recieve the same current
% Parallel cell are placed next to eachother as done in the audi modules
% Get relevant indices for mapping from electrical to thermal
thermHeatInds = zeros(thermDims); % Empty matrix
m = 1; % Module counter
l = conf.orient.l; w =conf.orient.w; t = conf.orient.t; % To save space
for lInd = 1:conf.modules.(l) % Loop through modules in cell length direction
    for wInd = 1:conf.modules.(w) % Loop through modules in cell width direction 
        for tInd = 1:conf.modules.(t) % Loop through modules in cell thickness direction  
            for s = 1:conf.cells.s % Loop through series cells in the module
                for p = 1:conf.cells.p % Loop through parallel cells in the module
                    % Get cell id
                    cellID = s + conf.cells.s*(m - 1) + (p - 1)*cellsS;
                    % Get cell positon in module based on series and parallel                
                    cellPos = (s - 1)*conf.cells.p + p;
                    % Get vector of volumes related to each cell
                    % Length and width don't account for cell position
                    % First term produces a vector for the number of regions
                    % Second adds values to skip previous modules
                    % Third accounts for layers of aluminium between modules
                    if strcmp(conf.crash, 'y')
                        mult = 3;
                    else
                        mult = 2;
                    end
                    if conf.orient.l == 'z' % Different approach for z modules, no crash
                        cellTemp.(l) = (1:conf.regions.l) + (lInd - 1)*conf.regions.l + lInd; 
                    else
                        cellTemp.(l) = (1:conf.regions.l) + (lInd - 1)*conf.regions.l + 1 + mult*(lInd - 1); 
                    end
                    if conf.orient.w == 'z' % Different approach for z modules, no crash
                        cellTemp.(w) = (1:conf.regions.w) + (wInd - 1)*conf.regions.w + wInd; 
                    else
                        cellTemp.(w) = (1:conf.regions.w) + (wInd - 1)*conf.regions.w + 1 + mult*(wInd - 1); 
                    end
                    % First and third terms are same
                    % Second term also accounts for multiple cells in the direction per module
                    % Additional term accouts for position of the cell in the module
                    if conf.orient.t == 'z' % Different approach for z modules, no crash
                        cellTemp.(t) = (1:conf.regions.t) + (tInd - 1)*conf.regions.t*conf.cells.tot + tInd + (cellPos - 1)*conf.regions.t;
                    else
                        cellTemp.(t) = (1:conf.regions.t) + (tInd - 1)*conf.regions.t*conf.cells.tot + 1 + mult*(tInd - 1) + (cellPos - 1)*conf.regions.t;
                    end
                    
                    % Add one to z values if bottom cooling
                    if bottomCool; cellTemp.z = cellTemp.z + 1; end
                    
                    % Get all volumes from the x,y,z values
                    volumes = combvec(cellTemp.x, cellTemp.y, cellTemp.z)';

                    % Assign cell id to all relevant volumes
                    thermHeatInds(volumes(:,1),volumes(:,2),volumes(:,3)) = cellID;
                end
            end
            % Increment module counter
            m = m + 1; 
        end
    end
end
% Remove zeros from heat inds to avoid indexing error
thermHeatIndsNo0 = thermHeatInds;
thermHeatIndsNo0(thermHeatInds == 0) = 1;
s = size(thermHeatIndsNo0); thermHeatIndsNo0 = reshape(thermHeatIndsNo0,s(1)*s(2)*s(3),[]); % Reshape cell positions to 1D for saturation use
thermCellVols = (thermHeatInds ~= 0); % Volumes that contain cells to remove aluminium heat generation
% Use previous matrix to get indices for mapping from thermal to electrical
thermTempInds = zeros(conf.regions.l*conf.regions.w*conf.regions.t,1,cellsS*cellsP); % Empty matrix
for i = 1:cellsS*cellsP; thermTempInds(:,1,i) =  find(thermHeatInds == i)'; end
% Multiplier for converting heat from cell to volume level
thermHeatCellToVol = 1/(conf.regions.l*conf.regions.w*conf.regions.t);
% Initial temperature condition
thermInitTemp = tempAmb*ones(thermDims)+ conf.bound.initBias; % Currently just ambient temp [degC]
% thermInitTemp = tempAmb - 20 + 40*rand(thermDims); % Currently just ambient temp [degC]
tempTarget = 25; % Battery temp target [degC]

% Calculate max coolant flow rate 
% This is so the flow through one volume at the timestep is less than the
rhoCool = matDict.Density(3);
maxFlowX = min((V*rhoCool)./(conf.time.step.*mdotDistX),[],'all');
maxFlowX(maxFlowX == inf) = [];
maxFlowY = min((V*rhoCool)./(conf.time.step.*mdotDistY),[],'all');
maxFlowY(maxFlowY == inf) = [];
maxFlow = min([maxFlowX maxFlowY]);
if isempty(maxFlow); maxFlow = 0; end 

%% Clear extra variables and define init variable
dataPath = [pwd '\Data\'];
load([dataPath 'keepVars.mat']);
clearvars('-except',keepVars{:});
init = true;

%% Model initialisation function
function [conf, ok] = GetInputs(s,inputLayoutS,inputThermS,inputBoundSimS)
    % Get model configuration from user and store in a config structure
    % Second output is true if function is completed succesfully
    ok = true;
    if ~s % Skip this part if called by a script
        % Variables for cycle input
        inputOptions = ['1 L_NEDC' newline '2 L_NEDC_AB' newline...
                 '3 L_NEDC_H' newline '4 L_NEDC_REV' newline '5 WLTC_1' newline...
                 '6 WLTC_2' newline '7 WLTC_3' newline '8 Zero vel'];
        % Get user inputs - split as there are too many options for one inputdlg
        inputLayout = string(inputdlg({['Cells per module' newline 'Series cells - per module [#]']...
            'Parallel cells - per module [#]' ['Modules in pack x and y direction' newline 'Module count x [#]']...
            'Module count y [#]' ['Cells regions for FVM in cell directions' newline 'Cell regions length [#]']...
            'Cell regions width [#]' 'Cell regions thickness [#]' ['Orientation of cells/modules in pack' newline...
            'Cell length dimension [x,y,z]'] 'Cell width dimension [x,y,z]' 'Include crash structure [y/n]'},...
            'Base pack layout', [1 50], {'4' '3' '9' '3' '5' '3' '1' 'y' 'z' 'y'}));
        inputTherm = string(inputdlg({['Ambient Temperature [' char(176) 'C]'] ['Cooling' newline...
            'Use unlisted input for no cooling e.g. none or no' newline 'Top cooling? [+/-xy]']...
            'Bottom cooling? [+/-xy]' ['Boundary conductivities' newline '0 min - 1 max' newline...
            'Cell-cell [#]'] 'Cell length-module [#]' 'Cell width-module [#]' 'Cell thickness-module [#]'...
            'Module-crash struct/module [#]' ['Init temperature bias [' char(176) 'C]']}, 'Thermal config', [1 50], {'25' 'none' '+y' '1' '0' '1' '1' '0' '0'}));
        inputBoundSim = string(inputdlg({['Boundaries - leave as "nan" for insulated' newline 'x start boundary ['...
            char(176) 'C]'] ['x end boundary [' char(176) 'C]'] ['y start boundary [' char(176) 'C]'] ['y end boundary ['...
            char(176) 'C]'] ['z start boundary [' char(176) 'C]'] ['z end boundary [' char(176) 'C]'] ['Cycle selection [#]'...
            newline inputOptions] ['Simulation options' newline 'Cycle count [#]'] 'Model time step - high values lead to errors [s]'...
            'Logging sample time multiple - rounded up [#]'},...
            'Boundaries and sim', [1 50], {'nan' 'nan' 'nan' 'nan' 'nan' 'nan' '7' '1' '0.1' '10'}));  
    else
        inputLayout = inputLayoutS;
        inputTherm = inputThermS;
        inputBoundSim = inputBoundSimS;
    end
    % Check inputs
    % Attempt to convert string to number
    disp('Inputs collected, checking...');
    conf = struct; % Structure for pack configuration
    % Layout inputs
    [conf.cells.s, cellSeriesCountStatus] = str2num(inputLayout(1)); % series cells
    [conf.cells.p, cellParallelCountStatus] = str2num(inputLayout(2)); % parallel cells
    [conf.modules.x, moduleCountXStatus] = str2num(inputLayout(3)); % x modules
    [conf.modules.y, moduleCountYStatus] = str2num(inputLayout(4)); % y modules
    conf.modules.z = 1; % Only one module in z direction
    [conf.regions.l, cellRegLStatus] = str2num(inputLayout(5)); % length modules
    [conf.regions.w, cellRegWStatus] = str2num(inputLayout(6)); % width modules
    [conf.regions.t, cellRegTStatus] = str2num(inputLayout(7)); % thickness modules
    conf.orient.l = char(inputLayout(8)); % length direction
    conf.orient.w = char(inputLayout(9)); % width direction
    conf.crash = char(inputLayout(10));
    % Thermal inputs
    [conf.bound.amb, ambStatus] = str2num(inputTherm(1)); % ambient temp
    coolingT = char(inputTherm(2)); % top cooling temp variable
    coolingB = char(inputTherm(3)); % bottom cooling temp variable
    [conf.cond.cellCell, cellCellStatus] = str2num(inputTherm(4)); % cell-cell boundary cond
    [conf.cond.cellLenMod, cellLenModStatus] = str2num(inputTherm(5)); % cell length-module boundary cond
    [conf.cond.cellWidMod, cellWidModStatus] = str2num(inputTherm(6)); % cell width-module boundary cond
    [conf.cond.cellThickMod, cellThickModStatus] = str2num(inputTherm(7)); % cell thickness-module boundary cond
    [conf.cond.modCrash, modCrashStatus] = str2num(inputTherm(8)); % module-crash structure boundary cond
    [conf.bound.initBias, initBiasStatus] = str2num(inputTherm(9)); % Init temperature bias
    % Boundaries and sim inputs
    [conf.bound.xa, xaStatus] = str2num(inputBoundSim(1)); % start x bound
    [conf.bound.xb, xbStatus] = str2num(inputBoundSim(2)); % end x bound
    [conf.bound.ya, yaStatus] = str2num(inputBoundSim(3)); % start y bound
    [conf.bound.yb, ybStatus] = str2num(inputBoundSim(4)); % end y bound
    [conf.bound.za, zaStatus] = str2num(inputBoundSim(5)); % start z bound
    [conf.bound.zb, zbStatus] = str2num(inputBoundSim(6)); % end z bound
    [conf.cycle.num, cycleCountStatus] = str2num(inputBoundSim(8)); % cycle count
    [conf.time.step, stepStatus] = str2num(inputBoundSim(9)); % time step
    [conf.time.sample, sampleStatus] = str2num(inputBoundSim(10)); % logging sample time 
    
    % Check the status of conversions and if inputs are in correct format
    % Layout input
    % Cell counts
    if ~cellSeriesCountStatus || conf.cells.s < 1 || rem(conf.cells.s,1) ~= 0
        disp('Cell series count must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    if ~cellParallelCountStatus || conf.cells.p < 1 || rem(conf.cells.p,1) ~= 0
        disp('Cell parallel count must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    conf.cells.tot = conf.cells.s*conf.cells.p; % Total cells in module
    % Module counts
    if ~moduleCountXStatus || conf.modules.x < 1 || rem(conf.modules.x,1) ~= 0
        disp('Module count x must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    if ~moduleCountYStatus || conf.modules.y < 1 || rem(conf.modules.y,1) ~= 0
        disp('Module count y must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    % Cell regions
    if ~cellRegLStatus || conf.regions.l < 1 || rem(conf.regions.l,1) ~= 0
        disp('Cell regions x must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    if ~cellRegWStatus || conf.regions.w < 1 || rem(conf.regions.w,1) ~= 0
        disp('Cell regions y must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    if ~cellRegTStatus || conf.regions.t < 1 || rem(conf.regions.t,1) ~= 0
        disp('Cell regions z must be an integer more than or equal to 1, aborting...')
        ok = false;
        return;
    end
    % Cell orientation
    if conf.orient.l ~= 'x' && conf.orient.l ~= 'y' && conf.orient.l ~= 'z'
        disp('Cell length dimension must be x, y or z, aborting...')
        ok = false;
        return;
    end
    if conf.orient.w ~= 'x' && conf.orient.w ~= 'y' && conf.orient.w ~= 'z'
        disp('Cell width dimension must be x, y or z, aborting...')
        ok = false;
        return;
    end
    % Assign third dimension
    if ~contains([conf.orient.l conf.orient.w],'x')
        conf.orient.t = 'x'; 
    elseif ~contains([conf.orient.l conf.orient.w],'y')
        conf.orient.t = 'y'; 
    else
        conf.orient.t = 'z';
    end
    
    % Thermal inputs
    % Ambient temperature
    if ~ambStatus || isnan(conf.bound.amb)
        disp('Ambient temp must be a number, aborting...')
        ok = false;
        return;
    end
    % Cooling
    switch coolingT
        case '+x' % Positive x direction top cooling
            conf.cooling.t.x = 1; % Mass flow multiplier x
            conf.cooling.t.y = 0; % Mass flow multiplier y
        case '-x' % Positive x direction top cooling
            conf.cooling.t.x = -1; % Mass flow multiplier x
            conf.cooling.t.y = 0; % Mass flow multiplier y
        case '+y' % Positive x direction top cooling
            conf.cooling.t.x = 0; % Mass flow multiplier x
            conf.cooling.t.y = 1; % Mass flow multiplier y
        case '-y' % Positive x direction top cooling
            conf.cooling.t.x = 0; % Mass flow multiplier x
            conf.cooling.t.y = -1; % Mass flow multiplier y
        otherwise
            conf.cooling.t.x = 0; % Mass flow multiplier x
            conf.cooling.t.y = 0; % Mass flow multiplier y
    end
    switch coolingB
        case '+x' % Positive x direction top cooling
            conf.cooling.b.x = 1; % Mass flow multiplier x
            conf.cooling.b.y = 0; % Mass flow multiplier y
        case '-x' % Positive x direction top cooling
            conf.cooling.b.x = -1; % Mass flow multiplier x
            conf.cooling.b.y = 0; % Mass flow multiplier y
        case '+y' % Positive x direction top cooling
            conf.cooling.b.x = 0; % Mass flow multiplier x
            conf.cooling.b.y = 1; % Mass flow multiplier y
        case '-y' % Positive x direction top cooling
            conf.cooling.b.x = 0; % Mass flow multiplier x
            conf.cooling.b.y = -1; % Mass flow multiplier y
        otherwise
            conf.cooling.b.x = 0; % Mass flow multiplier x
            conf.cooling.b.y = 0; % Mass flow multiplier y
    end 
    % Boundary conductivities
    if ~cellCellStatus || conf.cond.cellCell > 1 || conf.cond.cellCell < 0
       disp('Cell-cell conductivity incorrect input, must be between 0 and 1, aborting...');  
       ok = false;
       return;
    end
    if ~cellLenModStatus || conf.cond.cellLenMod > 1 || conf.cond.cellLenMod < 0
       disp('Cell length-module conductivity incorrect input, must be between 0 and 1, aborting...');  
       ok = false;
       return;
    end
    if ~cellWidModStatus || conf.cond.cellWidMod > 1 || conf.cond.cellWidMod < 0
       disp('Cell widh-module conductivity incorrect input, must be between 0 and 1, aborting...');  
       ok = false;
       return;
    end
    if ~cellThickModStatus || conf.cond.cellThickMod > 1 || conf.cond.cellThickMod < 0
       disp('Cell thickness-module conductivity incorrect input, must be between 0 and 1, aborting...');  
       ok = false;
       return;
    end
    if ~modCrashStatus || conf.cond.modCrash > 1 || conf.cond.modCrash < 0
       disp('Module-crash structure conductivity incorrect input, must be between 0 and 1, aborting...');  
       ok = false;
       return;
    end
    % Init bias
    if ~initBiasStatus
       disp('Init temperature bias input incorrect, must be a number, aborting...');  
       ok = false;
       return;
    end
    
    % Boundary and sim inputs
    % Boundary conditions
    if ~xaStatus || ~xbStatus || ~yaStatus || ~ybStatus || ~zaStatus || ~zbStatus 
        disp('Boundary conditions must be a number or nan, aborting...')
        ok = false;
        return;
    end
    % Cycle data
    dataPath = [pwd '\Data\'];
    load([dataPath 'combined_cycle_data.mat']);
    switch inputBoundSim(7) % cycle switch case
        case "1"
            cycle = L_NEDC;
        case "2"
            cycle = L_NEDC_AB;
        case "3"
            cycle = L_NEDC_H;
        case "4"
            cycle = L_NEDC_REV;
        case "5"
            cycle = WLTC_1;
        case "6"
            cycle = WLTC_2;
        case "7"
            cycle = WLTC_3;
        case "8"
            cycle(:,1) = (0:1000)';
            cycle(:,2) = zeros(1001,1);
        otherwise
            disp('Invalid cycle input, see input prompt for options, aborting...');
            ok = false;
            return;
    end
    if ~cycleCountStatus || conf.cycle.num < 0
        disp('Cycle number input must be more than 0, aborting...');
        ok = false;
        return;
    end
    % Resise cycle if needed
    if conf.cycle.num < 1 % Need to trim the cycle
        cycle(round(conf.cycle.num*height(cycle)):end,:) = [];
    elseif conf.cycle.num > 1 % Need to add more to the cycle
        step = cycle(2,1) - cycle(1,1); % Get time step, must be consistent
        integ = floor(conf.cycle.num); fract = conf.cycle.num - floor(conf.cycle.num);
        cycle = [repmat(cycle,integ,1); cycle(1:round(fract*height(cycle)),:)];
        cycle(:,1) = (cycle(1,1):step:step*(length(cycle(:,1)) - 1))';
    end
    conf.cycle.time = cycle(:,1);
    conf.cycle.vel = cycle(:,2);
    conf.cycle.stop = cycle(end,1);
    if ~stepStatus || conf.time.step < 0
        disp('Time step input must be more than 0, aborting...');
        ok = false;
        return;
    end
    if ~sampleStatus || conf.time.sample < 0
        disp('Sample time input must be more than 0 or equal to -1, aborting...');
        ok = false;
        return;
    end
    conf.time.sample = ceil(conf.time.sample)*conf.time.step;        
     
    disp('Inputs ok, continuing...');
end