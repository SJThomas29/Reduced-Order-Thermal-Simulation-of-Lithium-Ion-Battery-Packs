%% Model order reduction
%% Define state space model
% Using M for the test data - reshape into a vector for each time step
if ~exist('volTemp','var')
    % Data not included in github repo due to size
end
% Need to define the variables below from a simulation
tempAmb = data.tempAmb;
volTemp = data.volTemp;
volHeatGen = data.volHeatGen.Data;
coolInTemp = data.coolIn.Data;
time = data.Time;

s = size(volTemp);
M = reshape(volTemp,[s(1)*s(2)*s(3) s(4)])';
% Centre M by removing the mean
Mc = M - mean(M,1);

% Perform the PCA, will be no difference in centred or non-centred data as
% the PCA function does this
if ~exist('coeff','var')
    load([pwd '/Data/pca_results.mat']);
%     tic;[coeff,score,latent,tsquared,explained,mu] = pca(Mc,"NumComponents",10);toc; % Uncomment this line to run pca on new data
end

% Calculate error caused by PCA
R = score*coeff'; % Results from reduced temperatures
error = Mc - R;
disp(['Mean temp error from pca reduction is ' num2str(mean(abs(error),'all')) char(176) 'C']);

% Calculate state space form of the order reduction using linear regression
% xdot = Ax + Bu
% Where the states (x) are the temperatures, and inputs (u) are the total
% power generation across the pack and the input coolant temperatures
xdot = Mc([ 2:end end],:) - Mc([1:end-1 end-1],:); % Using centred values - doesnt make a difference to answer though
xdotr = xdot*coeff; % Reduced order xdot values using pca coefficients
totalHeatGen = squeeze(sum(volHeatGen,[1 2 3])); % Total heat generation across pack at each step
u = [totalHeatGen (coolInTemp - tempAmb)];
AB = [score u]\xdotr;
A = AB(1:10,:); % 10 modes considered
B = AB(11:12,:); % Two inputs

%% Test state space model
step = 1; % Time step in data
X = zeros(length(time),10); % Initial conditions
Xres = zeros(s);
Xres(:,:,:,1) = volTemp(:,:,:,1);
for i = 2:length(time)
    u = [totalHeatGen(i-1) (coolInTemp(i-1) - tempAmb)]; % Inputs from last step

    Xdot = u*B + X(i-1,:)*A; % Calculate last step gradient

    X(i,:) = X(i-1,:) + Xdot*step; % Discrete integrator for current step  

    XfullTemp = X(i,:)*coeff'; % Calculate the full data

    XfullTemp = reshape(XfullTemp,[s(1) s(2) s(3)]); % Reshape data
    
    Xres(:,:,:,i) = XfullTemp + tempAmb;
end

% Mean error
Xerror = Xres - volTemp;
disp(['Mean temp error from state space test ' num2str(mean(abs(Xerror),'all')) char(176) 'C']);

