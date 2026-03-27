%% This script was prepared by: Ally Vizcarra : Created March 19 2026
%                               University of Wisconsin-Madison

% Welcome! This is the code that documents how to run the linear inverse model 
% (LIM) and cyclostationary LIM (CSLIM). In addition, you will also learn
% how to run a simple forecast using these models.

% Clear your workspace and terminal of any variables before continuing
clear; clc; 

% Here is where you put your own directory of where to store the files you
% need in order to run this code. To change the directory, you can rewrite
% the string in 'dirin' as eg: dirin = '/Users/folder/where/you/keep/your/files';
dirin = '/Users/ally/Documents/Research/Data/Dan_Paper';

% The variables that we will put in the state vector are monthly sea surface
% temperature (SST) from ERSSTv5, and monthly sea surface height (SSH) from ORAS5
% from the time period 1958-2022 (780 months after getting rid of leap years) 

% From the .nc file, load in the principle components (PCs) of:
pcs_sst = ncread(fullfile(dirin,'state_vector_data.nc'),'pcs_sst');  % SST: leading 11   
size(pcs_sst) % this should be of size 780 x 11 (years x SST pcs)
pcs_ssh = ncread(fullfile(dirin,'state_vector_data.nc'),'pcs_ssh');  % SSH: leading 5
size(pcs_ssh) % this should be of size 780 x 5 (years x SSH pcs)
% These PCs represent the temporal evolution of the dominant modes

% NOTE: the number of leading PCs we chose are already pre-determined using
% a tau-test and cannot be changed from the .nc file. 
% In addition, the PCs were calculated using the bounds: 30N-30S latitude 
% and 100E-285E longitude, and also cannot be changed from the .nc file.
% So, the region we will be looking at is the Tropical Pacific! 

% Now we load in the scaling factor for the PCs. It contians the eigenvalues 
% for the SST and SSH PCs.
scaling = ncread(fullfile(dirin,'state_vector_data.nc'),'scaling');
scale_sst = scaling(1); % eigenvalue for SST
disp(['SST eigenvalue: ', num2str(scale_sst)]);
scale_ssh = scaling(2); % eigenvalue for SSH
disp(['SSH eigenvalue: ', num2str(scale_ssh)]);

% Now we build the state vector that's used to run the models. We use the
% PCs for both SST and SSH, and then scale the PCs by their eigenvalues
x = [pcs_sst/scale_sst pcs_ssh/scale_ssh]; 
size(x) % this should be of size 780 x 16 (time x total # of PCs)

% Here we will name the the values in the state vector. The first value is
% the total number of months (780), so its name is 'ntime'
% The second value is the total number of PCs (16), so its name is 'mkp'
[ntime, mkp] = size(x); 

% Congrats! Now we can construct the LIM and the CSLIM in the next section

%% Constructing the LIM

% Define a time lag between the covariance matrices. This value is in months
lag_use = 3; 

% This is calculating the zero-lag covariance matrix
C0 = x'*x / (ntime-1);

% Here we create lagged time indices
time_init = (1:(ntime - lag_use))';  % this represents x(t), an initial time 't'
time_final = ((lag_use+1):ntime)';   % this represents x(t+tau), a final time 't+tau'

% This is calculating tau-lag covariance matrix
Ct = (x(time_final,:)'*x(time_init,:))/(ntime - lag_use);

% The variable 'G' is called the propagator matrix that is being estimated
% using the zero-lag and tau-lag covariance matrices
G_LIM = Ct/C0;

% Then, we convert 'G' into the continuous dynamical operator 'L' 
% using matrix logarithms
L_LIM = logm(squeeze(G_LIM))/lag_use; 

% Here we are calculating estimating the stochastic forcing covariance
% matrix 'Q' using the Fluctuation-Dissipation Relation
Q_LIM = -L_LIM*C0 - C0*L_LIM';

disp("LIM Construction Complete");
% Thus the LIM has been constructed!

%% Constructing the CSLIM
% The construction of the CSLIM is different compared to the LIM. To put it
% shortly before you continue on to the rest of the code: the LIM does not
% take into account of seasonality, so it only has one 'G' based on an
% initial time and a final time. The CSLIM corrects that by calculating
% twelve 'G's for the calendar year, and is able to account for
% seasonality.

% Recall these variables we used before for the LIM, we are just using them
% again for the CSLIM
[ntime, mkp] = size(x); 
lag_use = 3;

nyear = ntime/12 ; 
% Januaries indices for each year
time_ind = 13:12:(ntime-12); 

% In this case the CSLIM is calculated after applying a three month running
% mean filter to the zero-lag and 1-lag covariance matrices.

% Here we are calculating the one-lag (one month) and zero-lag (zero month)
% covariance matrices
C1_out = nan(12, mkp, mkp); % one-lag covariance
C0_out = nan(12, mkp, mkp); % zero-lag covariance

% The covariance matrices are calculated for every month in the calendar
% year
for imon = 1:12  
    C0_out(imon,:,:) = x(time_ind+imon-1,:)'*x(time_ind+imon-1,:)/(nyear - lag_use);
    C1_out(imon,:,:) = x(time_ind+imon,:)'*x(time_ind+imon-1,:)/(nyear - lag_use); 
end

% A three-month "circular" filter is applied to the covariance matrices
C0_out = (C0_out([12 1:11],:,:) + C0_out + C0_out([2:12 1],:,:)) / lag_use; 
C1_out = (C1_out([12 1:11],:,:) + C1_out + C1_out([2:12 1],:,:)) / lag_use; 

% Now we compute 'G' and 'L'. First we initialize, and we note that they are
% also being computed per month
G_CSLIM = nan(12, mkp, mkp);
L_CSLIM = nan(12, mkp, mkp); 

% NOTE: 'x' is centered on mid-month values (monthly mean).
% The zero-lag covariance matrix is thus centered on mid-month, and the one-lag
% covariance matrix is mid-month to mid-month, meaning that 'G' takes
% the mid-month state to the next mid-month state. 'L' is representative of
% the dynamics that takes you from one mid-month to another mid-month. You
% can think of this as being representative of the end of one month through
% the beginning of the next month, but that may not be the best way to
% think of it. This will matter, though, when we try to estimate Q.
for imon = 1:12  
    C0 = squeeze(C0_out(imon,:,:)); 
    C1 = squeeze(C1_out(imon,:,:)); 
    G_CSLIM(imon,:,:) = C1/C0; 
    L_CSLIM(imon,:,:) = logm(squeeze(G_CSLIM(imon,:,:)))/1; 
end

% Recall that 'Q' is the stochastic forcing covariance matrix. It is also
% calculated per month
Q_CSLIM = nan(12, mkp, mkp); % This is the raw estimate of 'Q'
Q_adjusted = Q_CSLIM;        % This 'Q' is adjusted by eigenvalues
Q_adjusted_SPD = Q_CSLIM;    % Nearest semi-positive definite (SPD) matrix

% Variable to store negative eigenvalues
neg_eig_vals = nan(12, 1); 

% Now, we calculate 'Q' using the Fluctuation-Dissipation Relation. 
% There will be some negative eigenvalues ('neg_eig_vals') so we have to 
% track those. Then we will "fix" 'Q' using the nearest semi-positive definite matrix.
for imon = 1:12  
    % Recall that the CSLIM is taking account of seasonal cycle, so to take
    % into account of the calendar months, we have to be sure that the
    % months are of modulo 12
    month_next = mod(imon,12)+1;            % Next month index
    month_current = mod(imon-1,12)+1;       % This month index 
    month_prev = mod(imon-2,12)+1;          % Last month index

    % Now, calculate 'Q' using centered differencing. 
    % NOTE: 'L' is representative of "end of month", while all the "C0"s are
    % representative of "mid-month" values. We will calculate "Q" to be
    % consistent with "L" in that it is the noise that takes us from one
    % mid-month to another mid-month
    C_next = squeeze(C0_out(month_next,:,:));  % Next month
    C_current = squeeze(C0_out(month_current,:,:));  % THIS month
    
    % The next line centers 'C0' so that it is between 'month_next' and 
    % 'month_current', and the centered differencing does not introduce 
    % a phase lag. This also means that 'Q' will be "centered" at the end 
    % of the month, like 'L'. 
    C0 = squeeze(C0_out(month_current,:,:) + C0_out(month_next,:,:))/2; 
    Q_CSLIM(imon,:,:) = (C_next - C_current)/1 - (squeeze(L_CSLIM(imon,:,:))*C0 + C0*squeeze(L_CSLIM(imon,:,:))'); 

    % Now we check for negative eigenvalues
    [v, d] = eig(squeeze(Q_CSLIM(imon,:,:))); 
    % This counts how many eigenvalues are negative
    neg_eig_vals(imon) = sum(diag(d)<0);

    % The next two lines basically set the negative eigenvalues to zero,
    % and rescales the eigenvalues to preserve the total variance (trace)
    trd = trace(d); 
    d(d < 0) = 0 ; 
    d = d*trd/trace(d); 
    Q_adjusted(imon,:,:) = v*d*v'; 

    % Now, adjust using nearest semi positive definite matrix
    % NOTE: You will need the nearestSPD tool. If you don't have this
    % already, it will run an error here. Use the link below to download
    % this into your toolbox
    % https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd
    Q_adjusted_SPD(imon,:,:) = nearestSPD(squeeze(Q_CSLIM(imon,:,:))); 
end

% Congratulations! The CSLIM has also been constructed!

%% This section is how to run a simple forecast

% Now that we know how the models are constructed, let us then go through
% how to make a simple forecast. For simplicity, let's just focus on
% deterministic (predictable) forecasts (no noise).

% So the goal for this part of the script is to:
% 1. Define an initial state x(t)
% 2. Evolve that state forward in time using the LIM and CSLIM
% 3. Compare how the two models differ in their predictions by plotting
%    SST spatial maps 

% In order to get those spatial maps, we need the empirical orthogonal
% functions (EOFs) maps for SST - we'll come back to these later. So let's
% just load them in at the beginning.
% A short note on what EOFs are: they provide spatial patterns 
% corresponding to each PC.

% As a reminder, the latitude bounds are 30N-30S, and the longitude 
% bounds are from 100E-285E.

% Load in latitude ('lat'), longitude ('lon'), and the EOF SST map ('EOF_SST')
lat = ncread(fullfile(dirin,'EOFs.nc'),'lat');
lon = ncread(fullfile(dirin,'EOFs.nc'),'lon');
eof_sst = ncread(fullfile(dirin,'EOFs.nc'),'EOF_SST');

% We are just defining the length of the lats and lons
nlat = length(lat);
nlon = length(lon);

% Now let's first define an initial state. Recall that 'x' is of size 
% 780x16, which is time (month) x PCs. So the initial state x(1,:) 
% means that the first time point (January 1958) is the starting point 
% for the models to run.
% Feel free to change this value!
x0 = x(780,:)';

% The models will be run forward at different monthly lags. Let's set the
% max lag to be at 18 months. 
lead_max = 18;

%% LIM Forecast

% Here, the LIM assumes: x(t+1) = G_LIM * x(t)
% where x(t+1) is the final state, x(t) is the initial state,
% and G_LIM is "constant" in time 

% Initializing the LIM forecast
x_forecast_LIM = zeros(mkp,lead_max+1);
x_forecast_LIM(:,1) = x0;

% Forward integration - now we are making the forecast at different lags
% (0-18)
for ilead = 1:lead_max
    x_forecast_LIM(:,ilead+1) = G_LIM * x_forecast_LIM(:,ilead);
end

%% CSLIM Forecast

% Here the CSLIM assumes: x(t+1) = G(month) * x(t)
% where x(t+1) is the final state, x(t) is the initial state,
% and G varies with calendar month.

% Initializing the CSLIM forecast
x_forecast_CSLIM = zeros(mkp,lead_max+1);
x_forecast_CSLIM(:,1) = x0;

% Defining the starting month: January
start_month = 1;

% Forward integration - now we are making the forecast
for ilead = 1:lead_max
    % Cycle through months (1–12)
    month = mod(start_month + ilead -2,12) + 1;

    % Extract the twelve 'G's
    G = squeeze(G_CSLIM(month,:,:));

    % Propagate state forward
    x_forecast_CSLIM(:,ilead+1) = G * x_forecast_CSLIM(:,ilead);
end

%% Below we are separating out the SST PCs and the SSH PCs

sst_pc_LIM = x_forecast_LIM(1:11,:) * scale_sst;
ssh_pc_LIM = x_forecast_LIM(12:16,:) * scale_ssh;

sst_pc_CSLIM = x_forecast_CSLIM(1:11,:) * scale_sst;
ssh_pc_CSLIM = x_forecast_CSLIM(12:16,:) * scale_ssh;

sst_map_LIM = zeros(nlat,nlon,lead_max+1);
sst_map_CSLIM = zeros(nlat,nlon,lead_max+1);

% And then we are reconstructing the SST fields back into physical space
% The reason why we only convert SST is because it best shows ENSO
% signature. But if you were curious on the SSH field, feel free to do that
% on your own- but just as a note you would need to compute the EOFs for
% SSH.

% EOF expansion: SST(lat,lon,t) = sum_k PC_k(t) * EOF_k(lat,lon)
for ilead = 1:lead_max+1
    for mode = 1:11
        sst_map_LIM(:,:,ilead) = sst_map_LIM(:,:,ilead) + ...
            eof_sst(:,:,mode) * sst_pc_LIM(mode,ilead);

        sst_map_CSLIM(:,:,ilead) = sst_map_CSLIM(:,:,ilead) + ...
            eof_sst(:,:,mode) * sst_pc_CSLIM(mode,ilead);
    end
end

disp('Full Deterministic LIM and CSLIM forecast complete.')

%% Plot the Tropical Pacific at different lead times
% This is a nice visual of showcasing the differences between the models at
% different lead times.

% This is an array where you can change the lead times 
% NOTE: inputs have to be >= 1 and <= 19
% Eg: 'lead_to_plot = [1]' would be equate to zero-month lead, but you can
% also look at multiple plots at the same time by changing (eg:)
% 'lead_to_plot' = [1 3] to look at zero and two-month lead
% times, and plots will be generated for both
% As a reccomendation, I would keep the zero-month lag for comparison
% because it represents the "current" state of the system.
lead_to_plot = [1 19]; 

% In this for loop you can view the model differences from 'lead_to_plot'
% These plots will show the SST (C) anomaly field
for i = 1:length(lead_to_plot)

    % Lead time in 'lead_to_plot'
    lead = lead_to_plot(i);

    % Produce figure
    figure('Position',[200 200 1000 500]);

    % We are making three different subplots 
    % The first is the LIM forecast 
    subplot(1,3,1)  
    imagesc(lon,lat,sst_map_LIM(:,:,lead))  
    set(gca,'ydir','normal')  
    cb = colorbar;  
    xlabel('Longitude')  
    ylabel('Latitude')  
    ylabel(cb,'SST Anomaly (°C)')  
    title(['LIM Lead ' num2str(lead-1) ' months'])  
    
    % The second is the CSLIM forecast  
    subplot(1,3,2)  
    imagesc(lon,lat,sst_map_CSLIM(:,:,lead))  
    set(gca,'ydir','normal')  
    cb = colorbar;  
    xlabel('Longitude')  
    ylabel('Latitude')  
    ylabel(cb,'SST Anomaly (°C)')  
    title(['CSLIM Lead ' num2str(lead-1) ' months'])  
    
    % The third is taking the difference  
    % The plot should look like a solid color for a zero-month lead since
    % it should be showing what the state of the system looks like
    % "currently"
    subplot(1,3,3)  
    imagesc(lon,lat,sst_map_CSLIM(:,:,lead) - sst_map_LIM(:,:,lead))  
    set(gca,'ydir','normal')  
    cb = colorbar;  
    xlabel('Longitude')  
    ylabel('Latitude')  
    ylabel(cb,'ΔSST (°C)')  
    title('CSLIM - LIM')  
end

% By changing when the initial state 'x(0)' and viewing the different lags,
% you can clearly see the difference on how far in advance the LIM and CSLIM 
% are able to capture an El Niño (or La Niña structure).
