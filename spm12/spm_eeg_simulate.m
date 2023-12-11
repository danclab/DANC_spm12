function [Dnew, meshsourceind] = spm_eeg_simulate(D, varargin)

useind=1; % D to use

% Create an inputParser instance
p = inputParser;

% Define the required argument
addRequired(p, 'D');

% Define the optional arguments with default values
addOptional(p, 'prefix', 'sim');
addOptional(p, 'patchmni', []);
addOptional(p, 'simsignal', []);
addOptional(p, 'ormni', []);
addOptional(p, 'woi', []);  % Default will be set later based on D
addOptional(p, 'whitenoise', []);
addOptional(p, 'SNRdB', []);
addOptional(p, 'trialind', []);
addOptional(p, 'dipfwhm', 6);
addOptional(p, 'nAmdipmom', 1);
addOptional(p, 'noiseD', []);

% Parse the input arguments
parse(p, D, varargin{:});

% Assign the parsed values to variables
prefix = p.Results.prefix;
patchmni = p.Results.patchmni;
simsignal = p.Results.simsignal;
ormni = p.Results.ormni;
woi = p.Results.woi;
whitenoise = p.Results.whitenoise;
SNRdB = p.Results.SNRdB;
trialind = p.Results.trialind;
dipfwhm = p.Results.dipfwhm;
nAmdipmom = p.Results.nAmdipmom;
noiseD = p.Results.noiseD;

% Set the default value for woi if it was not provided
if isempty(woi)
    woi = [D{useind}.time(1) D{useind}.time(end)];
end

val = D{useind}.val;

% Simplified logical check for mutual exclusivity
if ~isempty(whitenoise) + ~isempty(SNRdB) + ~isempty(noiseD) ~= 1
    error('Must specify exactly one of white noise level, sensor level SNR, or noiseD');
end

% Efficient file parts handling
[~, b1, ~] = fileparts(D{useind}.fname);
newfilename = [prefix b1 '.mat'];

% Forcing overwrite of an existing file
Dnew = D{useind}.clone(newfilename);

% Default value assignment
if isempty(trialind)
    trialind=1:Dnew.ntrials;
end

modstr = deblank(modality(D{1}));
disp(['Simulating data on ' modstr ' channels only']);

% Source and signal size consistency check
if ~isempty(patchmni),
    Ndips=size(patchmni,1);
else
    Ndips=0;
end
if size(simsignal, 1) ~= Ndips
    error('Number of signals given does not match number of sources');
end

% Optimized distance calculation
meshsourceind = zeros(1, Ndips);
mnidist = zeros(1, Ndips);
for d = 1:Ndips
    distances = sqrt(sum(bsxfun(@minus, Dnew.inv{val}.mesh.tess_mni.vert, patchmni(d, :)).^2, 2));
    [mnidist(d), meshsourceind(d)] = min(distances);
end

disp(['Furthest distance from dipole location to mesh ' num2str(max(mnidist), '%3.2f') ' mm']);
if max(mnidist) > 0.1
    warning('Supplied vertices do not sit on the mesh!');
end

try
    chanind = Dnew.indchantype({'MEG', 'MEGPLANAR'}, 'GOOD');
catch
    chanind = Dnew.indchantype(modality(D{1}), 'GOOD');
end


sensorunits = Dnew.units; % Units of sensors (T or fT)
Ntrials = Dnew.ntrials;  % Number of trials

% Define period over which dipoles are active
startf1 = woi(1);  % Start time (sec)
endf1 = woi(2);    % End time
f1ind = find(Dnew.time > startf1 & Dnew.time <= endf1);

if length(f1ind) ~= size(simsignal, 2)
    error('Signal does not fit in time window');
end


% Scale calculation for sensor units
chan_idx=Dnew.indchantype('MEG');
switch sensorunits{chan_idx(1)}
    case 'T'
        simscale = 1e-15;  % Convert from fT to T
    case 'fT'
        simscale = 1.0;    % Sensors already in fT
    case 'uV'
        % Do nothing, already in micro volts
    case 'V'
        simscale = 1e-6;   % Convert from V to uV
        error('EEG units (V) not supported at the moment');
    otherwise
        error('Unknown sensor unit');
end

% Apply scale to whitenoise if it's not empty
if ~isempty(whitenoise)
    whitenoise = whitenoise * simscale;
end


if ~isempty(ormni) %%%% DIPOLE SIMULATION
    disp('SIMULATING DIPOLE SOURCES');

    if size(ormni) ~= size(patchmni)
        error('A 3D orientation must be specified for each source location');
    end

    % Transform positions and orientations to MEG space
    posdipmm = Dnew.inv{val}.datareg.fromMNI * [patchmni, ones(size(ormni, 1), 1)]';
    posdipmm = posdipmm(1:3, :)';

    M = Dnew.inv{val}.datareg.fromMNI * Dnew.inv{val}.mesh.Affine;
    ordip = [ormni, ones(size(ormni, 1), 1)] * inv(M')';
    ordip = bsxfun(@rdivide, ordip(:, 1:3), sqrt(sum(ordip(:, 1:3).^2, 2)));

    sens = Dnew.inv{val}.forward.sensors;
    vol = Dnew.inv{val}.forward.vol;
    sensorind = Dnew.indchantype({'MEG', 'MEGPLANAR', 'REFMAG', 'REFGRAD'});
    [vol, sens] = ft_prepare_vol_sens(vol, sens);

    % Preallocate tmp
    if length(size(simsignal)) == 2
        tmp = zeros(length(chanind), Dnew.nsamples);
    else
        tmp = zeros(length(chanind), Dnew.nsamples, Dnew.ntrials);
    end

    % Compute lead fields and gains
    for i = 1:Ndips
        gmn = ft_compute_leadfield(posdipmm(i, :) * 1e-3, sens, vol, 'dipoleunit', 'nA*m', 'chanunit', sensorunits(sensorind));
        gain = gmn(chanind-3, :) * ordip(i, :)' * nAmdipmom(i);

        if length(size(simsignal)) == 2
            tmp(:, f1ind) = tmp(:, f1ind) + gain * simsignal(i, :);
        else
            for j = 1:Dnew.ntrials
                tmp(:, f1ind, j) = tmp(:, f1ind, j) + gain * squeeze(simsignal(i, :, j));
            end
        end
    end



else %%% CURRENT DENSITY ON SURFACE SIMULATION
    disp('SIMULATING CURRENT DISTRIBUTIONS ON MESH');
    fprintf('Computing Gain Matrix: ');

    % Gain matrix computation
    [L, Dnew] = spm_eeg_lgainmat(Dnew);              
    if isfield(Dnew.inv{val}.forward, 'scale')
        L = L ./ Dnew.inv{val}.forward.scale; % Rescale lead fields
    end

    Nchans = size(L, 1);

    fprintf(' - done\n');

    nativemesh = Dnew.inv{val}.forward.mesh;
    Qe = []; % SNR may be defined by sensor level data

    base.FWHMmm = dipfwhm;
    base.nAm = nAmdipmom;

    % Directory for saving prior
    [~, b1, ~] = fileparts(Dnew.fname);
    priordir = fullfile(Dnew.path, ['simprior_' b1]);
    if ~exist(priordir, 'dir')
        mkdir(priordir);
    end
    fprintf('Saving prior in directory %s\n', priordir);

    [Qp, ~, ~] = spm_eeg_invert_setuppatches(meshsourceind, nativemesh,...
        base, priordir, Qe, L);

    % Preallocate fullsignal and tmp
    fullsignal = zeros(Ndips, Dnew.nsamples); 
    fullsignal(1:Ndips, f1ind) = simsignal;

    tmp = sparse(Nchans, Dnew.nsamples); % Simulated data
    X = zeros(size(full(Qp{1}.q)));

    % Efficient matrix operations for signal computation
    for j = 1:Ndips
        Lq = L * Qp{j}.q; % Lead field * prior source distribution
        X = X + full(Qp{j}.q);
        tmp = tmp + Lq * fullsignal(j, :);
    end

    tmp = tmp .* simscale; % Apply scale to sensor units


end % if ori


% Standard deviation across channels
allchanstd = std(tmp, [], 2);
if ismatrix(simsignal) % For 3D simsignal
    allchanstd = squeeze(mean(allchanstd, 3));
end
meanrmssignal = mean(allchanstd);

% Handling SNR and white noise
if ~isempty(SNRdB)
    whitenoise = meanrmssignal * (10^(-SNRdB / 20));
    disp(['Setting white noise to give sensor level SNR of ', num2str(SNRdB), 'dB']);
elseif ~isempty(noiseD)
    noiseLevel = mean(std(squeeze(mean(noiseD, 3)), [], 2));
    snr = -20 * log10(noiseLevel / meanrmssignal);
    disp(['Noise level=', num2str(noiseLevel, '%.3f'), ' fT RMS, SNR=', num2str(snr, '%.3f'), ' dB']);
end

% Simulating data for each trial
for i = 1:Ntrials
    if ismember(i, trialind) % Add signal to specific trials
        if ismatrix(tmp)
            Dnew(chanind, :, i) = full(tmp);
        else
            Dnew(chanind, :, i) = squeeze(tmp(:, :, i));
        end
    else
        Dnew(chanind, :, i) = zeros(size(tmp, 1), size(tmp, 2));
    end

    % Adding noise
    if isempty(noiseD)
        Dnew(:, :, i) = Dnew(:, :, i) + randn(size(Dnew(:, :, i))) * whitenoise;
    else
        noise_channel = randsample(setdiff(1:noiseD.ntrials, noiseD.badtrials), 1);
        Dnew(chanind, :, i) = Dnew(chanind, :, i) + noiseD(chanind, :, noise_channel);
    end
end




%% Plot and save

Dnew.save;

fprintf('\n Finish\n')

end
