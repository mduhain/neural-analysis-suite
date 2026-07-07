%% NEURODATA WITHOUT BORDERS TEST SCRIPT
%
% Basic script for playing around with this data structure with 2P Ca+ Data
%
% 2022-08-09 (mduhain) - Created script
%
% 2026-07-07 (mduhain) - Updated script with NWB version 2.1.0 core schema
%
% Follow along from:
% https://neurodatawithoutborders.github.io/matnwb/tutorials/html/ophys.html
%
% OVERALL ORGANIZATION PLAN:
% 
% Total of three imaging planes: 
%  1. Surface (vasculature)
%  2. Inhibitory (field, tdTom+ somas)
%  3. Experiment (field, GCaMP activity)
%
%  Each imaging plane has 2 optical channels (green and red)

%% Installing MatNWB (mac / bash)
% !git clone https://github.com/NeurodataWithoutBorders/matnwb.git
% cd matnwb
% addpath(genpath(pwd))
generateCore('2.1.0'); %Generate Matlab classes from NWB core schema files


%% Set up a NWB file
nwb = NwbFile( ...
    'session_description', 'mouse vibrotactile stimulation',...
    'identifier', 'Mouse0_Day0', ...
    'session_start_time', datetime(2022, 8, 9, 3, 37, 3), ...
    'general_experimenter', 'mduhain', ... % optional
    'general_session_id', 'session_0001', ... % optional
    'general_institution', 'University of Rochester' );


%% Subject Information
subject = types.core.Subject( ...
    'subject_id', '000', ...
    'age', 'P90D', ...
    'species', 'Mus musculus', ...
    'sex', 'M', ...
    'genotype', 'PV-Cre-tdTom');
nwb.general_subject = subject;


%% Optical Physiology

% Initialize: Device information
device = types.core.Device( ...
    'description', 'Haptics Lab two-photon microscope', ...
    'manufacturer', 'Sutter', ...
    'model_name', 'MOM-RES');

% Incorperate: Device to NWB
nwb.general_devices.set('Device', device);

% Initialize: Imaging Plane 1 (surface-level vasculature)

optical_channel_green = types.core.OpticalChannel( ...
    'name','GreenChannel',...
    'description', 'Green emission from GCaMP', ...
    'emission_lambda', 510);

optical_channel_red = types.core.OpticalChannel( ...
    'name','RedChannel',...
    'description', 'Red emission from tdTomato', ...
    'emission_lambda', 585);


imaging_plane_name = 'imaging_plane_surface';

imaging_plane = types.core.ImagingPlane( ...
    'optical_channel', optical_channel_green, ...
    'description', 'a very interesting part of the brain', ...
    'device', types.untyped.SoftLink(device), ...
    'excitation_lambda', 600., ...
    'imaging_rate', 5., ...
    'indicator', 'GFP', ...
    'location', 'Primary visual area');

nwb.general_optophysiology.set(imaging_plane_name, imaging_plane);

%% SpatialSeries and Position
position_data = [linspace(0,10,100); linspace(0,8,100)];
 
spatial_series_ts = types.core.SpatialSeries( ...
    'data', position_data, ...
    'reference_frame', '(0,0) is bottom left corner', ...
    'timestamps', linspace(0, 100)/200);

Position = types.core.Position('SpatialSeries', spatial_series_ts);

% create processing module
behavior_mod = types.core.ProcessingModule( ...
    'description',  'contains behavioral data');
 
% add the Position object (that holds the SpatialSeries object)
behavior_mod.nwbdatainterface.set(...
    'Position', Position);
 
% add the processing module to the NWBFile object, and name it "behavior"
nwb.processing.set('behavior', behavior_mod);

%% Test Write
nwbExport(nwb, 'ophys_tutorial1.nwb');

