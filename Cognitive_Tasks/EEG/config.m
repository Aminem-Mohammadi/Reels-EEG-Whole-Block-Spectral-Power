% config.m
%
% Set DATA_DIR and OUT_DIR once, before running any script in this folder.
%
% DATA_DIR must point to the folder of preprocessed EEGLAB .set/.fdt
% files for the N-back/SART/Full-session recordings (one pair of files
% per subject x task x block), including channel_locations_template.mat.
%
% OUT_DIR is the analysis-output root. Every script below writes its
% output under OUT_DIR (in a <Full|Nback|SART>/<NN_name>/<RUN_MODE>/
% structure, created automatically); nothing is ever written under
% DATA_DIR.

DATA_DIR = '<path to the preprocessed cognitive-task data>';   % set once before running
OUT_DIR  = fullfile(pwd, 'results');   % set once before running
