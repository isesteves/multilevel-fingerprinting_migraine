function config_tracker(config_path, pipeline_step, config)
%ARCHIVECONFIGURATION Archive and compare configuration files.
%
%   archiveConfiguration(config_path, pipeline_step, config) archives the
%   current configuration file, compares it with previous configurations,
%   and stores the new configuration if changes are detected.
%
%   INPUTS:
%   - config_path: The path to the folder containing configuration files.
%   - pipeline_step: A string indicating the pipeline step.
%   - config: The current configuration structure to be archived and compared.
%
%   Example:
%   archiveConfiguration('/path/to/config', 'preprocessing', current_config);
%
%   This function checks for previous configuration files related to the
%   specified pipeline step. If no previous files are found, the current
%   configuration is saved with a timestamp. If previous files exist, the
%   function compares the current configuration with the most recent
%   archived configuration. If changes are detected, the new configuration
%   is stored with a timestamp.

files = dir([config_path, '/archive/*', pipeline_step,'*config*']);
if isempty(files)
   disp('--- No previous configuration files found in this folder. Saving current configuration.')
   save([config_path, '/archive/', pipeline_step,'_config_',datestr(now, 'yyyymmdd-HHMMSS'),'.mat'], 'config');
else  
   [~, idx]= max([files.datenum]);
   previous_config = load([config_path, '/archive/', files(idx).name]);
   if isequal(previous_config.config, config)     
      disp('--- Current configuration matches previous one.')
   else
      save([files(idx).folder, '/', pipeline_step, '_config_',datestr(now, 'yyyymmdd-HHMMSS'),'.mat'], 'config');
      disp('--- Storing new configuration.')
   end
end