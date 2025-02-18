function logCustom(customMessage, customLogFileName)
%
%   logCustom(customMessage, customLogFileName)
%
% Description:
%   logCustom logs a custom message to a specified custom log file.
%
% Input:
%   - customMessage: A string containing the custom message to be logged.
%   - customLogFileName: A string specifying the path and filename of the
%                       custom log file. The log file is created or
%                       appended to if it already exists.
%
% Example:
%   logCustom('Custom log entry: Data processing completed.', 'logs/custom_logs/custom_log.txt');
%
% See also:
%   disp

    customLogMessage = [datestr(now,'dd-mm-yyyy HH:MM:SS FFF'), ' - ', customMessage];
    fid = fopen(customLogFileName, 'a');
    fprintf(fid, '%s\n', customLogMessage);
    fclose(fid);
end
