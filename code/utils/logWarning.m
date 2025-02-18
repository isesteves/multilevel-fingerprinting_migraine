function logWarning(warningMessage, warningLogFileName)
% 
%   logWarning(warningMessage, warningLogFileName)
%
% Description:
%   logWarning logs a warning message to a specified warning log file.
%
% Input:
%   - warningMessage: A string containing the warning message to be logged.
%   - warningLogFileName: A string specifying the path and filename of the
%                         warning log file. The log file is created or
%                         appended to if it already exists.
%
% Example:
%   logWarning('Variable x is negative.', 'logs/warning_logs/warning_log_2023-09-28.txt');
%
% See also:
%   warning
    warningLogMessage = [datestr(now,'dd-mm-yyyy HH:MM:SS FFF'), ' [WARNING] - ', warningMessage];
    fid = fopen(warningLogFileName, 'a+');
    fprintf(fid, '%s \n', warningLogMessage);
    fclose(fid);
end
