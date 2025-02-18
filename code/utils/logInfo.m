function logInfo(infoMessage, infoLogFileName)
%
%   logInfo(infoMessage, infoLogFileName)
%
% Description:
%   logInfo logs an information message to a specified info log file.
%
% Input:
%   - infoMessage: A string containing the information message to be logged.
%   - infoLogFileName: A string specifying the path and filename of the
%                     info log file. The log file is created or appended to
%                     if it already exists.
%
% Example:
%   logInfo('Data loaded and processed successfully.', 'logs/info_logs/info_log_2023-09-28.txt');
%
% See also:
%   disp
    infoLogMessage = [datestr(now,'dd-mm-yyyy HH:MM:SS FFF'), ' [INFO] - ', infoMessage];
    fid = fopen(infoLogFileName, 'a+');
    fprintf(fid, '%s \n', infoLogMessage);
    fclose(fid);
end