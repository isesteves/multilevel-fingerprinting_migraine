function logError(errorMessage, errorLogFileName)
%
%   logError(errorMessage, errorLogFileName)
%
% Description:
%   logError logs an error message to a specified error log file.
%
% Input:
%   - errorMessage: A string containing the error message to be logged.
%   - errorLogFileName: A string specifying the path and filename of the
%                       error log file. The log file is created or
%                       appended to if it already exists.
%
% Example:
%   logError('An error occurred: File not found.', 'logs/error_logs/error_log_2023-09-28.txt');
%
% See also:
%   try, catch, getReport

    errorLogMessage = [datestr(now,'dd-mm-yyyy HH:MM:SS FFF'), ' [ERROR] - ', errorMessage];
    fid = fopen(errorLogFileName, 'a+');
    fprintf(fid, '%s \n', errorLogMessage);
    fclose(fid);
end

