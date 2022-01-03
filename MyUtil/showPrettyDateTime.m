function prettyTimeinString = showPrettyDateTime(time, format)
%  Using datestr;
% 
%     Table 1: Standard MATLAB Date format definitions
%  
%     Number           Format                   Example
%     ===========================================================================
%        0             'dd-mmm-yyyy HH:MM:SS'   01-Mar-2000 15:45:17 
%        1             'dd-mmm-yyyy'            01-Mar-2000  
%        2             'mm/dd/yy'               03/01/00     
%        3             'mmm'                    Mar          
%        4             'm'                      M            
%        5             'mm'                     03            
%        6             'mm/dd'                  03/01        
%        7             'dd'                     01            
%        8             'ddd'                    Wed          
%        9             'd'                      W            
%       10             'yyyy'                   2000         
%       11             'yy'                     00           
%       12             'mmmyy'                  Mar00        
%       13             'HH:MM:SS'               15:45:17     
%       14             'HH:MM:SS PM'             3:45:17 PM  
%       15             'HH:MM'                  15:45        
%       16             'HH:MM PM'                3:45 PM     
%       17             'QQ-YY'                  Q1-96        
%       18             'QQ'                     Q1           
%       19             'dd/mm'                  01/03        
%       20             'dd/mm/yy'               01/03/00     
%       21             'mmm.dd,yyyy HH:MM:SS'   Mar.01,2000 15:45:17 
%       22             'mmm.dd,yyyy'            Mar.01,2000  
%       23             'mm/dd/yyyy'             03/01/2000 
%       24             'dd/mm/yyyy'             01/03/2000 
%       25             'yy/mm/dd'               00/03/01 
%       26             'yyyy/mm/dd'             2000/03/01 
%       27             'QQ-YYYY'                Q1-1996        
%       28             'mmmyyyy'                Mar2000        
%       29 (ISO 8601)  'yyyy-mm-dd'             2000-03-01
%       30 (ISO 8601)  'yyyymmddTHHMMSS'        20000301T154517 
%       31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
%  
%     Table 2: Date format symbolic identifiers (Examples are in US English)
%     
%     Symbol  Interpretation of format symbol
%     ===========================================================================
%     yyyy    full year, e.g. 1990, 2000, 2002
%     yy      partial year, e.g. 90, 00, 02
%     mmmm    full name of the month, according to the calendar locale, e.g.
%             "March", "April". 
%     mmm     first three letters of the month, according to the calendar 
%             locale, e.g. "Mar", "Apr". 
%     mm      numeric month of year, padded with leading zeros, e.g. ../03/..
%             or ../12/.. 
%     m       capitalized first letter of the month, according to the
%             calendar locale; for backwards compatibility. 
%     dddd    full name of the weekday, according to the calendar locale, e.g.
%             "Monday", "Tuesday". 
%     ddd     first three letters of the weekday, according to the calendar
%             locale, e.g. "Mon", "Tue". 
%     dd      numeric day of the month, padded with leading zeros, e.g. 
%             05/../.. or 20/../.. 
%     d       capitalized first letter of the weekday; for backwards 
%             compatibility
%     HH      hour of the day, according to the time format. In case the time
%             format AM | PM is set, HH does not pad with leading zeros. In 
%             case AM | PM is not set, display the hour of the day, padded 
%             with leading zeros. e.g 10:20 PM, which is equivalent to 22:20; 
%             9:00 AM, which is equivalent to 09:00.
%     MM      minutes of the hour, padded with leading zeros, e.g. 10:15, 
%             10:05, 10:05 AM.
%     SS      second of the minute, padded with leading zeros, e.g. 10:15:30,
%             10:05:30, 10:05:30 AM. 
%     FFF     milliseconds field, padded with leading zeros, e.g.
%             10:15:30.015.
%     PM      AM or PM
%  
    if(nargin>1)
        prettyTimeinString = datestr(time,format);
    else
        prettyTimeinString = datestr(time,'yy/mm/dd_hh:MM:SS');
    end
end