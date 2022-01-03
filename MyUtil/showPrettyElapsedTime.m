function prettyTimeinString = showPrettyElapsedTime(timeInSeconds, printTimeInSeconds)
    if(nargin==1)
        printTimeInSeconds = true;
    end
    if(timeInSeconds<=60)
        prettyTimeinString = sprintf('%2f sec', timeInSeconds );
    elseif(timeInSeconds <=3600)
        totElaspedT = timeInSeconds;
        min = floor(timeInSeconds / 60);
        sec = floor(timeInSeconds - 60*min);
        if(printTimeInSeconds)
            prettyTimeinString = sprintf('%2d min, %2d sec (%.1f sec in total)', min, sec, totElaspedT );
        else
            prettyTimeinString = sprintf('%2d min, %2d sec', min, sec);
        end
    elseif(timeInSeconds <=3600*24)
        totElaspedT = timeInSeconds;
        hr = floor(timeInSeconds / 3600);
        timeInSeconds = timeInSeconds - hr*3600;
        min = floor(timeInSeconds / 60);
        sec = floor(timeInSeconds - 60*min);
        if(printTimeInSeconds)
            prettyTimeinString = sprintf('%2d hr, %2d min, %2d sec (%.1f sec in total)', hr, min, sec, totElaspedT );
        else
            prettyTimeinString = sprintf('%2d hr, %2d min, %2d sec', hr, min, sec);
        end
    else
        totElaspedT = timeInSeconds;
        day = floor(timeInSeconds / 3600/24);
        timeInSeconds = timeInSeconds - day*3600*24;
        hr = floor(timeInSeconds / 3600);
        timeInSeconds = timeInSeconds - hr*3600;
        min = floor(timeInSeconds / 60);
        sec = floor(timeInSeconds - 60*min);
        if(printTimeInSeconds)
            prettyTimeinString = sprintf('%2d day, %2d hr, %2d min, %2d sec ( %.1f sec in total)', day, hr, min, sec, totElaspedT );
        else
            prettyTimeinString = sprintf('%2d day, %2d hr, %2d min, %2d sec', day, hr, min, sec );
        end
    end
    
end