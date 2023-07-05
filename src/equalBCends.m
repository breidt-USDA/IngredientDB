function newcurve = equalBCends(BCcurve)
%Trim BC curve so start and end BC values are similar
%Note: this is needed for accurate total BC (integration) function
%  see getTotalBC.m and/or plotABDarea.m
%-------------------------------------------------------------------------
   newcurve = BCcurve;                       %copy curve
   %trim curve to make ends same height (BC val)
   startval = newcurve(1,2);                 %start BC value of curve
   endval = newcurve(end,2);                 %end BC value of curve
   if startval < endval                      %make ends even on both sides
      if newcurve                            %make sure 1 row left!
         while (newcurve(end,2) > startval)  %start is lower than finish
            newcurve(end,:) = [];            %remove last row
         end
      end
   else                                      %finish val lower than start
      if newcurve                            %make sure 1 row left!
         while (newcurve(1,2) > endval)      %trim till even
          newcurve(1,:) = [];                %remove first row  
         end
      end
   end
   if newcurve                               %make sure curve exists
      newlen = size(newcurve,1);             %get new size
      minsize = 0.5 * BCcurve;               %minimum is half of original
      if newlen < minsize                    %check for severe reduction
         newcurve = BCcurve;                 %if so, return original
      end
   else                                      %no curve left!
      newcurve = BCcurve;                    %return original
   end
end