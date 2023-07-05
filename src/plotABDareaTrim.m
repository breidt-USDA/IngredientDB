function res = plotABDareaTrim(ABDvar)
%Plot ABD data with shaded total area under curve and area output
   IS = ABDvar.titration_IS;                 %get Ionic strength from ABD
   NaClpercent = 5.84*IS;                    %convert to equivalent NaCl%
   ABmat = ABDvar.buffer_table{:,1:2};       %extract AB (a_b not needed)
   pHvec = 2:0.05:12;                        %get pH vec, default range   
   BCcurve = BetaModel_AB(ABmat,0,pHvec');   %Nx2: pH and BC (M)
   Water = BetaModel_AB([0,0],NaClpercent,pHvec'); %make water curve
   %trim curve 
   temp = equalBCends(BCcurve);                           %copy curve
   %temp = BCcurve;
   SCBCdata = SCBCfit_area(temp,15);         %get SCBC fit with area calc
   %save data to output structure
   res.name = ABDvar.name;                   %ingredient name
   res.conc = ABDvar.conc;                   %ingredient concentration
   res.pH = ABDvar.pH;                       %pH
   res.adjC = ABDvar.adjC;                   %adjC
   res.TotalBC = SCBCdata.area;              %area under curve
   res.Nbuffers = size(ABmat,1);             %number of buffers
   res.minconc = min(ABmat(:,1));            %minimum conc for buffers
   res.maxconc = max(ABmat(:,1));            %maximum conc for buffers
   res.ABmat = ABmat;                        %ABmat (Nx2): Conc(M), pK
   res.area_crv = temp;                      %trimmed BC curve
   %plot
   figure('Name','BC Area Plot');            %name the figure
%    Ymax = temp(1,2);                         %initial Y value (after trim)
%    Ymax = Ymax + 0.1*Ymax;                   %add 10% to total graph height
   plot(temp(:,1),temp(:,2),'-k');           %plot the BC data (black line)
   hold on                                   %allow multiple plots        
   plot(Water(:,1),Water(:,2),'.k');         %plot water crv (black dots)
   rows = size(ABmat,1);                     %ABmatrix size
   if ABmat(1,1) ~= 0                        %If BC conc vals
      ABmat = Conc2Beta(ABmat);              %convert conc (M) to BC vals
      for i=1:rows                           %plot vertical pK lines
         line([ABmat(i,2),ABmat(i,2)],[0,ABmat(i,1)],...
            'color','black');
      end
   end
   
   area( temp(:,1),temp(:,2),'FaceColor',[0.65 0.65 0.65]); %area gray 
   area(Water(:,1),Water(:,2),'FaceColor',[0.95 0.95,0.95]); %light gray
   alpha( 0.5);                               %area is transparant
   xlim( [(BCcurve(1,1) -1) (BCcurve(end,1) + 1)]); %X-axis limits
   %ylim( [0,Ymax]);                           %Y-axis limits
   xlabel( 'pH');                             %X-axis label
   ylabel( 'Buffer Capacity');                %Y-axis label
   %titlestr = sprintf('Buffer Area = %f',res.TotalBC);
   %title( titlestr);
   try
      if ABDvar.titration == 1
         BCdata = ABDvar.titration_BCcurve;    %get BC data
         BCdata = equalBCends(BCdata);         %trim BC data for plotting
         plot( BCdata(:,1),BCdata(:,2),'ok');  %plot BC data (open circles)
         res.BCdata = BCdata;                  %add BC curve to output 
      end
   catch
      hold off
      close('BC Area Plot');                   %close plot
      %post error message that crvbool should be zero:
      error('No BC curve data available, please set crvbool = 0') 
   end
   hold off                                    %turn off hold
end