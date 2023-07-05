function res = plotABDarea_hBC(hBC, ABDvar, title_string, trimbool)
%Plot ABD data with shaded total area under curve and area output
   IS = ABDvar.titration_IS;                 %get Ionic strength from ABD
   NaClpercent = 5.84*IS;                    %convert to equivalent NaCl%
   ABmat = ABDvar.buffer_table{:,1:2};       %extract AB (a_b not needed)
   pHstart = 2;                              %default start of pH crv
   pHend = 12;                               %default end of pH crv
   pHvec = pHstart:0.05:pHend;               %get pH vec, default range   
   BCcurve = BetaModel_AB(ABmat,0,pHvec');   %Nx2: pH and BC (M)
   Water = BetaModel_AB([0,0],NaClpercent,pHvec'); %make water curve
   %trim curve 
   %copy curve and trim if trimbool == 1
   temp = BCcurve;
   if trimbool
      temp = equalBCends(BCcurve);           %make ends equal
      pHstart = temp(1,1);
      pHend = temp(end,1);
   end
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
   res.pHstart = pHstart;                    %pH (trimmed or not) pH start
   res.pHend = pHend;                        %pH (trimmed or not) pH end
   %plot
%   figure('Name','BC Area Plot');          %name the figure
%    Ymax = temp(1,2);                      %initial Y value (after trim)
%    Ymax = Ymax + 0.1*Ymax;                %add 10% to total graph height
   plot(hBC,temp(:,1),temp(:,2),'-k');       %plot the BC data (black line)
   hold(hBC,"on");                           %allow multiple plots        
   plot(hBC,Water(:,1),Water(:,2),'.k');     %plot water crv (black dots)
   rows = size(ABmat,1);                     %ABmatrix size
   if ABmat(1,1) ~= 0                        %If BC conc vals
      ABmat = Conc2Beta(ABmat);              %convert conc (M) to BC vals
      for i=1:rows                           %plot vertical pK lines
         line(hBC,[ABmat(i,2),ABmat(i,2)],[0,ABmat(i,1)],...
            'color','black');
      end
   end
   
   area(hBC, temp(:,1),temp(:,2),'FaceColor',[0.65 0.65 0.65],...
      'LineStyle','none'); %gray
   area(hBC, Water(:,1),Water(:,2),'FaceColor',[0.95 0.95,0.95], ...
      'LineStyle','none'); %light gray
   alpha(hBC, 0.5);                               %area is transparant
   xlim(hBC, [(BCcurve(1,1) -1) (BCcurve(end,1) + 1)]); %X-axis limits
   xlabel(hBC, 'pH');                             %X-axis label
   ylabel(hBC, 'Buffer Capacity');                %Y-axis label
   titlestr = sprintf(title_string);
   title(hBC, titlestr);
   try
      if ABDvar.titration == 1
         BCdata = ABDvar.titration_BCcurve;    %get BC data
         plot(hBC, BCdata(:,1),BCdata(:,2),'ok');  %plot BC data (open circles)
         res.BCdata = BCdata;                  %add BC curve to output 
      end
   catch
      hold(hBC,"off");
      close(hBC,'BC Area Plot');                   %close plot
   end
   hold(hBC,"off");                                    %turn off hold
end