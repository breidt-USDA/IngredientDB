function PlotData = plotABD(ABDobj,trimbool)
   %NOTE: USED BY ABT.show(trimbool)
   % custom plot for ABD showing tBeta areas with shading
   %by default this algorithm uses 'equalBCends' function to trim curve
   %--------------------------------------------------------------------
   water_IS = ABDobj.titration_IS;           %water needs to match
   crv_IS = 0;                               %no additional adjustment
   ABmat = ABDobj.buffer_table{:,1:2};       %basic AB matrix
   tbd = get_tBetaData(ABmat,2,12,water_IS,crv_IS); %get data to plot
   % tbd.waterNaClpct
   % tbd.crvNaClpct                          %the data struct members
   % tbd.pHvec
   % tbd.BCcrv
   % tbd.Watercrv
   % tbd.WatertBeta
   % tbd.ModeltBeta
   BCcrv = tbd.BCcrv;                        %prepare curve for trimming
   if trimbool
      BCcrv = equalBCends(BCcrv);               %trim curve
   end
   pHmin = BCcrv(1,1);                       %get min pH for graph
   pHmax = BCcrv(end,1);                     %get max pH for graph
   %redo water curve with same ends
   pHvec = pHmin:0.05:pHmax;
   watersalt = water_IS/5.84;
   Watercrv = BetaModel_AB([0,0],watersalt,pHvec');
   %get tBeta for BCcrv
   waterData = SCBCfit_area(Watercrv,15);
   modelData = SCBCfit_area(BCcrv,15);
   tBeta = 100*(modelData.area - waterData.area);
 %plot
   figure('Name','BC Area Plot');            %name the figure
   area(BCcrv(:,1),BCcrv(:,2),'FaceColor',[0.8 0.8 0.8], ...
      'LineStyle','-');                      %BC data dark gray, black line
   hold on                                   %allow multiple plots  
   HI = area(Watercrv(:,1),Watercrv(:,2),'FaceColor',[1,1,1], ...
      'LineStyle',':');      %for alpha = 0.5 below, with dotted line
   rows = size(ABmat,1);                     %ABmatrix size
   if ABmat(1,1) ~= 0                        %If BC conc vals
      ABmat = Conc2Beta(ABmat);              %convert conc (M) to BC vals
      for i=1:rows                           %plot vertical pK lines
         line([ABmat(i,2),ABmat(i,2)],[0,ABmat(i,1)],...
            'color','black');
      end
   end
   %aphap 0.5 helps if additional data added to graph (BC original crv)
   alpha(HI, 0.5);                           %area is somewhat transparant
   xlim( [1.8 12.2]);                        %X-axis limits
   %ylim( [0,Ymax]);                         %Y-axis limits
   xlabel( 'pH');                            %X-axis label
   ylabel( 'Buffer Capacity');               %Y-axis label
   titlestr = sprintf('Buffer Area = %f', tBeta);
   th = title( titlestr);
   titlePos = [6.9250 0.0363 0];             %Default = [6.9250 0.0353 0];
   set(th, 'position', titlePos);
   if ABDobj.titration == 1 
      BCdata = ABDobj.titration_BCcurve;    %get BC data
      if trimbool
         BCdata = equalBCends(BCdata);               %trim curve
      end
      plot( BCdata(:,1),BCdata(:,2),'ok');  %plot BC data (open circles)
      PlotData.BCpts = BCdata;
   end
   hold off                                  %turn off hold 
   %save all needed data to a structure:
   PlotData.BCcrv = tbd.BCcrv;  
   PlotData.Watercrv = tbd.Watercrv;
   PlotData.tBeta = tbd.ModeltBeta;
   PlotData.WaterIS = tbd.waterNaClpct/5.84; %MW of NaCl = 5.84
   PlotData.CrvIS = tbd.crvNaClpct/5.84;
   PlotData.pHmin = pHmin;
   PlotData.pHmax = pHmax;
   PlotData.ABtable = ABDobj.buffer_table; 
end