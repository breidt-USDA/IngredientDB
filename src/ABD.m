classdef ABD   %Value class by default (assignment is deep copy)
   %BufferCapacity database class with buffer table and metadata.
   %Version 230316
   %-----------------------------------------------------------------------
   %     functions: Return ABD objects (except get_BCcurve)
   %     set_abletter - sets a specific ab letter
   %     setpH_TotBC - sets ABD pH from table, 
   %     add_buffermatrix - adds buffers from Nx2: conc (M), pH matrix
   %        uses default vals for a or b (based on pH)
   %     reset_ingredientconc - adjusts buffers, adjC based on conc change
   %        to the new concentration
   %     update_pKs - updates all pK values based on given IS
   %     get_BCcurve - gets smooth BC curve from existing BC data (200 pts)
   %     add_ABD - combines buffer tables, adjC vals
   %     set_adjC - resets adjC
   %     show - calls showABD, gives print of data and graph
   %-----------------------------------------------------------------------
   %     NOTE: add_buffertable - adds buffers from table: Conc, pK, a_b
   %     NOTE: pH automatically updated for all modifications affecting it
   %     NOTE: all internal CalcpH_ABT use IS = 0, assumes IS in pK vals...
   %-----------------------------------------------------------------------
   
   properties
      buffer_table = table;      %buffer table: Conc(M),pK,a_b('a' or 'b')
      name = 'ingredient'        %name of ingredient
      conc = 100;                %concentraiton of ingredient in solution
      adjC = 0;                  %adjC for pH calc
      pH = -1;                   %pH value from buffer table
      totalBC = -1;              %area under BC curve
      pH_cal = 0;                %curve adjustment +/- 0.5 for pH cal
      titration = 0;             %is titration data present? 0 or 1
      titration_name = 'ingredient'; %text name of ingredient
      titration_conc = 100;      %current concentration of ingredient (%)
      titration_IS = 0;          %original titration IS (molar)
      titration_files = struct;  %empty struct: format in BufferCapacity
      titration_ModelCurve = []; %SXP model curve
      titration_BCcurve = [];    %actual data curve (datapoints)
      titration_SCBCcurve = [];  %smooth data fit curve NOT model
      titration_buffertable = table; %Results table
      titration_SCBC = struct;   %SCBC results structure
      titration_SPX = struct;    %simplex results structure
      titration_SSE = 0;         %fit of BC model from titration
      titration_obspH = 0;       %observed pH from titration
      titration_volume = 0;      %titration volume
      titration_HClconc = 0;     %titration HCl conc (N)
      titration_NaOHconc = 0;    %titration NaOH conc (N)
      parameters = table;        %parameter table: format in BufferCapacity
      ver = '230316';            %version of ABD class (YYMMDD)
   end
   
   methods
   %CONSTRUCTOR takes 0 or up to 4 arguments in order:
   %  0. none, switch does nothing
   %  1. buffer_table
   %  2. buffer_table, adjC
   %  3. buffer_table, adjC, conc
   %  4. buffer_table, adjC, conc, name
   % Default values are set in properties above. With no arguments, a default
   % buffer table with 1 row (0,0,'a') is returned.  
   function obj = ABD(varargin) %buffer_table,adjC,name,conc) 
      errorbool = 0;                      %set default, no error                         
      switch nargin                       %nargin has number of args
         case 0                           %no args
            Conc = 0;                     %use default conc, pK and a_b
            pK = 0;
            a_b = {'a'};                  
            obj.buffer_table = table(Conc,pK,a_b); %create table
            obj = obj.setpH_totBC();      %calc pH and total buffering
         case 1                           %buffer table only
            buffer_table = varargin{1};   %assign data
            try
               obj.buffer_table = buffer_table; %ok, assign
               obj.buffer_table.Properties.VariableNames={'Conc',...
                  'pK','a_b'};            %make sure var names correct
               obj = obj.setpH_totBC();   %calc pH and total BC
            catch
               errorbool = 1;             %data not valid set errorbool
            end
         case 2                           %BC table and adjC
            buffer_table = varargin{1};   %table 
            adjC = varargin{2};           %adjC val
            try
               obj.buffer_table = buffer_table; %set table
               obj.buffer_table.Properties.VariableNames={'Conc',...
                  'pK','a_b'};            %set names
               obj.adjC = adjC;           %set adjC
               obj = obj.setpH_totBC();   %calc data
            catch
               errorbool = 1;             %throw error below
            end
         case 3                           %set BC table, adjC and conc
            buffer_table = varargin{1};   %BC table
            adjC = varargin{2};           %adjC
            conc = varargin{3};           %conc
            try                           
               obj.buffer_table = buffer_table; %set table
               obj.buffer_table.Properties.VariableNames={'Conc',...
                  'pK','a_b'};            %set var names
               obj.adjC = adjC;           %assing adjC
               obj.conc = conc;           %conc
               obj = obj.setpH_totBC();   %calc
            catch
               errorbool = 1;             %throw error below
            end
         case 4                           %BC table, adjC, conc, and name
            buffer_table = varargin{1};   %table
            adjC = varargin{2};           %adjC
            conc = varargin{3};           %conc
            name = varargin{4};           %name
            try
               obj.buffer_table = buffer_table; %set table
               obj.buffer_table.Properties.VariableNames={'Conc',...
                  'pK','a_b'};            %assure names correct
               obj.adjC = adjC;           %assing adjC
               obj.name = name;           %name
               obj.conc = conc;           %conc
               obj = obj.setpH_totBC();   %calc
            catch
               errorbool = 1;             %throw error below
            end 
         otherwise                        %don't assign
            errorbool = 1;                %throw error below
      end %end switch statement
      if errorbool %throw error and post message
         error('Check ABD arguments or try: ABD()');
      end
   end

   function obj = add_buffertable(obj,buffer_table)
      %returns new object with added buffer table. Does not change adjC!
      if size(buffer_table,2) == 3        %3 cols needed
         %make sure of correct variable names for combining tables
         buffer_table.Properties.VariableNames={'Conc','pK','a_b'};
         newtable = sortrows([obj.buffer_table;buffer_table],'pK'); %add, sort
         newtable(ismember(newtable.Conc,0),:)= []; %remove 0 conc rows
         %assign to new ABD with adjC. NOTE: this calculates pH, tot BC
         if isempty(newtable)
            newtable = buffer_table;
         end
         obj.buffer_table = newtable;
         obj = obj.setpH_totBC();
      else
         error('buffer_table must have: Conc, pK, and a_b columns');
      end
   end

   function obj = setpH_totBC(obj) 
      %returns ABD object with pH and totalBC set
      %set pH (note: IS = 0 as included in titration...)
      obj.pH = CalcpH_ABT(obj.buffer_table,0,obj.adjC);
      ABmat = obj.buffer_table{:,1:2};         %extract AB (a_b not needed)
      pHvec = 2:0.05:12;                       %set pH vec 2 to 12
      BCcurve = BetaModel_AB(ABmat,0,pHvec');  %Nx2: pH and BC (M)
      SCBCdata = SCBCfit_area(BCcurve,15);     %get SCBC fit with area calc
      obj.totalBC = SCBCdata.area;             %NO Correction for water
   end

   function obj = setpH_totBCtrim(~,obj,IS)
      %returns ABD object with pH (total BC set for trimmed curve)
      NaClpct = IS*5.84;                  %NaCl percent from IS
      %set pH, (note: IS = 0 as included in titration...)
      obj.pH = CalcpH_ABT(obj.buffer_table,IS,obj.adjC); 
      ABmat = obj.buffer_table{:,1:2};    %get buffer table
      pHvec = 2:0.05:12;                  %pH vec
      BCcurve = BetaModel_AB(ABmat,NaClpct,pHvec'); %Nx2: pH and BC (M)
      tempBC = equalBCends(BCcurve); %May Fail: trim curve to equal ends
      SCBCdata = SCBCfit_area(tempBC,15); %get SCBC fit with area calc
      obj.totalBC = SCBCdata.area;        %set area under curve (total BC)
   end

   function obj = add_buffermatrix(obj,buffer_matrix)
      %takes an Nx2 matrix of conc (M), pK values 
      %NOTE: this function does not update adjC!
      if size(buffer_matrix,2) == 2       %check for 2 cols
         temptable = array2table(buffer_matrix,'VariableNames',...
            {'Conc', 'pK'});              %set first two cols
         rows = height(temptable);        %get number of items
         a_b = string.empty(rows,0);      %make empty string
         for i=1:rows                     %for each item
            if temptable.pK(i) > 7        %letters assigned by pK
               a_b(i) = 'b';              %assign base letter 'b'
            else
               a_b(i) = 'a';              %assign acid letter 'a'
            end
         end 
         a_b = a_b';                      %make cols from rows
         temptable.a_b = a_b;
         %add the new table to the existing buffer table 
         obj = obj.add_buffertable(temptable);
         obj = obj.setpH_totBC();
      else
         error('buffer_matrix must have 2 cols (conc [M] and pK')
      end
   end

   function BCcurve = get_BCcurve(obj)
      %get smooth curve, assumes IS pK adjust is in buffers already
      ABmat = obj.buffer_table{:,1:2}; %extract AB (a_b not needed)
      pHvec = 2:0.05:12; %200 data points
      BCcurve = BetaModel_AB(ABmat,0,pHvec'); %generate curve
   end
   
   function obj = reset_ingredientconc(obj, newconc)
      %Update all individual buffer concentrations based on newconc
      %if old conc = 0 this does nothing. 
      if (obj.conc > 0) && (newconc > 0)  %make sure conc is 0 or positive      
         factor = newconc/obj.conc;       %generate ratio
         obj.conc = newconc;              %assign conc
         obj.buffer_table.Conc(:) = ...
            obj.buffer_table.Conc(:) * factor; %adjust all table conc vals
         obj.adjC = obj.adjC * factor;    %adjust adjC
         obj = obj.setpH_totBC();         %recalculate pH and total BC
      end
   end
   
   function newABD = add_ABD(obj,ABDvar)
      %RETURNS: NEW ABD object with the result... combined tables & adjC
      %NOTE: name = 'ingredient' and conc = 100%. 
      %NOTE: assumption is pKs in table were for the same IS
      %NOTE: combines adjC values
      %note:
            % add_ABD makes new ABD variable for conbined data!
            % add_ABD combines buffer tables, sorts by pK and eliminates 
            %  zero conc buffers
            % add_ABD sums up adjC values and uses the combined adjC to
            %  calculate the pH
            % add_ABD calculated pH with IS = 0 (all of the pKs are already
            %  ajusted)!
            % add_ABD calculated totalBC is done with NaClpercent = 0;
            % This all makes sense as new conbined ABD has buffers pK 
            % values already adjusted for IS
      newABD = ABD;                       %create new ABD object
      oldtable = obj.buffer_table;        %get existing table
      temptable = ABDvar.buffer_table;    %get table from input var
      %new table with sorted buffers and no zero conc rows
      newABD.buffer_table = sortrows([oldtable;temptable],'pK');
      newABD.buffer_table(ismember(newABD.buffer_table.Conc,0),:) = [];
      newABD.adjC = ABDvar.adjC + obj.adjC; %combine adjC values
      newABD = newABD.setpH_totBC();      %set pH and total BC (IS=0)
   end
   
   function obj = set_adjC(obj,adjC)
      %function recalculates pH and total BC
      obj.adjC = adjC;                    %reset adjC
      obj = obj.setpH_totBC();            %recalc pH and total BC
   end
   
   function obj = update_pKs(obj, IS)
      %use with caution... adjusting IS when buffers are already adjusted
      %for IS means you should only take the difference from existing IS
      %i.e. the is used for titration, and the new IS!
      tempAB = obj.buffer_table{:,1:2};
      modAB = AdjpKa_AB(tempAB,IS); %adjust all pKs
      obj.buffer_table{:,1:2} = modAB;
      obj = obj.setpH_totBC();            %recalc pH and total BC
   end

   function obj = set_abletter(obj,number,letter_string)
      %change a or b in buffer table, will update adjC so pH is same
      obj = obj.setpH_totBC();            %make sure pH is calculated
      currentpH = obj.pH;                 %save pH as target for optim
      nbuffers = size(obj.buffer_table,1); %get number of rows 
      if (number > 0) && (number <= nbuffers) %check number for range
         if (letter_string == "a") || (letter_string == "b") %check a/b
            obj.buffer_table.a_b(number) = cellstr(letter_string); %set
         end 
      end
      res = GetAdjCT(currentpH,obj.buffer_table,0); %reset adjC
      obj.adjC = res.adjC;
      obj = obj.setpH_totBC();            %set pH and tot BC
   end

   function res = show(obj,trimbool)
      %trimbool == 1 -> autotrims curve so beta matches for ends
      %else plots all data
      res = plotABD(obj,trimbool);
   end

function res = write_csv(obj, filename)
%Write a .csv spreadsheet file with selected data from the ABD struct
%input is the filename, which should have a .csv extension
%function returns a true or false indicating if the file was written
%-------------------------------------------------------------------------
%write name line: 
ABDvar = obj;
res = 1;
try
   writelines(ABDvar.name,filename);
   %prepare and write additional lines:
   concstr = strcat("Conc",",",num2str(ABDvar.conc));
   writelines(concstr,filename,'WriteMode','append');
   adjCstr = strcat("adjC",",",num2str(ABDvar.adjC));
   writelines(adjCstr,filename,'WriteMode','append');
   ISstr = strcat("Ionic Str",",",num2str(ABDvar.titration_IS));
   writelines(ISstr,filename,'WriteMode','append');
   pHstr = strcat("pH",",",num2str(ABDvar.pH));
   writelines(pHstr,filename,'WriteMode','append');
   tBetastr = strcat("tBeta",",",num2str(ABDvar.totalBC));
   writelines(tBetastr,filename,'WriteMode','append');
   %if titration data, add acid/base conc and filenames
   if ABDvar.titration
      NaOHline = strcat("NaOH Conc",",", ...
         num2str(ABDvar.titration_NaOHconc),",",...
         ABDvar.titration_files.NaOH_titration_filename);
      writelines(NaOHline,filename,"WriteMode","append");
      HClline = strcat("HCl Conc",",", ...
         num2str(ABDvar.titration_HClconc),",", ...
         ABDvar.titration_files.HCl_titration_filename);
      writelines(HClline,filename,"WriteMode","append");
   else
      writelines("No titraton data",filename,...
         "WriteMode","append");
      writelines("BC model generated from the buffer table",filename,...
         "WriteMode","append");
   end
   %write buffer table
   writelines("",filename,'WriteMode','append'); %blank line
   headerstr = strcat("Conc",",","pK",",","a/b");
   writelines(headerstr,filename,'WriteMode','append');
   nbuff = size(ABDvar.buffer_table,1);
   outstr = strings;
   for i=1:nbuff
      for j=1:2
         outstr = strcat(outstr,num2str(ABDvar.buffer_table{i,j}),",");
      end
      outstr = strcat(outstr,ABDvar.buffer_table{i,3});
      writelines(outstr,filename,'WriteMode','append');
      outstr = strings;
   end
catch
   res = 0;
end %end of try block
end %end of function

   end %end of methods
end

