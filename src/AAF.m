classdef AAF < handle
   %Acid Acidified Food data class
   %Class to hold and process ABD data from ingredient lists and an acid 
   % Acid ABD vars are all start with 0 conc and 0 salt-of-acid. 
   % Acid ABD var pK values are adjusted for IS (based on ingABDvar salt)
   %   when combined with ingABDvar for returning allABD var. 
   % Acid ABD vars have adjC based on salt-of-acid
   %Public Properties:
   % Access is available to all ABD lists, combined ABDs and tables
   %  NOTE that the ABDs, ABD lists, and tables are kept up to date with
   %  any public function access change (see private function set_allABDs)
   %Public Functions:
   % NOTE: all public functions call set_allABD to update public variables
   % 1. Constructor (no args)
   %     formula = water with no salt by default. Initial post through
   %     set_allABD
   % 2. set_acidCS(Name_acid2set, acidconc%, saltofacid = 0,1,2,3)
   %     based on name of acid, the concentration (in percent) and salt of  
   %     acid is set. Salt of acid value can be 1, 2, or 3 for triprotic, 1 
   %     or 2 for diprotic, and 1 for monoprotic. Updates posted through a
   %     call to set_allABD
   % 3. add_ingABD(ingredientABD)
   %     add an ingredient to the formula, all titration_IS must be same.
   %     Updates posted through call to set_allABD
   % 4. remove_ingABD(index)
   %     remove an ingredient based on index in the ingredienttable.
   %     Updates posted through call to set_allABD
   % 5. reset_ingABDconc(name, newconc) NOTE: for ingdients NOT acids...
   %     changes the concentration of ingredient by name (conc always > 0).
   %     Updates posted through call to set_allABD
   % 6. reset_trimbool(obj,trimbool)
   %     sets the trimbool variable and resets all tBeta values through
   %     call to set_allABD
   % 7. reset_IS(IS)
   %     calls set_allABD with the new IS value. 
   % 8. clear_ingABDvar()
   %     removes all ingredients from the ingABDvar list, resetting the
   %     list to an empty list. Updates posted through call to set_allABD
   %Private Functions:
   % Most private functions support posting to tables and are called by the
   % public functions
   % NOTE: set_allABD(IS) is the call that updates all tables, lists and 
   %  public variables.
   %----------------------------------------------------------------------
   properties
      allABDobj = ABD;        %summary ABD
      ingABDobj = ABD;        %ingredient ABD
      acidABDobj = ABD;       %acid ABD
      ingABDvar = ABD.empty;  %empty list of ingredient ADB class variables
      acidABDvar = ABD.empty; %empty list for acid ABD class variables
      acidtable = table;      %table with all acids 
      ingredienttable = table; %table with all ingredients
      bufferinfotable = table; %buffer table for mixture with named buffers
      tBetatable = table;     %table for total buffering
      ingnamelist = {};       %list of ingredient names
      acidnamelist = {};      %list of acid names (even if conc = 0)
      acidconcvec = [];       %acid conc vector   
      ingconcvec = [];        %ingredient conc vector
      ingredient_conc = 0;    %conc of ingedients, percent: NOT allABD.conc
      acid_conc = 0;          %conc of acids in percent
      name = '';              %set by default to 'formulation'
      conc = 100;             %default conc for allABD
      adjC = 0;               %total adjC
      IS = 0;                 %IS of ingredients (not including acid salts)
      pH = 7;                 %pH of mixture
      tBeta = 0;              %calculated tBeta
      acidpercent = 0;        %percent of formulation tBeta for acids   
      ingpercent = 0;         %percent of formulation tBeta for ingredients
      trimbool = false;       %trim bool for tBeta calc and tBeta data
   end

   properties (Access = private)
      %private to the program
      acidMW_adjC = [];       %acid mw list
      acid_index = 7;         %number of acids
      ver = '1.000'           %current version of AAF
   end

   methods
      function obj = AAF()
      %--------------------------------------------------------------------   
      %CONSTRUCTOR, no args
      % only contains acids but all at 0 percent conc
      %NOTE: Conc must be update manually, since ABD.reset_ingredientconc
      %does nothing if conc is initially zero!
      %NOTE: polyprotic acids can be textbook values or measured...
      %--------------------------------------------------------------------
         %build acid ABD list: set conc, pK values, and a_b (all acids)
         %LACTIC:
         abmat = table(0,3.86,{'a'}); %each acid at zero conc
         obj.acidABDvar(1,1) = ABD(abmat,0,0,'Lactic');
         %ACETIC:
         abmat = table(0,4.76,{'a'});
         obj.acidABDvar(2,1) = ABD(abmat,0,0,'Acetic');
         %CITRIC:
         %Cit05pct, AdjC 0.0012, TI_00042 and 43, rpt 7 (CS 10282022)
         abmat = table([0;0;0],[3.0190;4.4745;5.7336],{'a';'a';'a'}); 
         obj.acidABDvar(3,1) = ABD(abmat,0,0,'Citric');
         %textbook
            % abmat = table([0;0;0],[3.128;4.761;6.396 ...
            %    ],{'a';'a';'a'}); %6.396
            % obj.acidABDvar(3,1) = ABD(abmat,0,0,'Citric');
         %MALIC:
         %Malic045pct, AdjC 8.125e-4, Ti_00046 and 47, (CS 10282022)
         abmat = table([0;0],[3.3246;4.7331],{'a';'a'});
         obj.acidABDvar(4,1) = ABD(abmat,0,0,'Malic');
         %textbook:
            % abmat = table([0;0],[3.40;5.20],{'a';'a'});
            % obj.acidABDvar(4,1) = ABD(abmat,0,0,'Malic');
         %ASCORBIC:
         %adjust for pK measured
         % abmat = table(0,4.291,{'a'});
         %textbook: NOTE: 11.6 is for acidic hydroxyls, not base, leave
         %out?
         abmat = table(0,4.17,{'a'});
         %abmat = table([0;0],[4.17;11.6],{'a';'a'});
         obj.acidABDvar(5,1) = ABD(abmat,0,0,'Ascorbic');
         %GLUCNIC:
         abmat = table(0,3.86,{'a'});
         obj.acidABDvar(6,1) = ABD(abmat,0,0,'Gluconic');
         %PHOSPHORIC:
         abmat = table([0;0;0],[2.14;7.20;12.37],{'a';'a';'a'});
         obj.acidABDvar(7,1) = ABD(abmat,0,0,'Phosphoric');
         %set MW list
         obj.acidMW_adjC = [90.08, 0; 60.052, 0; 192.124, 0; ...
            134.0874, 0; 176.12, 0; 196.16, 0; 97.994, 0];
         %finally, prepare the ABD formulation with all acid conc = 0 and
         % no other ingredients (argument is IS = 0)
         obj.set_allABD(0); 
      end

      function namelist = get_ingnamelist(obj)
      %--------------------------------------------------------------------
      % return a col vec of ingredient names (cell array)
      %--------------------------------------------------------------------
         namelist = {};                         %default return val
         if ~isempty(obj.ingABDvar)             %check for ingredients
            namelist = {obj.ingABDvar(:).name}'; %col vec of names
         end
      end

      function obj = set_acidCS(obj,acidname, ...
            acidconc_percent,saltofacid)
      %--------------------------------------------------------------------
      % Set acid conc in percent and adjC for salt of an acid in molar
      %    units. The acidname variable indicates which acid to set
      % NOTE: updates acid conc percent, buffer table conc(s), pH, totBC
      % NOTE: automatically updates allABD 
      %--------------------------------------------------------------------
         if saltofacid == 0                  %for setting adjC value
            mult = 0;                        %default adjC will be 0
         else
            mult = -1;                       %negative val for cation                     
         end
         switch acidname                     %index and mult by acidname
            case 'Lactic'
               index = 1;                    %from order of acids in list
            case 'Acetic'
               index = 2;
            case 'Citric'                    %triprotic, get salt number
               index = 3;
               if ismember(saltofacid,[1,2,3]) %must be in list
                  mult = -1*saltofacid;      %negative of selected val
               end
            case 'Malic'
               index = 4;
               if ismember(saltofacid,[1,2]) %diprotic 
                  mult = -1*saltofacid;      %negative of selected val
               end
            case 'Ascorbic'
               index = 5;
            case 'Gluconic'
               index = 6;
            case 'Phosphoric'
               index = 7;
               if ismember(saltofacid,[1,2,3]) %triprotic
                  mult = -1*saltofacid;      %negative of selected val
               end
         end
         %NOTE: index to acidMW_adjC list and acidABDvar list is same:
         obj.acidABDvar(index).conc = acidconc_percent; %new percent conc
         %so not zero,because reset does not work if zero...
         obj.acidABDvar(index).reset_ingredientconc(acidconc_percent);
         MW = obj.acidMW_adjC(index,1);      %get MW for acid
         buffer_conc = acidconc_percent/(MW * 0.1); %precent to Molar
         obj.acidABDvar(index).buffer_table.Conc(:) = buffer_conc; %assign
         obj.acidABDvar(index) = ... %NOTE: setting adjC updates pH, BC
            obj.acidABDvar(index).set_adjC(mult*buffer_conc); %update adjC
         obj.set_allABD(0);      %reset allABD

      end

      function obj = add_ingABD(obj, ingABD)
      %--------------------------------------------------------------------
      % Add new ingredient ABD to formulation
      %--------------------------------------------------------------------
         if isempty(obj.ingABDvar)           %if empty just assign
            obj.ingABDvar = ingABD;          %assign new ABD
         else 
            leng = length(obj.ingABDvar);    %get next index value
            obj.ingABDvar(leng+1,1) = ingABD; %assign new ABD
         end
         obj.set_allABD(0);                   %reset allABD formulation
      end

      function obj = reset_ingABDconc(obj, name, newconc)
      %--------------------------------------------------------------------
      % reset the concentration (percent) value of an ingredient ABD
      % NOTE: newconc > 0, if conc = 0, use delete instead
      %--------------------------------------------------------------------
         len = length(obj.ingABDvar);           %get curent ing list length
         if (len > 0) && (newconc > 0)          %only do if conc is not 0
            for i=1:len                         %ingredients have conc > 0
               found = strcmp(obj.ingABDvar(i).name, name);
               if found
                  obj.ingABDvar(i) = ...        %update ingABD
                     obj.ingABDvar(i).reset_ingredientconc(newconc);
                  obj.set_allABD(0);            %reset allABD formulation
               end
            end
            
         end
      end

      function obj = remove_ingABD(obj, name)
      %--------------------------------------------------------------------
      % remove an ingredient ABD by name
      % Assume input is text string
      %--------------------------------------------------------------------
         len = length(obj.ingABDvar);        %get curent ing list length
         if len > 0                          %for each item check name
            for i=1:len
               found = strcmp(obj.ingABDvar(i).name, name);
               if found                      %name found
                  obj.ingABDvar(i) = [];     %remove ingredient ABD data 
                  obj.set_allABD(0);         %post data
                  return
               end
            end 
         end
      end

      function obj = reset_IS(obj, IS)
      %--------------------------------------------------------------------
      % reset ionic strength
      %--------------------------------------------------------------------
         obj.set_allABD(IS); %sets all ABDs using IS value
      end

      function obj = clear_ingABDvar(obj)
      %--------------------------------------------------------------------
      % remove all ingredient ABDs (not acids)
      %--------------------------------------------------------------------
         while ~isempty(obj.ingABDvar)       %for each itme
            obj.ingABDvar(1) = [];           %pop from list
            obj.set_allABD(0);               %reset all ABDs with zero IS
         end
      end

      function obj = reset_trimbool(obj,trimbool)
         obj.trimbool = trimbool;
         obj.set_allABD(0);
      end

   end %<<<<<<<<<<<<<<<<< End of public methods >>>>>>>>>>>>>>>>>>>>>>>>>>>

   methods (Access = protected)

      function obj = set_allABD(obj,IS)
      %--------------------------------------------------------------------
      %Returns one ABD with all buffers, IS from ingredients,
      %  with acid buffers updated for IS, and the pH and total buffering 
      %NOTE: get formulation from ingABDvar and acidABDvar based on conc 
      %  of each ingredient and acid
      %NOTE: Ingredients are added together only if IS is the same for all
      %NOTE: All added ingredients must have positive percent conc
      %NOTE: The formulation is built from a new ABD class variable
      %  add_ABD combines buffer tables, sorts by pK and eliminates 
      %     zero conc buffers
      %  add_ABD sets the ABD formulation conc variable to 100%
      %  add_ABD sums up adjC values and uses the combined adjC to
      %     calculate the pH
      %  add_ABD calculated pH with IS = 0 (assume the pKs are already
      %     ajusted) and the calculated totalBC is also done with IS=0
      %  This all makes sense as new conbined ABD has buffers pK 
      %     values already adjusted for IS
      %--------------------------------------------------------------------
         allABDtemp = ABD;                   %default ABDs
         ingABDobjtemp = ABD;
         acidABDobjtemp = ABD;
         %to build buffer table for ingredients, check IS is same for all
         samebool = 1;                       %assume all IS the same
         titration_IS = IS;                  %defalut titration IS is 0
         if ~isempty(obj.ingABDvar)          %make sure all IS same
            diffmat = obj.ingABDvar(:).titration_IS; %get titration IS list
            %all titrations agree on IS value or no ing added to allABD
            samebool = find(all(~diff(diffmat))); %not same, samebool = 0;
         end
         ing_index = length(obj.ingABDvar);  %get number of ingredients
         if samebool && ing_index > 0        %IS is same for all in list
            for i = 1:ing_index              %for each ingredient 
               if obj.ingABDvar(i).conc > 0  %ingredients must have conc
                  %the add_ABD function (from ABD class) posts pH, tBeta 
                  allABDtemp = allABDtemp.add_ABD(obj.ingABDvar(i)); %add 
               end
            end
            %here, ingABD is the same as the allABD since acids not added
            ingABDobjtemp = allABDtemp;      %save combined ingredient ABD 
            %now set temp titration_IS, use first index:
            titration_IS = obj.ingABDvar(1).titration_IS;
         end
         for j = 1:obj.acid_index            %scroll through list of acids
            if obj.acidABDvar(j).conc > 0    %check conc not zero
               acidABDtemp = obj.acidABDvar(j); %get ABD if has conc
                  %reset pK values for titration_IS from ingredients:
               acidABDtemp = acidABDtemp.update_pKs(titration_IS);
                  %load acidABD obj and also finish building allABD 
               acidABDobjtemp = obj.acidABDobj.add_ABD(acidABDtemp);
               allABDtemp = allABDtemp.add_ABD(acidABDtemp); %add new acid
            end
         end
         %after adding all acids or ingredients: buffer table is done,
         % conc is set to 100%, adjC values are summed, totalBC is set, pH
         % is set. All that remains is to set name and titration_IS
         %titration value is left zero
         allABDtemp.name = 'formulation';    %default allABD name
         ingABDobjtemp.name = 'ingredients';
         acidABDobjtemp.name = 'acids';
         allABDtemp.titration_IS = titration_IS; %set after all adds
         ingABDobjtemp.titration_IS = titration_IS;
         acidABDobjtemp.titration_IS = titration_IS;
         %now set all object properties:
         obj.allABDobj = allABDtemp;         %set allABDs
         obj.ingABDobj = ingABDobjtemp;       
         obj.acidABDobj = acidABDobjtemp;
         obj.ingredient_conc = obj.get_ingconc(); %total ingredient conc
         obj.acid_conc = obj.get_acidconc(); %get total acid conc
         obj.name = allABDtemp.name;         %set the name  
         obj.adjC = allABDtemp.adjC;         %set adjC val
         obj.pH = allABDtemp.pH;             %set pH
         obj.IS = allABDtemp.titration_IS;   %set IS
         obj.reset_tBetas(obj.trimbool);     %post tBeta val & table
         obj.bufferinfotable = obj.get_bufferinfotable(); %post buffers
         obj.acidtable = obj.get_acidtable(); %post acid table
         obj.ingredienttable = obj.get_ingredienttable(); %post ingredients
         obj.ingnamelist = obj.get_ingnamelist(); %post ingredient names
         obj.acidnamelist = obj.get_acidnamelist(); %post acid names
         obj.ingconcvec = obj.get_ingconcvec(); %post ingredient concs
         obj.acidconcvec = obj.get_acidconcvec(); %post acid concs
      end

      function obj = reset_tBetas(obj,trimbool)
         %reset tBeta to account for trimming of ingredients
         obj.trimbool = trimbool;            %post the trim bool status
         obj.tBeta = obj.get_tBeta(obj.allABDobj,trimbool);  %post tBeta
         obj.tBetatable = obj.get_tBetatable(trimbool); %post tBeta table
         obj.set_percentages(trimbool);      %to set tBeta percentages
      end

      function tBeta = get_tBeta(obj,ABDobj, trimbool)
      %--------------------------------------------------------------------
      % returns tBeta from allABD with or without trim
      % NOTE: if trimbool set to 1, recalculate the total BC area
      %--------------------------------------------------------------------
         totalBC = ABDobj.totalBC;              %total BC from allABD
         if totalBC > 0                         %if acids or ingreds
            if trimbool == 1 && ~isempty(obj.ingABDvar) %trim and low acids
               ABmat = ABDobj.buffer_table{:,1:2}; %get AB matrix
               pHvec = 2:0.05:12;               %get default pH vec   
               BCcurve = BetaModel_AB(ABmat,0,pHvec'); %Nx2: pH and BC (M)
               temp = equalBCends(BCcurve);     %make ends equal
               SCBCdata = SCBCfit_area(temp,15); %SCBC fit with area
               totalBC = SCBCdata.area;         %assign new totalBC
            end
            tBeta = round(100*(totalBC-0.02),3); %convert totalBC to tBeta
         end
      end

      function tBetatable = set_percentages(obj,trimbool)
      %--------------------------------------------------------------------
      % set the acid and ingredient percent boxes (using trim)
      %--------------------------------------------------------------------
         %get all ingredient names
         ingtBetatable = obj.get_nametBetatable(obj.ingABDvar,trimbool);
         %get all acid names
         acidtBetatable = obj.get_nametBetatable(obj.acidABDvar,0);       
         tBetatable = [ingtBetatable;acidtBetatable]; %make new tBeta table
         %separate acid and ingredient total tBeta values
         acidtotaltBeta = obj.get_tBetatotal(acidtBetatable);
         ingtotaltBeta = obj.get_tBetatotal(ingtBetatable);
         totaltBeta = acidtotaltBeta + ingtotaltBeta; %sum for percent
         %default percent settings
         acidpercentval = 0;
         ingpercentval = 0;
         if totaltBeta > 0   %only adjust percentages if acids, ing present 
            acidpercentval = 100*acidtotaltBeta/totaltBeta;
            ingpercentval = 100*ingtotaltBeta/totaltBeta;
         end %post the calculated values
         obj.acidpercent = acidpercentval;
         obj.ingpercent =ingpercentval;
      end

      function tBetatotal = get_tBetatotal(~,nametab)
      %--------------------------------------------------------------------
      % helper for set_percentages function, get the total tBeta from
      % ingredients or acids based on the name table with both (if present)
      %--------------------------------------------------------------------
         tBetatotal = 0;                     %default return
         if ~isempty(nametab)                %check for items in table
            tBetatotal = sum(nametab.tBeta); %sum all tBeta field values
         end
      end


      function tBetatable = get_tBetatable(obj,trimbool)
      %--------------------------------------------------------------------
      % return a table with Name, tBeta for ingredients and acids
      %--------------------------------------------------------------------
         ingtBetatable = obj.get_nametBetatable(obj.ingABDvar,trimbool);
         %only gets acids with some conc val set
         acidtBetatable = obj.get_nametBetatable(obj.acidABDvar,0);       
         tBetatable = [ingtBetatable;acidtBetatable];      
         totaltBeta = sum(tBetatable.tBeta); %total tBeta
         tBetatable.percent = 100*tBetatable.tBeta/totaltBeta; %percentages
      end

      function nametBeta = get_nametBetatable(obj,ABDvar,trimbool)
      %--------------------------------------------------------------------
      % get name table from either acid or ingredient ABDvar
      % table has variables: Name, percent, and tBeta
      % NOTE: does not set the percent col (remains zeros)
      %--------------------------------------------------------------------
         nametBeta = table;                  %default empty table
         len = length(ABDvar);               %number of items in list    
         if len > 0   %even if only acids there are items with conc
            nametBeta = table('Size',[len,3],'VariableTypes', ...
               {'string','double','double'},'VariableNames', ...
               {'Name','percent','tBeta'});
            Name = {ABDvar(:).name}';        %name col vec
            tBetaval = zeros(len,1);         %blank tBeta col vector
            tBetapercent = zeros(len,1);     %blank percent col
            for i=1:len
               %now get each tBeta based on trimbool
               ABDobj = ABDvar(i).setpH_totBC(); %totalBC is correct
               tBetaval(i) = obj.get_tBeta(ABDobj,trimbool); %calc tBeta
            end
            nametBeta.Name = Name;           %set name list
            nametBeta.tBeta = tBetaval;      %set tBeta list
            nametBeta.percent = tBetapercent;  %placeholder (zeros)
            %remove all zero tBeta entries (i.e. from acid list)
            nametBeta(nametBeta.tBeta == 0,:) = []; %sort and remove empty
         end
      end 

      function bufferinfo = get_bufferinfotable(obj)
      %--------------------------------------------------------------------
      % retun buffer info table with all buffers in the formulation, 
      %  including: N, Conc, pK, a_b, and Name of each buffer
      %--------------------------------------------------------------------
         ing_tabledata = obj.get_bufferdata(obj.ingABDvar); %ing table
         tempABDacid = obj.acidABDvar;       %copy acidABDvar
         len = length(tempABDacid);          %length, even if conc = 0    
         for i=1:len                         %update pKs for IS
            tempABDacid(i) = tempABDacid(i).update_pKs(obj.IS);
         end
         acid_tabledata = obj.get_bufferdata(tempABDacid); %get acid table
         bufferinfo = [ing_tabledata;acid_tabledata]; % merge tables
         bufferinfo(ismember(bufferinfo.Conc,0),:) = []; %no 0 needed
         bufferinfo = sortrows(bufferinfo,'pK'); %sort by pK
         finallen = size(bufferinfo,1);      %number of buffers total
         % N = 1:finallen;                   %make vector 
         % N = N';                           %change to col vector
         bufferinfo.N = (1:finallen)';       %add to table
         bufferinfo = bufferinfo(:, ... 
            {'N','Conc','pK','a_b','Name'}); %move N to first col
      end

      function ingredienttable = get_ingredienttable(obj)
      %--------------------------------------------------------------------
      % return a table with N, Name, Conc for each ingredient
      %--------------------------------------------------------------------
         ingredienttable = table;             %for default of empty ingABD
         if ~isempty(obj.ingABDvar)  
            leng = length(obj.ingABDvar);     %get number of ingredients
            %preallocate table
            ingredienttable = table('Size',[leng,3],'VariableTypes', ...
               {'double','string','double'},'VariableNames', ...
               {'N','Name','Conc'});
            ingredienttable.N = (1:leng)';     %numbered ingredients
            ingredienttable.Name = {obj.ingABDvar(:).name}'; %names
            ingredienttable.Conc = [obj.ingABDvar(:).conc]'; %conc vals
            %sort by concentration
            ingredienttable = sortrows(ingredienttable,3,'descend');
            ingredienttable.N = (1:leng)';     %renumber ingredients
         end
      end

      function acidinfo = get_acidtable(obj)
      %--------------------------------------------------------------------
      % return a table with N, Name, Conc for each ingredient
      %--------------------------------------------------------------------
         len = obj.acid_index;               %number of acids
         %preallocate table   
         acidinfo = table('Size',[len,6],'VariableTypes', ...
            {'double','string','logical','double','double','double'}, ...
            'VariableNames', {'N','Name','Salt','Conc','mM','MW'});
         acidinfo.N = (1:len)';              %numbered ingredients
         acidinfo.Name = {obj.acidABDvar(:).name}'; %names
         acidinfo.Salt = logical([obj.acidABDvar(:).adjC]');
         acidinfo.Conc = [obj.acidABDvar(:).conc]'; %conc vals
         acidinfo.mM = 1000*(obj.get_acidmMvec());  %colvec
         acidinfo.MW = obj.acidMW_adjC(:,1); %MW 
      end

      function conc = get_ingconc(obj)
      %--------------------------------------------------------------------
      % return a sum of the concentrations for all ingredients
      %--------------------------------------------------------------------
         conc = 0;                           %default return
         if ~isempty(obj.ingABDvar)          %check for ingABD
            concvec = [obj.ingABDvar(:).conc]; %get vector of conc
            conc = sum(concvec);             %sum vector
         end
      end

      function conc = get_acidconc(obj)
      %--------------------------------------------------------------------
      % return a sum of the concentrations for all acids
      % NOTE: acidABDvar exists by default
      %--------------------------------------------------------------------
         concvec = [obj.acidABDvar(:).conc]; %acids may have zero conc
         conc = sum(concvec);                %sum conc
      end

      function acidmMvec = get_acidmMvec(obj)
      %--------------------------------------------------------------------
      % return a vector of acid mM values from existing percentages
      % NOTE: acidABDvar exists by default, but conc may be zero
      %--------------------------------------------------------------------
         len = obj.acid_index;               %get number of acids
         acidmMvec = zeros(len,1);           %preallocate vector
         for i=1:len
            acidmMvec(i) = obj.acidABDvar(i).conc/ ... %convert percent 
               (obj.acidMW_adjC(i) * 0.1);             % to mM
         end
      end

      function conc = get_adjC(obj)
      %--------------------------------------------------------------------
      % return the sum total for all ingredient adjC values for ingredients
      % as well as acids. 
      % NOTE: acidABDvar exists by default
      %--------------------------------------------------------------------
         conc = 0;                           %default return
         if ~isempty(obj.ingABDvar)          %check for ingrdients
            conc = sum(obj.ingABDvar(:).conc); %if so, add conc vals
         end
         conc = conc + sum(obj.acidABDvar(:).conc); %add acid conc vals
      end

      function nbuffers = get_nbuffers(~,ABDvar)
      %--------------------------------------------------------------------
      % return the sum of all buffers for either ingABDvar or acidABDvar
      %--------------------------------------------------------------------
         nindex = length(ABDvar);            %get number of ABDs
         nbuffers = 0;                       %default return val
         for i=1:nindex                      %for each ABD
            tempABD = ABDvar(i);             %get the ABD from the list
            btab = tempABD.buffer_table;     %extract buffer table
            len = size(btab,1);              %determine rows in table
            nbuffers = nbuffers + len;       %sum row val
         end
      end
      
      function bufferdata = get_bufferdata(obj,ABDvar)
      %--------------------------------------------------------------------
      % Given an ingABDvar or acidABDvar, return buffer table with names
      % buffer data table has four variables:
      %  Conc, pK, a_b, Name
      %--------------------------------------------------------------------
         nindex = length(ABDvar);            %number of ABD in vec
         %setup blank data table
         nbuffers = obj.get_nbuffers(ABDvar);
         bufferdata = table('Size',[nbuffers,4],'VariableTypes',... 
            {'double','double','cellstr','cellstr'}, 'VariableNames', ... 
            {'Conc', 'pK', 'a_b', 'Name'});
         dtabindex = 1;                      %first empty row
         for i=1:nindex                      %for each ABD in ABDvar
            temptable = ABDvar(i).buffer_table; %get the table (
            name_temp = ABDvar(i).name;      %get the name 
            len = size(temptable,1);         %number of buffers
            Name = [repelem({name_temp},len)]'; %col vec with cell str name 
            temptable.Name = Name;           %assign names
            %add temporary table to datatable at appropriate index
            start = dtabindex;               %first empty row
            finish = dtabindex + len - 1;    %n more, including first
            bufferdata(start:finish,:) = temptable(1:len,:); %add table data
            dtabindex = dtabindex + len;     %set for next time...
         end   
      end

      function namelist = get_acidnamelist(obj)
      %--------------------------------------------------------------------
      % return a col vec of acid names (cell array)
      %--------------------------------------------------------------------
         namelist = {obj.acidABDvar(:).name}'; %col vec of names 
      end

      function ingconcvec = get_ingconcvec(obj)
      %--------------------------------------------------------------------
      % return a vector of ingredient concentrations 
      %--------------------------------------------------------------------
         ingconcvec = 0;                        %default return
         if ~isempty(obj.ingABDvar)             %check for ingredients
            ingconcvec = [obj.ingABDvar(:).conc]; %sum conc values
         end
      end 

      function acidconcvec = get_acidconcvec(obj)
      %--------------------------------------------------------------------
      % return a vector of ingredient concentrations 
      %--------------------------------------------------------------------
         acidconcvec = [obj.acidABDvar(:).conc]; %acid vec always present
      end    
      
   end %end of methods

end