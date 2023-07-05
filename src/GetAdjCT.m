function res = GetAdjCT(targetpH, ABtable, IS)
%FB version 03222019
%Optimize targt pH solution from observed AB mat by adjusting adjC. 
%PARAMETERS
%   targetpH: scaler, mean of init pH for acid and base titrations
%   ABtable: Matlab Table, table headings of conc, pK, ab.
%       Conc - concentration of buffer in molar (M)
%       pK - pK value for buffer 
%       ab - char array of 'a' or 'b' or 'x' (acid, base, placeholder)
%   IS: ionic strength (scalar)
%RETURNS
%   Matlab structure including:
%   res.adjC: new adjC value
%   res.Pred_pH: predicted pH
%   res.sqerr: error on optimized adjC
%   res.ExitFlag = exitflag (scalar)            
%   res.output = matlab fminsearch stats (array)
%DEPENDENCIES: 
% pH = CalcpH_AB(ABmatrix, NaClpercent, AdjC)
%   ABmatrix = Nx2, conc, pH
%   NaClpercetn = percent NaCl (scalar) 
%   AdjC = M/L HCl or NaCl (scalar)
%   RETURNS: pH of solution (scalar)
%*************************************************************************

%set up parameters for fmincon optimization
adjC = 0;                   %initialize return value

%function handle for objective func
fhan = @(adjC)CalcError(adjC,targetpH,ABtable,IS); 
%set optimization options
% optionvals = optimset('Display','none','MaxFunEvals', 1e9, ...
%    'TolFun', 1e-9, 'TolX', 1e-9);
optionvals = optimset('Display','none');
%Run fminsearch
[new_adjC, sqerr, exitflag, output] = fminsearch(fhan,adjC,optionvals); 
%calculate pH value for the table with new adjC value
temp = CalcpH_ABT(ABtable,IS,new_adjC);


%output results
res.adjC = new_adjC;
res.pred_pH = temp;
res.sqerr = sqerr; 
res.exitflag = exitflag;
res.output = output;

%calc error for minimization in fminsearch
function sqerr = CalcError(adjC,target_pH,ABtable,IS)
   phpred = CalcpH_ABT(ABtable,IS,adjC); %predicted pH
   sqerr = (target_pH - phpred)*(target_pH - phpred); %squared error
end
      
end

