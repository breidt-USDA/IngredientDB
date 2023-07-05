function tBetaData = get_tBetaData(ABmat,pHmin,pHmax,waterIS,crvIS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
waterNaClpct = waterIS*5.84;
crvNaClpct = crvIS*5.84;
pHvec = pHmin:0.05:pHmax;
BCcrv = BetaModel_AB(ABmat,crvNaClpct,pHvec');
Watercrv = BetaModel_AB([0,0],waterNaClpct,pHvec');
waterData = SCBCfit_area(Watercrv,15);
modelData = SCBCfit_area(BCcrv,15);
tBetaData.tBeta = 100*(modelData.area - waterData.area);
tBetaData.waterNaClpct = waterNaClpct;
tBetaData.crvNaClpct = crvNaClpct;
tBetaData.pHvec = pHvec;
tBetaData.BCcrv = BCcrv;
tBetaData.Watercrv = Watercrv;
tBetaData.WatertBeta = 100*waterData.area;
tBetaData.ModeltBeta = 100*modelData.area;
% for debugging
% plot(pHvec,Watercrv(:,2),'--k',pHvec,BCcrv(:,2),'-k');
% hold on;
% area( pHvec,BCcrv(:,2),'FaceColor',[0.8 0.8 0.8]); %area light gray
% alpha( 0.5);                               %area is transparant
% area( pHvec,Watercrv(:,2),'FaceColor',[1 1 1]);
% xlim([2,12]);
% alpha( 0.5);
% hold off
end