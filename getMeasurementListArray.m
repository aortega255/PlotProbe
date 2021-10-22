function [ml,rho] = getMeasurementListArray( data, probe )

nMeas = length(data.measurementList);

ml = zeros(nMeas,4);
rho = zeros(nMeas,1);

for iM = 1:nMeas
    ml(iM,1) = data.measurementList(iM).sourceIndex;
    ml(iM,2) = data.measurementList(iM).detectorIndex;
    if ~isempty(data.measurementList(iM).dataTypeIndex)
        ml(iM,3) = data.measurementList(iM).dataTypeIndex;
    else
        ml(iM,3) = 0;
    end
    ml(iM,4) = data.measurementList(iM).wavelengthIndex;
    
    ps = probe.sourcePos3D( ml(iM,1), : );
    pd = probe.detectorPos3D( ml(iM,2), : );
    rho(iM) = norm( ps-pd);
end 
