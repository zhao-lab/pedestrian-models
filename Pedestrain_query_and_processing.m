 clear all
query=['SELECT  A.*, B.Latitude, B.Longitude,', ...
' B.GpsHeading,B.Ax,B.Ay,B.Speed', ...
' FROM [SpFot].[dbo].[DataFrontTargets] A', ...
' JOIN [SpFot].[dbo].[DataDas] B', ...
' ON A.Device = B.Device', ...
' AND A.Trip = B.Trip', ...
' AND A.Time = B.Time', ...
' WHERE A.TargetType = 3', ...
' AND B.Latitude BETWEEN 42.278112 AND 42.278582', ...
' AND B.Longitude BETWEEN -83.735712 AND -83.735188'];
% ' AND B.Latitude BETWEEN 42.278333 AND 42.278749', ...
% ' AND B.Longitude BETWEEN -83.736307 AND -83.735681'];%, ...

DBName = 'SP_Ding';
Server = 'Tri-spdb1';

PedestrianData=SQLQuery_Ding(DBName,Server,query);

PedBus = PedestrianData;
Ped = PedestrianData;
PedLon = Ped.Range.*sin(deg2rad(Ped.GpsHeading))+Ped.Transversal.*cos(deg2rad(Ped.GpsHeading));
PedLat = Ped.Range.*cos(deg2rad(Ped.GpsHeading))-Ped.Transversal.*sin(deg2rad(Ped.GpsHeading));

GainLon = 1e-6/(1e3*deg2km(distance(42.278506, -83.735964,42.278506, -83.735964+1e-6)));
GainLat = 1e-6/(1e3*deg2km(distance(42.278506, -83.735964,42.278506+1e-6, -83.735964)));

Ped.PedLon = Ped.Longitude+PedLon*GainLon;
Ped.PedLat = Ped.Latitude+PedLat*GainLat;


%% Add event index
PedSort = sortrows(Ped,{'Device','Trip','Time'});
Nline = height(PedSort);
PedSort.EventCombine = zeros(Nline,1);
nEventCombine = 1;
iEC = 0;
DeviceVec = unique(PedSort.Device)';
MinTimeSeq = 2; % [s]
clear EventCombineTable;
for nDevice = DeviceVec
    PedTemp = PedSort(PedSort.Device==nDevice,:);
    TripVec = unique(PedTemp.Trip)';
    for nTrip = TripVec
        PedTemp2 = PedTemp(PedTemp.Trip==nTrip,:);
        iJump = find(diff(PedTemp2.Time)>100*MinTimeSeq);
        iSeq = [[1;(iJump+1)],[iJump;height(PedTemp2)]];
        NSeq = size(iSeq,1);
        for nSeq = 1:NSeq
            PedSort.EventCombine(iEC+(iSeq(nSeq,1):iSeq(nSeq,2)),1) = nEventCombine;
            EventCombineTable(nEventCombine,1:2) = [iEC+iSeq(nSeq,1),iEC+iSeq(nSeq,2)];
            nEventCombine = nEventCombine+1;
        end
        iEC = iEC+iSeq(NSeq,end);
    end
end

nEvent = 1;
Event = zeros(Nline,1);
for nEventComb = 1:PedSort.EventCombine(end)
    nEventVec = unique(PedSort.ObstacleId(PedSort.EventCombine==nEventComb,:))';
    for ntempEvent = nEventVec
        Event(PedSort.EventCombine==nEventComb&PedSort.ObstacleId==ntempEvent) = nEvent;
        nEvent = nEvent+1;
    end
end

PedSort.Event = Event;
PedData= PedSort(:,[1,end,end-1, 2:end-2]);

Ped_SortEvent=sortrows(PedData,{'Event','Time'});
%column 14-15: Bus GPS coordinate; column 20-21: Pedestrain coordinate