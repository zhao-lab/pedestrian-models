clear;clc;close all;
FPlot = 0;
CurrentPath = fileparts(mfilename('fullpath'));
PedestrianData = readtable('Pedestrian_CCTC_north.csv');

Ped = PedestrianData(PedestrianData.GpsHeading>180&PedestrianData.Device>17000&PedestrianData.Longitude<-83.73543,:);

PedLon = Ped.Range.*sin(deg2rad(Ped.GpsHeading))-Ped.Transversal.*cos(deg2rad(Ped.GpsHeading));
PedLat = Ped.Range.*cos(deg2rad(Ped.GpsHeading))+Ped.Transversal.*sin(deg2rad(Ped.GpsHeading));

GainLon = 1e-9/(1e3*deg2km(distance(42.278506000, -83.735964000,42.278506000, -83.735964000+1e-9))); 
GainLat = 1e-9/(1e3*deg2km(distance(42.278506000, -83.735964000,42.278506000+1e-9, -83.735964000)));

Ped.PedLon = Ped.Longitude+PedLon*GainLon;
Ped.PedLat = Ped.Latitude+PedLat*GainLat;

ZeroPoint = [-83.735580, 42.278414];
cost = cos(deg2rad(-46.5));
sint = sin(deg2rad(-46.5));

Ped.BusX = cost.*(Ped.Longitude-ZeroPoint(1))./GainLon+sint.*(Ped.Latitude-ZeroPoint(2))./GainLat;
Ped.BusY = -sint.*(Ped.Longitude-ZeroPoint(1))./GainLon+cost.*(Ped.Latitude-ZeroPoint(2))./GainLat;
Ped.PedX = cost.*(Ped.PedLon-ZeroPoint(1))./GainLon+sint.*(Ped.PedLat-ZeroPoint(2))./GainLat;
Ped.PedY = -sint.*(Ped.PedLon-ZeroPoint(1))./GainLon+cost.*(Ped.PedLat-ZeroPoint(2))./GainLat;

Ped.TTC = (Ped.BusX-Ped.PedX)./(Ped.Speed);


Ped = Ped(Ped.PedX>-20&Ped.PedY>-10&Ped.BusY>-5&Ped.TargetId==1,:);%only one target


%% Add event index
PedSort = sortrows(Ped,{'Device','Trip','Time'});
Nline = height(PedSort);
PedSort.EventCombine = zeros(Nline,1);
nEventCombine = 1;
iEC = 0;
DeviceVec = unique(PedSort.Device)';
MinTimeSeq = 2; % [s]
clear EventCombineTable;
for nDevice = DeviceVec  % loop device
    PedTemp = PedSort(PedSort.Device==nDevice,:);
    TripVec = unique(PedTemp.Trip)';
    for nTrip = TripVec   % loop trip
        PedTemp2 = PedTemp(PedTemp.Trip==nTrip,:);
        iJump = find(diff(PedTemp2.Time)>10*MinTimeSeq);
        iSeq = [[1;(iJump+1)],[iJump;height(PedTemp2)]];
        NSeq = size(iSeq,1);  
        for nSeq = 1:NSeq  % loop pedestrian
            PedSort.EventCombine(iEC+(iSeq(nSeq,1):iSeq(nSeq,2)),1) = nEventCombine;
            EventCombineTable(nEventCombine,1:2) = [iEC+iSeq(nSeq,1),iEC+iSeq(nSeq,2)];
            nEventCombine = nEventCombine+1;  % 6113 events
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
%% screen short events
Duration_EventCombine = (EventCombineTable(:,2)-EventCombineTable(:,1))/10; % Time [s/10]
vpedestrian = zeros(length(Duration_EventCombine),1); %pedestrian speed screen
% for i = 1:length(Duration_EventCombine)
%     PedData.vp(EventCombineTable(i,1):EventCombineTable(i,2),1)=ones(EventCombineTable(i,2)-EventCombineTable(i,1)+1,1).*(PedData.PedY(EventCombineTable(i,2))-PedData.PedY(EventCombineTable(i,1)))./Duration_EventCombine(i);
%     vpedestrian(i) = PedData.vp(EventCombineTable(i,1));
% end
PedData.vp = zeros(size(PedData,1),1);
for i = 1:nEvent-1
    PedTemp = PedData(PedData.Event==i,:);
    vptemp = (PedTemp.PedY(size(PedTemp,1))-PedTemp.PedY(1))/(PedTemp.Time(size(PedTemp,1))-PedTemp.Time(1))*100;
    PedData.vp(PedData.Event==i)=ones(size(PedTemp,1),1)*vptemp;
end

Duration_Event = arrayfun(@(n) sum(PedSort.Event==n)/10,1:PedSort.Event(end));

nEventTooShort = find(Duration_Event>=0.5)'; %I think this means events acceptable

PedDataScreen = PedData(ismember(PedSort.Event,nEventTooShort),:);
sprintf('%d passing events, %d pedestrians',[length(unique(PedDataScreen.EventCombine)),length(unique(PedDataScreen.Event))])

PedDataScreen.Transversal = -PedDataScreen.Transversal;

%TranReverse = find(PedDataScreen.vp<0);

%PedDataScreen.cd=PedDataScreen.Transversal+PedDataScreen.vp.*PedDataScreen.TTC;

%PedDataScreen.Transversal(TranReverse) = -PedDataScreen.Transversal(TranReverse);

%PedDataScreen.cd(TranReverse) = -PedDataScreen.cd(TranReverse);


%PedDataScreen.vp = abs(PedDataScreen.vp);



%% machine learning
%PedDataScreen = PedDataScreen(PedDataScreen.Range>2&PedDataScreen.Range<30&PedDataScreen.Transversal>-15&PedDataScreen.Transversal<15,:);
%PedDataScreen = PedDataScreen(PedDataScreen.Transversal>-15&PedDataScreen.Transversal<15,:);

 
%writetable(PedDataScreen,'FitPEDnorthCORD.xlsx')
%% mean value
% MeanTable=zeros(nEventComb,4);
% for ievent = 1:nEventComb
%     DataTemp = PedDataScreen(PedDataScreen.EventCombine==ievent,:);
%     MeanTable(ievent,:) = [mean(DataTemp.TTC),mean(DataTemp.Range),mean(DataTemp.Transversal),abs(mean(DataTemp.vp))];
% end
% EventCombine = (1:nEventComb)';
% TTC = MeanTable(:,1);
% Range=MeanTable(:,2);
% Transversal = MeanTable(:,3);
% vp= MeanTable(:,4);
% MeanTable=table(EventCombine,TTC,Range,Transversal,vp);

%eventa=find(MeanTable.TTC>5&MeanTable.TTC<10&MeanTable.Range>3&MeanTable.Range<10&MeanTable.Transversal>0&MeanTable.Transversal<5&MeanTable.vp>1&MeanTable.vp<1.5);
% eventa=find(MeanTable.TTC>10&MeanTable.TTC<15&MeanTable.Range>3&MeanTable.Range<10&MeanTable.Transversal>-5&MeanTable.Transversal<0&MeanTable.vp>0&MeanTable.vp<1);
% for i = eventa'
%     DataTemp = PedDataScreen(PedDataScreen.EventCombine==i,:);
%     Time = DataTemp.Time-ones(height(DataTemp),1)*DataTemp.Time(1);
%     plot(Time,DataTemp.Ax);
%     hold on;
% end

%% plot the events within in a combine
EventCombineVec_All = unique(PedDataScreen.EventCombine)';
NEventCombineVec_All = length(EventCombineVec_All);


%% Pedestrians are nearest objects
