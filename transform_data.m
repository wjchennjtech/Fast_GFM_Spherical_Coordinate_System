function[data]=transform_data(datalow,dataup,deltarho)

data(:,[1,2])=dataup(:,[1,2]);
data(:,3)=datalow(:,3)*1e+3;
data(:,4)=dataup(:,3)*1e+3;
data(:,5)=(2.67-deltarho(:,3))*1e+3;

