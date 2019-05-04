In = rand(1,256);       %%%%%%Normalization

In = (((In - min(In(:))).*(100 - 10))./(max(In(:))-min(In(:))))+10;

Imin = min(In(:));    %%%%%Reference Dataset generation
Imax = max(In(:)); 
dev = round((Imax - Imin)./3);
dev1 =2.*dev;
pv = []; tv=[];
pv1(1,:) = Imin:dev+Imin-1;
tv1(1:length(pv1)) = 1; 
pv2(1,:) = dev+Imin:dev1+Imin-1;
tv2(1:length(pv2)) = 2; 
pv3(1,:) = dev1+Imin:Imax-1;
tv3(1:length(pv3)) = 2; 

pv = [pv pv1 pv2 pv3];
tv = [tv tv1 tv2 tv3];

 s1 = 100; s2 = 1;   %%%%%%%%NN Creation and training  
 
 tf1 = 'tansig'; tf2 = 'purelin';
 
 btf = 'trainrp'; blf = 'learngd';
 
 pf = 'mse';
 
 net = newff(pv,tv,[s1 s2],{tf1 tf2},btf,blf,pf);
 
 net.trainParam.epochs = 100;
 net.trainParam.goal = 0.00001;  
 
 net = train(net,pv,tv);
