N = 1020;
m_max = 2*10*N;
rep = 10000;

K = zeros(m_max,rep);

for m=1:m_max
    for r=1:rep
        n = sum(randi(2,[1 m])==1);
	%n=m;
        K(m,r) = numel(unique(randi(N,[1 n])));
    end
    disp(m);
end

save(['Half Sim 2 N' num2str(N) ' m_max' num2str(m_max) ' rep' num2str(rep) '.mat'],'-v7.3');
% quit

%%

Kt = cell(N+1,1);
Km = zeros(N+1,1);
Ks = zeros(N+1,1);
ba = zeros(N+1,1);

parfor k=0:N
    [m ~] = find(K==k);
    Kt{k+1,1} = m;
    Km(k+1,1) = mean(m);
    Ks(k+1,1) = std(m);
    ba(k+1,1) = (log(1 - (k/N))/log(1 - 1/N));
end

%%

ave = zeros(1,N+1);
cv = zeros(1,N+1);
lcb = zeros(1,N+1);
ucb = zeros(1,N+1);


for k=0:N
    rng = min(Kt{k+1}):max(Kt{k+1});
    df = hist(Kt{k+1},rng);
    df = df./sum(df);
    cudf = cumsum(df);
    
    ave(k+1) = sum(rng.*df);
    cv(k+1) = sqrt(abs(sum(df.*rng.^2) - (sum(df.*rng))^2))/(sum(df.*rng));
    
    indl = rng(find(cudf>0.025));
    if(indl(1) == 1)
        l = indl(1);
    else
        l = indl(1)-1;
    end;
    
    indu = rng(find(cudf<0.975));
    if(length(indu) == 0)
        u = l;
    elseif(indu(length(indu))==m_max)
        u = indu(length(indu));
    else
        u = indu(length(indu))+1;
    end;
    
    lcb(k+1) = l;
    ucb(k+1) = u;
end

ba = ba';

%%

Out = [0:N;2.*ba;ave;cv;lcb;ucb]';
outid=fopen(['N0-' num2str(N) '_m_max' num2str(m_max) '_stats.csv'],'w+');
fprintf(outid,'%s',['k,BA,E(m),CV,LCB,UCB']);
fclose(outid);
dlmwrite(['N0-' num2str(N) '_m_max' num2str(m_max) '_stats.csv'],Out,'roffset',1,'-append');

%%
close all
clear
load('HalfSampleStats.mat')
Have = ave;
Hcv = cv;
Hlcb = lcb;
Hucb = ucb;
load('1020ChamberStats.mat')


plot(1:N,Have,'k.',1:N,Hlcb,'b.',1:N,Hucb,'.')
set(gca,'YLim',[0 21000])
set(gca,'XLim',[0 1050])
xlabel('Number of hits k');
ylabel('Number of molecules E(m)')
axis square

%%
close all
plot(1:N,Hcv,'k.',1:N,cv,'b.')
axis square
xlabel('Number of hits k');ylabel('CV')
set(gca,'XLim',[0 1050])

