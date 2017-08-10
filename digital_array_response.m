clear all;%clc

N = 765;
m_max = 10*N;

% disp(N);
% disp(m_max);

pkm_arr = zeros(m_max,N);
pmk_arr = pkm_arr';

% resolve = max(1E4,ceil(5*m_max));
resolve = 1E4;
th_arr = linspace(-pi,pi,resolve);
th_arr(resolve) = [];

d = N+1;

if m_max < N
    error('Number of possible molecules must be larger than number of chambers in this constrained (but practical) implementation')
end

repint = 100;
ndis = round(m_max/repint);
t = zeros(1,ndis);
to = 1;
    
tickity = tic;

parfor m=1:m_max
    
%      if ~mod(m,repint)
%          t(to) = toc(tickity);
%          left = (mean(diff(t(1:to))).*(ndis-to))/60;
%          to = to+1;
%          disp([num2str(100*(m/m_max)) '% Est. ' num2str(left) ' min'])
%      end
    
    f = ones(size(th_arr)); % initialize
    omega = zeros(1,N);
    
    for k=1:N
        if m<k
            break
        end
        a = (N-k+1)./(d*exp(1i*th_arr)-k);
        f = f.*a;
        omega(k) = sum(f.*exp(1i*m*th_arr));
    end
    omega = real(omega);
    pkm_arr(m,:) = omega./sum(omega);
end

pkm_arr = abs((1E-9.*round(1E9.*pkm_arr)));


pkm_arr = sparse(pkm_arr);

for count=1:N
    % count
    pmk_arr(count,:) = pkm_arr(:,count)/sum(pkm_arr(:,count));
    % clc;
end
toc(tickity)

save(['N' num2str(N) ' m_max' num2str(m_max) '.mat'],'-v7.3');
% quit

%% calculate statistical properties
% This code can then be appended to calculate the relevant statistical measures discussed
% earlier:
vals = 1:m_max;
cv = zeros(1,N);
ave = zeros(1,N);
lcb = zeros(1,N);
ucb = zeros(1,N);
ba = zeros(1,N);
for(count = 1:N);
    cv(count) = sqrt(abs(sum(pmk_arr(count,:).*(vals).^2) - (sum(pmk_arr(count,:).*(vals)))^2))/(sum(pmk_arr(count,:).*(vals)));
    ave(count) = sum(pmk_arr(count,:).*vals);
    ba(count) = (log(1 - (count/N))/log(1 - 1/N));
    
    cudf = cumsum(pmk_arr(count,:));    
    
    indl = find(cudf>0.025);
    if(indl(1) == 1)
        l = indl(1);
    else
        l = indl(1)-1;
    end;
    
    indu = find(cudf<0.975);
    if(length(indu) == 0)
        u = l;
    elseif(indu(length(indu))==m_max)
        u = indu(length(indu));
    else
        u = indu(length(indu))+1;
    end;
    
    lcb(count) = l;
    ucb(count) = u;
end;

Out = [1:N;ba;ave;cv;lcb;ucb]';
outid=fopen(['N' num2str(N) '_m_max' num2str(m_max) '_stats.csv'],'w+');
fprintf(outid,'%s',['k,BA,E(m),CV,LCB,UCB']);
fclose(outid);
dlmwrite(['N' num2str(N) '_m_max' num2str(m_max) '_stats.csv'],Out,'roffset',1,'-append');
