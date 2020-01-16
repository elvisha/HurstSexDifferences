function [alpha, FL_all] = dfa_fast(vdata, istart, iend, L_all)
%function takes in your time series, the start and end time points, and the
%different L values you want to use in the implementation of DFA

%takes the cumulative sum
vdata=cumsum(vdata);

FL_all=zeros(numel(L_all),size(vdata,2));

%iterating through the L values you want to use
for il=1:numel(L_all)
    L=L_all(il);
    X=[[0:L-1]' ones(L,1)];

    %nice thing about this approach is if your data isn't an integer
    %multiple of the length L, it will just average as many windows as can fit
    c=0;
    FL=zeros(1,size(vdata,2));
    for i = istart:L:min(iend,size(vdata,1))-L+1
        vtmp=vdata(i+(0:L-1),:);

        %b=X\vtmp;
        %y=X*b;
        %r=vtmp-X*(X\vtmp);
        %calculates rms for that window
        rms=sqrt(mean((vtmp-X*(X\vtmp)).^2,1));
        FL=FL+rms;
        c=c+1;
    end

    FL=FL/c;
    FL_all(il,:)=FL;
    
end

logFL=log(FL_all);
X=L_all(:);
logX=[log(X) ones(size(X,1),1)];
b=logX\logFL;
alpha=b(1,:);