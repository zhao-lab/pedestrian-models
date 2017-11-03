function pattern = em_censor_pattern(x, xcb)
% finds the censoring pattern of x 
% 
% x         N*d     censored observations
% xcb       2*d     censoring bounds [lower ; upper]
% 
% pttn      N*d(0)  1: observed, 0: censored
% cmpl      N*1     1: complete
% xR        N*2(c)  ranges of censored elements
% 


[N,d] = size(x);


% lower/upper censoring bounds
if size(xcb,1)==1,  xcb = xcb';             end
if size(xcb,2)<d,   xcb = xcb*ones(1,d);    end
xl = xcb(1,:);
xu = xcb(2,:);


% censored coordinates 
% 1: uncensored, 0: censored
censored = bsxfun(@gt, x, xl) & bsxfun(@lt, x, xu);
% completely observed instances
cmpl = all(censored,2);


% integration regions
xR = cell(N,2);
for n=1:N
    mn = ~censored(n,:);
    
    if any(mn)
        regl = xu(mn);
        regu = xl(mn);

        regl(x(n,mn) == xl(mn)) = -Inf;
        regu(x(n,mn) == xu(mn)) =  Inf;

        xR{n,1} = regl;
        xR{n,2} = regu;
    else
        xR{n,1} = [];
        xR{n,2} = [];
    end
end


% TODO
% for each coordinate, only uni-directional censoring is considered.
% bi-directional censoring is not handled 


% re-order the censoring pattern 
uniq_pttn = unique(censored, 'rows');
same_pttn = false(N, size(uniq_pttn,1));

for k=1:size(uniq_pttn,1)
    same_pttn(:,k) = all(bsxfun(@eq, censored, uniq_pttn(k,:)), 2);
end
pttn_cnt = sum(same_pttn);

[pttn_cnt, idx] = sort(pttn_cnt, 'descend');
uniq_pttn = uniq_pttn(idx,:);
same_pttn = same_pttn(:,idx);


% output
pattern.censored = censored;
pattern.complete = cmpl;
pattern.range = xR;
pattern.uniq_censored = uniq_pttn;
pattern.count = pttn_cnt;
pattern.xpttn = same_pttn;

end

