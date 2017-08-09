clear valid_event
a=table2array(Ped_SortEvent);
n=0;   %initial number of valid event
for k=1:max(a(:,2))
    b=find(a(:,2)==k);
    temp=a(b,:);   %the k-th event
    max_tran=max(temp(:,10));
    min_tran=min(temp(:,10));
    ped_vel_max=max(abs(diff(temp(:,10)*10)));   
    if min_tran<=-2&max_tran>=2&ped_vel_max<=8    %then this is a valid crossing
        n=n+1;
        valid_event(n).alldata=temp;
        d=sqrt(temp(:,8).^2+temp(:,10).^2);
        [min_d,min_index]=min(d);
        valid_event(n).min_distance=min_d;
        valid_event(n).veh_velocity=temp(min_index,19);
        s_t=size(temp,1);
        if min_index==1
           valid_event(n).ped_velocity=(temp(2,10)-temp(1,10))*10;
        end
        if min_index==s_t
           valid_event(n).ped_velocity=(temp(s_t,10)-temp(s_t-1,10))*10;
        end
        if min_index>1&min_index<s_t
           valid_event(n).ped_velocity=(temp(min_index+1,10)-temp(min_index-1,10))*5;
        end
    end
end
            
        