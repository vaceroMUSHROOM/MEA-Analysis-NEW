function [timeVector,dataout] = getCorrectTStamps(ts,ds)

%tsrx=[];

%stampadd=ones(size(ts));
%tsx=ts(:);
%dsx=ds(1,:);
dsxo=ds(:);
dataout=dsxo; %dataout=[dataout dsxo];

% if isempty(tsrx)
%     for i=2:512
%         % ibite=i*31.25;
%         %%% offset=32.50*i; % I think this works with a different sampling
%         %%% frequency
%         offset = 
%         timenow=ts+stampadd.*offset;
%         timenow=timenow(:);
%         tsx=[tsx;timenow];
%     end
%     [reordered,~]=sort(tsx);
%     tsrx=reordered;
% end
% timeVector=tsrx;


if(length(dataout)/length(ts) == 512)
    timeVector = zeros(size(dataout));
    %count = 1;
    for i=1:(length(dataout)/512)
        if(i<length(dataout)/512)
            tp = linspace(ts(i),ts(i+1),513)';
        else
            tp = linspace(ts(i),ts(i)+1600,513)';
        end
        
       % timeVector(count:count+512-1,1) = tp(1:length(tp)-1);
       if(i==1)
          timeVector(i:512,1) = tp(1:512); 
       else
       timeVector(((i-1)*512)+1:i*512,1) = tp(1:512);
      % count = count + 512-1;
       end
       
       
    end
    
else
    disp('Problems with the data structures, please verify data files'); 
end



end