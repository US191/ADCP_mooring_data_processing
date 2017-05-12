function [zbins,zadcp1,offset,x_null]=adcp_surface_fit(zadcp,ea,surface_bins,blen,blnk,nbin);
    
    % Bin depth matrix
    dpt1 = repmat(zadcp,nbin,1);
    binmat = repmat((1:nbin)',1,length(dpt1));   
    z = dpt1-(binmat-0.5)*blen-blnk; 
    
    % Loop over time, determine bin of maximum ea in surface bin range and 
    % do quadratic fit over 2 neighbouring and center bins   
    easurf=ea(surface_bins,:);
    for ii=1:length(ea)
        [eamax,maxind]=max(easurf(:,ii));
        
        if length(eamax)>1; maxind=maxind(1); end
        
        if maxind>1 & eamax>80
            if surface_bins(maxind)==nbin; 
                xtofit(1:2) = easurf(maxind-1:maxind,ii);
                xtofit(3) = 0;
            else           
                xtofit = easurf(maxind-1:maxind+1,ii);
            end

            for jj=1:3
                A(jj,:)=[(surface_bins(maxind)+jj-2)^2, surface_bins(maxind)+jj-2, 1];
            end 
            coef(:,ii)= A\xtofit;
        else
            coef(:,ii)=NaN;
        end
    end
    
    % Find maximum of quadratic fit ax^2+bx+c: 2ax+b=0
    x_null = -coef(2,:)./2./coef(1,:);

    offset=-round(nmedian((x_null-0.5)*blen+blnk)-nmedian(zadcp));
    disp('-------------------------------');
    disp(['Depth offset is ' num2str(offset) ' m']);
    disp('-------------------------------');
    
    
    % Plot histogram of differences
    dz=((x_null-0.5)*blen+blnk)-zadcp;
    count=[-100:1:100];
    ncount=hist(-dz,count);
    figure(1);
    bar(count,ncount);
    title('Histogram of differences between original and reconstructed depth record');
    ylabel('Frequency');
    xlabel('Difference of Depth (m)');
    
    figure(2);
    plot(zadcp,'b');
    hold on;
    plot((x_null-0.5)*blen+blnk,'r');
    legend('Original','Reconstructed from surface reflection');
    
    if abs(offset)>15
       reply = input('Do you want to overwrite offset? 1/0 [0]:');
       if isempty(reply)
          reply = 0;
       end
       
       if reply==1
           offset=input('Enter new offset:');
       end
    end
    
    disp(['Offset of ' num2str(offset) ' m is applied']);

    % Apply offset to get correct bin depth and instrument depth:
    zbins=z-offset;
    zadcp1=zadcp-offset;
    
    figure(2);
    plot(zadcp1,'y');
    text(300, max(zadcp),['Offset applied: ' num2str(offset) ' m']);
    %,'fonts',12,'fontw','bold','backgroundc','w');
    legend('Original','Reconstructed from surface reflection','Offset applied');
    ylabel('Depth (m)');
    xlabel('Time index');
    
    %print -dpng surface_fit;
    
    figure(3);
    pcolor([1:length(x_null)],zbins,ea); shading flat;
    set(gca,'ydir', 'reverse');
    ylabel('Corrected Depth (m)');
    xlabel('Time index');
    title('Amplitude of the bins with corrected depth');

    %print -dpng surface_ea;

end