function [zbins,zadcp1,offset,x_null]=adcp_surface_fit(zadcp,ea,surface_bins,blen,tlen,lag,blnk,nbin,fbind)
    
    % Bin depth matrix
    dpt1   = repmat(zadcp,nbin,1);
    binmat = repmat((1:nbin)',1,length(dpt1));   
    
    z(1,:) = dpt1(1,:)-fbind;
    for ii = 2:length(binmat(:,1))
        z(ii,:) = z(1,:)-(binmat(ii,:)-1.5)*blen; 
    end
    
    % Loop over time, determine bin of maximum ea in surface bin range and 
    % do quadratic fit over 2 neighbouring and center bins   
    easurf = ea(surface_bins,:);
    for ii=1:length(ea)
        [eamax,maxind] = max(easurf(:,ii));
        
        if length(eamax)>1 
            maxind = maxind(1); 
        end
        
        if maxind>1 && eamax>80
            if surface_bins(maxind)==nbin 
                xtofit(1:2) = easurf(maxind-1:maxind,ii);
                xtofit(3)   = 0;
            else           
                xtofit      = easurf(maxind-1:maxind+1,ii);
            end

            for jj=1:3
                A(jj,:) = [(surface_bins(maxind)+jj-2)^2, surface_bins(maxind)+jj-2, 1];
            end 
            coef(:,ii) = A\xtofit;
        else
            coef(:,ii) = NaN;
        end
    end
    
    % Find maximum of quadratic fit ax^2+bx+c: 2ax+b=0
    x_null = -coef(2,:)./2./coef(1,:);

    %% Calculate offset
    offset = round(((x_null-1)*blen+blen/2+fbind)-(zadcp));     
    disp('-------------------------------');
    disp(['Depth offset is around ' num2str(round(nanmean(offset))) ' m']);
    disp('-------------------------------');
    % offset over time cleaned (median filter)
    [offset_clean,~] = clean_median(offset,20,2.8,[0.5 5],2,NaN);
    lin_offset       = linspace(1,length(offset),length(offset));
    offset           = interp1(lin_offset(~isnan(offset_clean)),...
        offset_clean(~isnan(offset_clean)),lin_offset,'linear','extrap');

    % offset median
    %offset = nanmedian(offset);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Plot histogram of differences
    dz     = ((x_null-1)*blen+blen/2+fbind)-zadcp;
    count  = [-100:1:100];
    ncount = hist(-dz,count);
    ncount = ncount/length(dz)*100;

    figure;
    bar(count,ncount);  
    axis([-50 50 0 100])
    xlabel('Depth difference [m]')
    ylabel('Occurrence percentage [%]')
    grid on
    title('Histogram of differences between depth record calculated with pressure sensor and from surface detection');
    
    figure(6);
    plot(-zadcp,'b');
    grid on
    hold on;
    plot(-((x_null-1)*blen+blen/2+fbind),'r');
    legend('Original','Reconstructed from surface reflection');
    
    if abs(offset)>15
       reply = input('Do you want to overwrite offset? 1/0 [0]:');
       if isempty(reply)
          reply = 0;
       end
       
       if reply==1
           offset=input('->Enter new offset:');
       else
            disp(['->Cleaned median filter offset is applied']);          
       end
    end
    
    %% Apply offset to get correct bin depth and instrument depth:
    zbins  = z+offset;
    zadcp1 = zadcp+offset;
    
    figure(6);
    plot(-zadcp1,'y');
    legend('Original','Reconstructed from surface reflection','Offset applied');   
    xlabel('Time Index')
    ylabel('Depth [m]')
    %print -dpng surface_fit;
    
    figure;
    pcolor([1:length(x_null)],-zbins,ea); shading flat;
    title('Amplitude'); colorbar; ylabel('Depth [m]');xlabel('Time index');
    %print -dpng surface_ea;

end
