%more graphs
ix=0;
for i=sensi1
    phytom2=yy.(sprintf('case%d',i))(end-2*365+1:end,1:length(p.xgrid));
    ix=ix+1;
    phyto(ix,:)=sum(phytom2,2);
end
    %%
    hold on
    plot(t(1:2*365), phyto(1,:), 'k','linewidth',4)
       plot(t(1:2*365), phyto(2,:),'linewidth',4)
       plot(t(1:2*365), phyto(3,:),'linewidth',4)
       plot(t(1:2*365), phyto(4,:),'linewidth',4)
        hold on
        
    hold off
    xlabel('time (days)')
    ylabel('phytoplankton concentration (μMN m^-^2)')
   
   
    %hleg=legend('0','1','8','13','19','24', '29', '35','40','45','51','56','62');
    hleg=legend('0','7','42','426')
    htitle = get(hleg,'Title');
set(htitle,'String','g m^-^2 fish')
    set(gca,'fontsize',20)