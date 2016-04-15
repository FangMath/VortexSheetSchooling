function plotset(wdth, hght)
set(0,'DefaultTextFontSize',12);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontname', 'Arial')

%wdth = 5; hght = 3.5;
set(gcf,'position',[100,100,wdth*100,hght*100]);
set(gcf, 'PaperPosition', [0 0 wdth hght]); %Position plot at left hand corner with width 2 and height 6.5.
set(gcf, 'PaperSize', [wdth hght]); %Set the paper to have width 2 and height 6.5.
end
