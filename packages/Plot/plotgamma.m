function plotgamma(Str, Tk)
close all;
global wing Para GG1 GG2

str = Str;
load( [str,'parameter.mat'],'Para');
str=[str,'Ddata/'];

load( [str,'tem.mat'], 'tem_at');
st = Tk;

h=figure(1);
plotset(3.5,3);
for tk=st:1:5500%tem_at-1
    tk
    load( [str,'T',num2str(tk+1),'.mat'], 'wing');

    GG1(tk) = wing(1).gammaMax;
    GG2(tk) = wing(2).gammaMax;
    %plot(Para.s, wing(1).nu./sqrt(1-Para.s.^2),'b.-'); hold on;
    %plot(Para.s, wing(2).nu./sqrt(1-Para.s.^2),'r.-'); hold on;
    %%plot(Para.s, wing(1).nu,'b.-'); hold on;
    %%plot(Para.s, wing(2).nu,'r.-'); hold on;

    % plot the up down velocity on the wing, used for boundary layer calculation
%%    [vup, vdown]=boundaryspeed(1, tk);
%%    plot(Para.s, vup,'b--'); hold on;
%%    plot(Para.s, vdown,'r--'); hold on;
%%        [vup, vdown]=boundaryspeed(2, tk);
%%    plot(Para.s, vup,'b-'); hold on;
%%    plot(Para.s, vdown,'r-'); hold on;
%%    plot([-1, 1], [0, 0], 'k--');
%%
%%    xlabel('Wing');
%%    ylabel('Bound vortex sheet strength');
%%    title(['Bound vortex sheet strength, t=', num2str(Tk/100)]);
%%    ylim([-20, 20]);
%%    legend({'1st wing','2nd wing'}, 'Location', 'best', 'FontSize', 8);
%%    hold off;
%%    ps = .1;
%%    pause(ps)

%fig_name = [Str, 'gamma', num2str(tk),]
%saveas(h, [fig_name, '.fig']);
%eval(['print -dpng -r300 ' fig_name '.png']);
end
plot(GG1,'r'); hold on;
plot(GG2,'b')
grid on;

end
