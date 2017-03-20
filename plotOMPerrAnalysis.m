
% plotOMPerrAnalysis.m
%
%
% This plots OMPerrAnalysis.m outputs as a function of epsilon
% producing separate plots for each and saving them.
%
% Created: October 29, 2011 Ra Inta
% Last modified: October 29, 2011 R.I.

for i=1:7; 
    epsilon = 10^(i -6); 
    subplot(2,1,1)
    plot(1:150,s,'-r+',1:length(sRecon{i}),sRecon{i},'-ko')
    title(['OMP reconstruction error analysis (\epsilon = 10^{', num2str(i-6),'})'])
    legend('Sparse signal vector','OMP reconstruction')
    ylabel('Amplitude (arbitrary)')
    subplot(2,1,2)
    plot(Rrecon{i})
    ylabel('Residual error')
    saveas(gcf,['OMPerrAnalysis_epsilon1E', num2str(i-6)],'fig')
end



