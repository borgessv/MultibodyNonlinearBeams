function [] = stability_plot(A)
[Vvec,Dvec] = eig(A);
stability = sort(diag(Dvec));
plot(stability,'ro');
axis([-10 1 -400 400]);
%ax = gca;
%ax.YAxisLocation = 'origin';
%ax.XAxisLocation = 'origin';

for i = 1:length(stability)
    if real(stability(i)) > 0
        fprintf('Instability! \n')
        break
    end
end
end