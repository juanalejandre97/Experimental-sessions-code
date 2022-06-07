function [valores,localizaciones] = limpiapicos(data,corte)

n = length(data);
data_limpiao = zeros(1,n);

for i =1:n
    if (data(i) > corte)
        data_limpiao(i) = data(i);
    else
        data_limpiao(i) = 0;
    end
end
[valores,localizaciones] = findpeaks(data_limpiao);

end