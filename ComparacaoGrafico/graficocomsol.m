data = readtable('dados.txt', 'HeaderLines',0); % Pula a primeira linha do cabeçalho
frequency = data.Frequence;
value = data.Value;

data2 = readtable('output.txt', 'HeaderLines',0); % Pula a primeira linha do cabeçalho
frequency2 = data2.freq;
value2 = data2.auto_valeur;

% Filtrar linhas em que a coluna 'Frequence' é menor ou igual a 1000
filteredData = data(data.Frequence <= 1000, :);

% Gerar gráfico com os dados filtrados
figure;
plot(filteredData.Frequence, filteredData.Value, 'o');
title('Gráfico de Frequência vs. Valor (Filtrado)');
xlabel('Frequência');
ylabel('Valor');
grid on;
hold on;
plot(frequency2,value2, 'x');

