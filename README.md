# Resumo

Em estudos de sobrevivência, é comum a presença de indivíduos que não experimentam o evento de 
interesse após um longo do período de acompanhamento, configurando a chamada fração de cura. Além disso, 
em aplicações biomédicas, os tempos de sobrevivência são frequentemente observados em unidades discretas 
e inteiras, como dias, semanas ou meses, o que motiva o uso de modelos formulados diretamente nessa escala 
temporal. Neste trabalho, propõe-se um modelo de sobrevivência com fração de cura em tempo discreto baseado 
em riscos competitivos latentes. Assume-se que cada indivíduo possui um número latente de causas potenciais 
para a ocorrência do evento, sendo considerado curado quando esse número é nulo. O número de riscos latentes 
é modelado por uma distribuição Binomial Negativa, permitindo capturar heterogeneidade não observada entre 
indivíduos, enquanto os tempos de ativação associados a cada risco seguem uma distribuição Lognormal Discreta, 
proposta como alternativa às distribuições discretas mais usuais em análise de sobrevivência em tempo discreto, 
como aquelas baseadas na Weibull discreta, oferecendo maior flexibilidade na descrição de diferentes padrões de 
risco. A inferência estatística é conduzida via máxima verossimilhança, considerando censura à direita. 
O desempenho do modelo é avaliado por meio de um estudo de simulação de Monte Carlo sob diferentes cenários de 
tamanho amostral e níveis de censura. Por fim, o modelo é aplicado a um conjunto de dados reais de pacientes com 
infarto agudo do miocárdio, ilustrando sua aplicabilidade e interpretabilidade em contextos clínicos com tempos de 
sobrevivência naturalmente discretos.

    
### Palavras-chave: 
Análise de sobrevivência em tempo discreto; Modelo de fração de cura;
Distribuição binomial negativa; Riscos competitivos latentes; Distribuição Lognormal Discreta; 
Modelos de Regressão; Estimação por máxima verossimilhança; Simulação de Monte Carlo.
