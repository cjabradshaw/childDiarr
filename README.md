# Predictors of diarrhoea in children under five years old
Determine the most important predictors of diarrhoea in children under five in Southeast Asia by exploring the spatiotemporal association between diarrhoeal incidence and various behavioural, socio-demographic, and environmental factors.
<img align="right" src="www/stomachache.jpg" alt="" width="400" style="margin-top: 20px">
<br>
<br>
Dr <a href="https://globalecologyflinders.com/people/#DIRECTOR">Syeda Hira Fatima</a><br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2024<br>
<a href=mailto:syeda.fatima@flinders.edu.au>e-mail</a> <br>
<br>
and<br>
<br>
Prof <a href="https://globalecologyflinders.com/people/#DIRECTOR">Corey J. A. Bradshaw</a> <br>
<a href="http://globalecologyflinders.com" target="_blank">Global Ecology</a> | <em><a href="https://globalecologyflinders.com/partuyarta-ngadluku-wardli-kuu/" target="_blank">Partuyarta Ngadluku Wardli Kuu</a></em>, <a href="http://flinders.edu.au" target="_blank">Flinders University</a>, Adelaide, Australia <br>
September 2024<br>
<a href=mailto:corey.bradshaw@flinders.edu.au>e-mail</a> <br>
<br>

## <a href="https://github.com/cjabradshaw/childDiarr/tree/main/scripts">Scripts</a>
- <code>DHSDiarrAnalysis.R</code>: R code to reproduce the resampled boosted regression tree analysis for determining the relationships between probability of diarrhoea, and socio-economic, maternal, child, climate data.

## <a href="https://github.com/cjabradshaw/childDiarr/tree/main/data/brtdata">Data</a>
- <em>DHSclusterLevelDiarrData.csv.zip</em>: <a href="https://dhsprogram.com/data/">Demographic Health Survey</a> data summarised by cluster with central parameter (mean, proportion, etc.) and variance per cluster. Overlaid (cluster-level) climate data derived from <a href="https://www.worldclim.org/">WorldClim</a> <a href="https://www.worldclim.org/data/bioclim.html">bioclimatic variables</a> (mean annual temperature, temperature annual range, total annual precipitation, precipitation seasonality, and precipitation of the driest quarter).

## Required R libraries
- <code>dismo</code>, <code>dplyr</code>, <code>ggplot2</code>, <code>ggpubr</code>, <code>gbm</code>, <code>spatstat.random</code>, <code>tidyr</code>, <code>truncnorm</code>, <code>usdm</code>

<br>
<br>
<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Horizontal_RGB_Master.png" alt="Flinders University" width="150" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="Global Ecology Lab" width="85" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.uwa.edu.au/"><img align="bottom-left" src="www/UWA.png" alt="UWA" width="100" style="margin-top: 20px"></a> &nbsp; &nbsp; <a href="https://www.thekids.org.au"><img align="bottom-left" src="www/TheKids-Logo.png" alt="The Kids Research Institute" width="90" style="margin-top: 20px"></a>
