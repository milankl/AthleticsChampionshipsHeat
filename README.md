[![DOI](https://zenodo.org/badge/262036919.svg)](https://zenodo.org/badge/latestdoi/262036919)
# Heat indices during European and World Championships in Athletics
Universal thermal heat indices during some World and European championships in athletics.

This repository contains meteorological data and analysis scripts for

K Hollander, M Kloewer, A Richardson, L Navarro, S Racinais, V Scheer, A Murray, P Branco, T Timpka, A Junge and P Edouard, 2021. *Apparent temperature and heat-related illnesses during international athletic championships: A prospective cohort study*. **Scandinavian Journal of Medicine and Science in Sports**, accepted.

and is used for daily heat forecasts at the Tokyo Olympic Games, tweeted daily via [twitter.com/HeatForecast](https://twitter.com/HeatForecast)
## Events and stadium coordinates

| City | Year | Date | Location | Championship | UTC |
|------|------|------|----------|--|--|
|Berlin | 2009 | 15-23 Aug | 52.5°N, 13.2°E|World | +2 |
|Daegu  | 2011 | 27 Aug - 4 Sep | 35.8°N, 128.7°E|World | +9 | 
|Helsinki | 2012 | 27 Jun - 1 Jul | 60.2°N, 24.9°E|European | +3 |
|Moscow   | 2013 | 10-18 Aug | 55.7°N, 37.6°E|World | [+4](https://en.wikipedia.org/wiki/Moscow_Time) |
|Zürich   | 2014 | 12-17 Aug | 47.4°N, 8.5°E|European | +2 |
|Amsterdam| 2016 | 6-10 Jul  | 52.3°N, 4.9°E|European | +2 |
|Berlin   | 2018 | 6-12 Aug  | 52.5°N, 13.2°E|European | +2 |
  
## Data 

From [ERA5-Land hourly data](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview) at 0.1° resolution

 - 2m temperature (t2m)
 - 2m dew point temperature (d2m)
 - 10m u-component of wind (u10)
 - 10m v-component of wind (v10)
 - Surface net solar radiation (ssr)
 - Surface net thermal radiation (str)
 - Surface pressure (sp)
 - Surface solar radiation downwards (ssrd)
 - Surface thermal radiation downwards (strd)
 - Total precipitation (tp)

From [ERA5 hourly data](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview) at 0.25° resolution

  - High, medium and low cloud cover (hcc, mcc, lcc)
  - Total sky direct solar radiation at surface (fdir)

The 0.25° ERA5 data was bilinearly interpolated on ERA5-Land's 0.1°-grid.
The spatial resolution is only used for uncertainty estimation.
The closest grid point to the stadium is used.

### The Universal Thermal Climate Index

The [universal thermal climate index (UTCI)](http://utci.org/) is calculated with the polynomial approximation from the [pythermalcomfort](https://github.com/CenterForTheBuiltEnvironment/pythermalcomfort) package. The 2m temperature and wind speed (=sqrt(u^2+v^) are directly used from ERA5-Land, the 2m dew point temperature is converted to relative humidity using the [August-Roche-Magnus formula](https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation#August-Roche-Magnus_approximation). The [mean radiant temperature](https://en.wikipedia.org/wiki/Mean_radiant_temperature) is calculated following [Di Napoli, 2020](https://link.springer.com/article/10.1007%2Fs00484-020-01900-5) from surface solar and thermal radiation (net and downward) from ERA5-Land and ERA5.

### Example

![meteogram](plots/berlin2009.png?raw=true "Meteogram")
