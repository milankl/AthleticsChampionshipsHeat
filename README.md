# AthleticsChampionshipsHeat
Universal thermal heat indices during some World and European championships in ahtletics.

## Events and stadium coordinates

| City | Year | Date | Location | Championships |
|------|------|------|----------|--|
|Berlin | 2009 | 15-23 Aug | 52.5°N, 13.2°E|World |
|Daegu  | 2011 | 27 Aug - 4 Sep | 35.8°N, 128.7°E|World |
|Helsinki | 2012 | 27 Jun - 1 Jul | 60.2°N, 24.9°E|European |
|Moscow   | 2013 | 10-18 Aug | 55.7°N, 37.6°E|World |
|Zürich   | 2014 | 12-17 Aug | 47.4°N, 8.5°E|European |
|Amsterdam| 2016 | 6-10 Jul  | 52.3°N, 4.9°E|European |
|Berlin   | 2018 | 6-12 Aug  | 52.5°N, 13.2°E|European |
  
## Data 

From [ERA5-Land hourly data](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=overview) at 0.1° resolution

 - 2m temperature
 - 2m dew point temperature
 - 10m u-component of wind
 - 10m v-component of wind
 - 2m dewpoint temperature
 - 2m temperature
 - Surface net solar radiation
 - Surface net thermal radiation
 - Surface pressure
 - Surface solar radiation downwards
 - Surface thermal radiation downwards
 - Total precipitation

From [ERA5 hourly data](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview) at 0.25° resolution

  - High cloud cover
  - Low cloud cover
  - Medium cloud cover
  - Total sky direct solar radiation at surface

The 0.25° ERA5 data was bilinearly interpolated on ERA5-Land's 0.1°-grid. The spatial resolution is only used for uncertainty estimation. The closest grid point the stadium is used.


### The Universal Thermal Climate Index (UTCI)

![meteogram](plots/berlin2009.png?raw=true "Meteogram")
