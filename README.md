# Mapping urban economic centers worldwide
Shuai Pang, Junlong Zhang and Lei Dong

## Abstract

Urban centers serve as engines of regional development, yet accurately defining and identifying the socioeconomic centers of cities globally remains a big challenge. Existing mapping efforts are often limited to large cities in developed regions and rely on data sources that are unavailable in many developing countries. This data scarcity hinders the establishment of consistent urban indicators, such as accessibility, to assess progress towards the United Nations Sustainable Development Goals (SDGs). Here, we develop and validate a global map of the socioeconomic centers of cities for 2020 by integrating nighttime light and population density data within an advanced geospatial modeling framework. Our analysis reveals that monocentric cities -- the standard urban model -- still dominate our planet, accounting for over 80% of cities worldwide. However, these monocentric cities encompass only approximately 20% of the total urbanized area, urban population, and nighttime light intensity; this 80/20 pattern underscores significant disparities in urban development. Further analysis, combined with socioeconomic datasets, reveals a marked difference between developed and developing regions: high-income countries exhibit greater polycentricity than low-income countries, demonstrating a positive correlation between urban sprawl and economic growth. Our global dataset and findings provide critical insights into urban structure and development, with important implications for urban planning, policymaking, and the formulation of indicators for urban sustainability assessment.

## Codes

- Parameters: Constants for convenience.
    * `consts.py`: Model parameters and predefined terms, etc.
    * `plot_utils.py`: Customer plot constants and functions.

- Method: Codes to generate the dataset.
    - Preprocess:
        * `pcca.py`: Module to generate the percolation graphs.
        * `contours_generation.py`: Module to generate the contours map.
    - Main:
        * `urban_areas_generator.py`: Class to generate the urban areas system under the optimal threshold.
        * `centers_generator.py`: Class to identify the socioeconomic centers within the urban areas. 
        * `execute.py`: Module to generate the dataset by calling above classes.

- Graphs: Codes to plot some of the figures in the paper.
    * `average_coverage_area.py`: Figure 2c: *The average coverage area of each center*.
    * `polycentricity.py`: Figure 3b: *Socioeconomic characteristics of monocentric and polycentric cities*.
    * `polycentricity_economics.py`: Figure 4 & Extended Data Figure 5: *Urban center and economic development*. 
    * `pcca_threholds.py`: Extended Data Figure 2: *Parameter estimation in PCCA*.
    * `population_decay.py`: Extended Data Figure 5: *Population density decay curve from urban center to periphery*.
    * `brightest_centroid.py`: Extended Data Figure 6: *Comparison of the identified centers with the brightest grids and centroids*.

### Instruction 
For the purpose of replication, users are recommended to flow this workflow: 
- Run `pcca.py` to generate the percolation graph and determine the optimal threshold for the country.
- Run `contours_generation.py` to generate the contours map of the country.
- Run `execute.py` to generate the cities and centers of the country.

## Dataset
The urban centers are saved in CSV format, with each row representing one center, encoded in UTF-8. The meanings of the fields are given below.
- iso: The ISO alpha-3 code of the country that the center belonged to.
- country_name: The name of the belonged country.
- cluster_id: The id of the belonged cluster, reindexed for each country.
- center_id: The id of the center, reindexed for each urban area.
- is_mc: Whether the center is the main center of the urban area, 0 for false, 1 for true.
- center_name: The geoname of the center.
- latitude: The latitude of the center under the WGS84 coordinate reference system.
- longitude: The longitude of the center under the WGS84 coordinate reference system.

## Contact

If you have any questions, feel free to contact us through email (pwenss2004@gmail.com).
