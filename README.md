# Mapping urban economic centers worldwide
Shuai Pang, Junlong Zhang and Lei Dong

## Abstract

Urban centers are key engines of regional development, yet accurately defining and identifying their economic cores at the global scale remains a major challenge. Existing mapping efforts largely focus on major cities in developed regions and depend on data sources unavailable in many developing countries. Such data limitations hinders the establishment of consistent urban indicators capturing accessibility, urban form, and economic concentration patterns. Here, we develop and validate a global map of city economic centers for the year 2020 by integrating nighttime light and the Global Human Settlement Layer within an advanced geospatial modeling framework. Our analysis reveals that monocentric cities still dominate our planet, accounting for over 80% of all cities. However, these cities encompass only about 36.1% of the total urbanized area, 29.5% of the urban population, and 26% of the nighttime light intensity, revealing large disparities in urban development. By examining the spatial distribution of centers across multiple scales, we find that the average coverage area of economic centers remains remarkably stable across countries. This scaling pattern may imply that city expansion is primarily accompanied by an increase in the number of centers. Consequently, the relationship between urban area and center count yields a near-constant level of average accessibility within cities -- contrasting sharply with monocentric models that predict declining accessibility with increasing city size. Our global dataset and findings offer new insights into the structure and evolution of cities.

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
We provide two types of center: the point location and its contour extent.
### Point
The urban centers are saved in CSV format, with each row representing one center, encoded in UTF-8. The meanings of the fields are given below.
    - iso: The ISO alpha-3 code of the country that the center belonged to.
    - country_name: The name of the belonged country.
    - cluster_id: The id of the belonged cluster, reindexed for each country.
    - center_id: The id of the center, reindexed for each urban area.
    - is_mc: Whether the center is the main center of the urban area, 0 for false, 1 for true.
    - center_name: The geoname of the center.
    - latitude: The latitude of the center under the WGS84 coordinate reference system.
    - longitude: The longitude of the center under the WGS84 coordinate reference system.
### Contour


## Contact

If you have any questions, feel free to contact us through email (pwenss2004@gmail.com).
