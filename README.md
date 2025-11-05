# Mapping urban economic centers worldwide
Shuai Pang, Junlong Zhang and Lei Dong

## Abstract

Urban centers are key engines of regional development, yet accurately defining and identifying their economic cores at the global scale remains a major challenge. Existing mapping efforts largely focus on major cities in developed regions and depend on data sources unavailable in many developing countries. Such data limitations hinders the establishment of consistent urban indicators capturing accessibility, urban form, and economic concentration patterns. Here, we develop and validate a global map of city economic centers for the year 2020 by integrating nighttime light and the Global Human Settlement Layer within an advanced geospatial modeling framework. Our analysis reveals that monocentric cities still dominate our planet, accounting for over 80% of all cities. However, these cities encompass only about 36.1% of the total urbanized area, 29.5% of the urban population, and 26% of the nighttime light intensity, revealing large disparities in urban development. By examining the spatial distribution of centers across multiple scales, we find that the average coverage area of economic centers remains remarkably stable across countries. This scaling pattern may imply that city expansion is primarily accompanied by an increase in the number of centers. Consequently, the relationship between urban area and center count yields a near-constant level of average accessibility within cities -- contrasting sharply with monocentric models that predict declining accessibility with increasing city size. Our global dataset and findings offer new insights into the structure and evolution of cities.

## Codes

- Parameters: Constants for convenience.
    * `consts.py`: Model parameters and predefined terms, etc.
    * `plot_utils.py`: Customer plot constants and functions.

- Method: Codes to generate the dataset.
    - Centers identification:
        * `contours_generator.py`: Class to generate the contours maps.
        * `centers_generator.py`: Class to identify economic centers within urban areas, based on a specified contour parameter value.. 
        * `centers_identification.py`: Module to generate economic centers for each country/region under different contour paramters.
    - Majority voting:
        * `majority_voting.py`: Module to obtain the voted results among different contour paramter values.
    - Attribute enhancement:
        * `attribute_enhancement.py`: Module to provide toponyms and functional categories for centers.
   
- Graphs: Codes to plot some of the figures in the paper.
    * `polycentricity.py`: Figure 3b & 3c: *Characteristics of monocentric and polycentric cities*.
    * `spatial_pattern.py`: Figure 4b, 4c & 4d: *Spatial distribution of centers within the city*. 
    * `average_distance.py`: Figure 5b: *Population-weighted average distance to the center*.

### Instruction 
For the purpose of replication, users are recommended to flow this workflow: 
- Run `centers_identification.py` to generate the results under multiple contour parameters.
- Run `majority_voting.py` to obatin the voted results across different contour paramters.
- Run `attribute_enhancement.py` to enrich the centers with additional attributes.

## Dataset
We provide two types of center: the point location and its contour extent.
### Point
The economic centers are saved in CSV format, with each row representing one center, encoded in UTF-8. The meanings of the fields are given below.
- iso: The ISO alpha-3 code of the country that the center belonged to.
- country_name: The name of the belonged country.
- cluster_id: The id of the belonged GHSL urban area.
- is_mc: Whether the center is the main center of the urban area, 0 for false, 1 for true.
- name: The geoname of the center.
- category: The functional category of the center.
- count: The vote count of the center during the majority voting process.
- is_refined: Whether the center location has been manually refined, 0 for false, 1 for true.
- latitude: The latitude of the center under the WGS84 coordinate reference system.
- longitude: The longitude of the center under the WGS84 coordinate reference system.
### Contour
The contours corresponding to economic centers are saved in CSV format, with each row representing one contour, encoded in UTF-8. The meanings of the fields are given below.
- iso: The ISO alpha-3 code of the country that the center belonged to.
- country_name: The name of the belonged country.
- parameter: The value of the minimum contour area under which the contour is identified.
- area: The enclosed area of the contour, in the unit of square kilometers.
- geometry: The geographic shape of the contour under the WGS84 coordinate reference system.

## Contact

If you have any questions, feel free to contact us through email (pwenss2004@gmail.com).
