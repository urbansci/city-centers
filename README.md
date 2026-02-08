Intra-urban center data
Shuai Pang, Junlong Zhang, Yu Liu and Lei Dong

## Abstract

Urban space is highly heterogeneous, with population and human activities concentrating in localized centers. However, the global organization of such intra-urban centers remains poorly understood due to the lack of consistent, comparable data. Here we develop a scalable geospatial framework to identify intra-urban activity centers worldwide using nighttime light observations. Applying this approach to more than 9,500 cities, we construct a high-resolution global dataset of over 15,000 centers. We uncover a striking regularity: despite vast differences in city size, regional development, and population density, the built-up area associated with individual centers remains remarkably consistent. Across cities, total urban area scales proportionally with the number of centers, yielding a stable mean spatial footprint. This regularity holds at the micro-scale, where Voronoi-based service areas exhibit a characteristic size that is persistent across countries and independent of local population concentration. As a geometric consequence, this polycentric multiplication maintains stable average distances to the nearest center as cities expand, preventing the accessibility decay inherent in monocentric growth. These findings reveal a universal organizing principle whereby urban expansion is accommodated through the replication of activity centers with a consistent spatial extent, providing a new empirical foundation for understanding the nature of urban growth.

## Codes

- Parameters: Constants for convenience.
    * `consts.py`: Model parameters and predefined terms, etc.
    * `plot_utils.py`: Customer plot constants and functions.

- Method: Codes to generate the dataset.
    - Centers identification:
        * `contours_generator.py`: Class to generate the contours maps.
        * `centers_generator.py`: Class to identify urban centers within urban areas, based on a specified contour parameter value.. 
        * `centers_identification.py`: Module to generate urban centers for each country/region under different contour paramters.
    - Majority voting:
        * `majority_voting.py`: Module to obtain the voted results among different contour paramter values.
    - Attribute enhancement:
        * `attribute_enhancement.py`: Module to provide toponyms and functional categories for centers.
   
- Graphs: Codes to plot some of the figures in the paper.
    * `area_vs_center_number.py`: Figure 2b & 2c: *Geographic distribution of global urban centers*.
    * `spatial_pattern.py`: Figure 3b, 3c, 3d & 3e: *Spatial distribution of centers within the city*. 
    * `average_distance.py`: Figure 4b: *Population-weighted average distance to the center*.
    * `population_decay.py`: Figure 5: *Population density decay from urban center to peripheries*.

### Instruction 
For the purpose of replication, users are recommended to flow this workflow: 
- Run `centers_identification.py` to generate the results under multiple contour parameters.
- Run `majority_voting.py` to obatin the voted results across different contour paramters.
- Run `attribute_enhancement.py` to enrich the centers with additional attributes.

## Dataset
We provide two types of center: the point location and its contour extent.
### Point
The point data is available in two versions, distinguished by whether they have been manually refined.
   * `centers_raw.csv`: The direct output of the detection algorithm without human intervention.
   * `centers_refined.csv`: The version after human main center refinement procedure to enhance usability in applied contexts.

Both are saved in CSV format, with each row representing one center, encoded in UTF-8. The meanings of the fields are given below.
- iso: The ISO alpha-3 code of the country that the center belonged to.
- country_name: The name of the belonged country.
- cluster_id: The id of the belonged GHSL urban area.
- is_mc: Whether the center is the main center of the urban area, 0 for false, 1 for true.
- name: The geoname of the center.
- category: The functional category of the center.
- count: The vote count of the center during the majority voting process.
- is_refined: Whether the center location has been manually refined, 0 for false, 1 for true. Only for `centers_refined.csv`.
- latitude: The latitude of the center under the WGS84 coordinate reference system.
- longitude: The longitude of the center under the WGS84 coordinate reference system.
### Contour
The contours corresponding to urban centers are saved in CSV format, with each row representing one contour, encoded in UTF-8. The meanings of the fields are given below.
- iso: The ISO alpha-3 code of the country that the center belonged to.
- country_name: The name of the belonged country.
- parameter: The value of the minimum contour area under which the contour is identified.
- area: The enclosed area of the contour, in the unit of square kilometers.
- geometry: The geographic shape of the contour under the WGS84 coordinate reference system.

### Citation

- Shuai Pang, Junlong Zhang, Yu Liu, and Lei Dong. Global evidence for a consistent spatial footprint of intra-urban centers. https://arxiv.org/abs/2503.06445

## Contact

If you have any questions, feel free to contact us through email (pwenss2004@gmail.com).
