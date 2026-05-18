# Urban center data

## Codes

- Method: Codes to generate the dataset.
    - Centers identification:
        * `contours_generator.py`: Class to conduct Gaussian smoothing of NTL rasters and generate contour maps.
        * `centers_generator.py`: Class to identify urban centers within urban areas based on contour topology.
        * `centers_identification.py`: Module to generate urban centers for global urban areas under different contour area thresholds.
    - Attribute enhancement:
        * `attribute_enhancement.py`: Module to assign toponyms (via GeoNames) and functional categories (via Foursquare POI) to centers.

- Graph: Codes to plot the figures in the paper.
    * `scaling.py`: Figure 2: *Linear scaling relationship between city area and center number*.
    * `decomposition.py`: Figure 3: *Consistency between center–population and area–population scaling laws*.
    * `spatial_pattern.py`: Figure 4b, 4c & 4e: *Spatial distribution of centers within the city*.
    * `average_distance.py`: Figure 5b: *Population-weighted average distance to the center*.

### Instruction
For the purpose of replication, users are recommended to follow this workflow:
1. Run `centers_identification.py` to generate center results under multiple contour area thresholds.
2. Run `attribute_enhancement.py` to enrich the centers with toponyms and functional categories.
3. Run the scripts in `Graph/` to reproduce the figures in the paper.

## Dataset
We provide the center point locations and their corresponding contour extents identified under the minimum contour area threshold of 4 km², which is the version used in the main analysis of the paper.

### Centers
The center data is saved in CSV format (`centers.csv`), with each row representing one center, encoded in UTF-8. 

The results of the multimodal large language model (MLLM)-based post-filtering are recorded in the two additional fields (`is_fp` and `fp_type`). The meanings of the fields are given below.
- cluster_id: The ID of the GHSL urban area the center belongs to.
- center_id: The unique identifier for the center.
- latitude: The latitude of the center under the WGS84 coordinate reference system.
- longitude: The longitude of the center under the WGS84 coordinate reference system.
- is_main: Whether the center is the main center of the urban area, 0 for false, 1 for true.
- name: The toponym of the center, assigned from GeoNames.
- category: The functional category of the center, classified from Foursquare POI.
- is_fp: Whether the center is a false positive, 0 for false, 1 for true.
- fp_type: The type of false positive, if applicable.

### Contours
The contours corresponding to urban centers are saved in CSV format (`contours.csv`), with each row representing one contour, encoded in UTF-8. The meanings of the fields are given below.
- cluster_id: The ID of the GHSL urban area.
- center_id: The ID of the corresponding center.
- area: The enclosed area of the contour, in the unit of square kilometers.
- geometry: The geographic shape of the contour under the WGS84 coordinate reference system.

## Citation

Shuai Pang, Junlong Zhang, Yu Liu, and Lei Dong. Global evidence for a consistent spatial footprint of intra-urban centers. https://arxiv.org/abs/2503.06445

## Contact

If you have any questions, feel free to contact us through email (pwenss2004@gmail.com).
