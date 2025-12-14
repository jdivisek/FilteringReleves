# `resample`: Technical Description of the Function

The `resample` function performs a thinning (resampling) of vegetation-plot data based on a combination of geographic proximity and species composition similarity. The goal is to reduce sampling density in areas where many plots are located close to each other and record very similar (or the same) vegetation.

The plot removal process is iterative: the function first identifies the most similar pair of plots within the specified distance, removes one of them based on a selected rule, and repeats this process until no pairs that violate the specified distance and similarity thresholds remain in the dataset.

The function is built on the highly optimized `data.table` package and enables processing of large datasets (thousands and lower hundreds of thousands of plots) even on standard hardware.

![](images/clipboard-1370335396.png){width="164"}

## Data Requirements

The function requires two input data objects in `data.table` format.

1.  `coord`

    A `data.table` containing coordinates and other metadata for each vegetation plot.

    -   **Required columns:**

        -   `PlotObservationID`: A unique identifier for each plot. It can be in any format (numeric, character).

        -   Second column: X-coordinate or Longitude. Must be of type numeric.

        -   Third column: Y-coordinate or Latitude. Must be of type numeric.

    -   **Rules:**

        -   The coordinate columns must not contain any `NA` values.

        -   Each `PlotObservationID` must be unique.

2.  `spec`

    A `data.table` containing vegetation data in "long format".

    -   **Required columns:**

        -   `PlotObservationID` : The plot identifier, which corresponds to the IDs in `coord`.

        -   `Taxon_name` : The name of the recorded taxon (species).

        -   `cover`: The species cover value (usually in %). Must be of type `numeric`.

    -   **Rules:**

        -   The `cover` column must not contain any `NA` and zero values.

        -   All `PlotObservationID` present in `spec` must also be present in `coord` (and vice versa).

## Function Parameter Descriptions

-   `coord`: A `data.table` containing plot coordinates and metadata. Must meet the requirements listed above.

-   `spec`: A `data.table` containing species data. Must meet the requirements listed above.

-   `longlat`: A logical value (`TRUE`/`FALSE`). If `TRUE`, coordinates are treated as latitude/longitude, and distances are calculated in kilometers. If `FALSE` (default), a projected coordinate system (e.g., UTM) is assumed, and distances are in meters.

-   `dist.threshold`: A numeric value. The threshold for geographic distance. From a pair of plots closer than this value, one will be removed (if they also meet the `sim.threshold`). Units (meters/kilometers) depend on the `longlat` parameter. Default is `1000`.

-   `sim.threshold`: A numeric value (0-1). The threshold for species composition similarity. From a pair of plots more similar than this value, one will be removed (if they also meet the `dist.threshold`). Default is `0.8`.

-   `sim.method`: The method for calculating similarity. Options: `"bray"` (Bray-Curtis), `"simpson"` (Simpson), `"sorensen"` (Sørensen), `"jaccard"` (Jaccard). For all methods except `"bray"`, cover data is automatically converted to presence/absence.

-   `remove`: The rule that decides which plot from a conflicting pair will be removed. Default is `"random"`.

    -   `"random"`: Randomly removes one of the two plots.

    -   `"less diverse"`: Removes the plot with the lower number of species. Ties are broken by `"random"` order.

    -   `"more diverse"`: Removes the plot with the higher number of species. Ties are broken by `"random"` order.

    -   `"lower var.value"`: Removes the plot with the lower value in the column defined by the `var.value` parameter.

    -   `"higher var.value"`: Removes the plot with the higher value in the column defined by the `var.value` parameter.

-   `var.value`: A character string. The name of a column in `coord` used for decision-making with the `"lower var.value"` and `"higher var.value"` methods.

-   `strata`: A character string. The name of a column in `coord` that defines plot stratification. If provided, resampling is performed separately within each stratum (group).

-   `seed`: A number. The seed value for the random number generator, which ensures reproducibility of results for the `"random"` method. Default is `1234`.

## Technical Description of Internal Workflow

The function operates in the following steps:

1.  **Data Preparation and Validation:** Libraries (`data.table`, `spdep`, `Matrix`, `vegan` and `igraph`) are loaded, and a series of checks are performed to verify that the input data meets all requirements regarding format, type, and the absence of NA values.

2.  **Data Preparation and Sorting:** This key step ensures reproducibility and efficiency.

    -   First, the `coord` data.table is sorted according to the rule defined in `remove`. For `"random"`, the rows are randomly shuffled. For other methods, they first randomly shuffled and then sorted by diversity or the value in `var.value`. This serves as a **universal tie-breaking rule** in subsequent steps.

    -   Based on this final order, a new, internal numeric identifier `id` (1, 2, 3...) is created and joined to the `spec` table. All further operations work with this `id`.

3.  **Neighbor Identification:** Using the `spdep::dnearneigh` function, a list of neighbors within the `dist.threshold` is found for each plot. If stratification is used, neighbors are only searched for within the same stratum.

4.  **Splitting into Groups:**

    -   A graph is created from the list of neighbors (using `igraph`).

    -   The `igraph::components` function analyzes this graph and divides all plots into independent, geographically contiguous groups (components).

5.  **The Filtering Core (`filtering_task`):** The actual thinning process is performed for each group individually.

    -   For the given group, similarity between neighboring plots is calculated and a list of all unique pairs of plots that are both geographically close (i.e. neighboring) and exceed the similarity threshold is created.

    -   This list of "conflicting" pairs is sorted in descending order based on their compositional similarity.

    -   The script then iterates through this sorted list, starting with the most similar pair. For the first pair in the list, it decides which plot to remove based on the `remove` rule and adds it to a "blacklist". Subsequently, it removes ALL pairs from the list that contained this just-removed plot. The process is repeated on the reduced list until no conflicting pairs remain.

6.  **Result:** The blacklists of removed plots from all groups are combined. Based on this final blacklist, the original `coord` data.table is filtered, and the thinned and cleaned `coord.filtered` is returned.

## Performance notes

The function was tested with a dataset containing 468,341 vegetation plots from the European Vegetation Archive. On an older PC with 8 GB RAM and Intel Core i5-9400F 2.8 GHz processor, resampling took 20 hours and 53 minutes. However, the actual performance of the function critically depends on the `dist.threshold` value. The larger the value, the bigger the geographically contiguous groups of plots must be processed, and larger distances are also not ecologically very meaningful. It is therefore recommended to set the `dist.threshold` value no higher than 5,000 m (with `dist.threshold = 1000`, the largest group of plots had more than 29,000 plots). Although the `dnearneigh` function, which is used to identify neighboring plots, can handle geographical coordinates in degrees, it is also highly recommended to provide coordinates in a projected coordinate system such as ETRS89. In this case, the function uses Euclidean distance, instead of Great Circle distance, which speeds up identification of neighboring plots.

The function was also tested with a smaller dataset of ca. 114,000 plots from the Czech Vegetation Database. With `dist.threshold = 1000`, this number of plots was processed in 8 minutes on a laptop with 32 GB RAM and Intel Core i7-11850H 2.5 GHz processor. Although faster implementations of this filtering procedure in R are certainly possible, they require calculation of pairwise similarity matrices, which becomes very memory-demanding and inefficient with large datasets.
