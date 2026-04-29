# Yield Monitor Data Cleaning Pipeline

A comprehensive R script for cleaning and filtering yield monitor data from agricultural machinery, designed for precision agriculture analysis.

## Overview

This pipeline performs multi-stage quality control on spatial yield monitor data to remove erroneous readings caused by:
- Equipment calibration issues
- Edge effects at field boundaries
- Operator behavior anomalies
- GPS positioning errors
- Sensor malfunctions

The cleaned data is suitable for yield mapping, spatial analysis, and precision agriculture decision-making.

## Features

### Core Cleaning Steps

1. **Yield Range Filtering** - Removes physically impossible yield values
2. **Missing Data Removal** - Eliminates points with no yield readings
3. **Edge Effects Removal** - Filters points near field boundaries (buffer zone)
4. **Speed/Velocity Filtering** - Removes points where equipment speed was abnormal
5. **Swath Width Validation** - Identifies incomplete harvester passes
6. **Global Statistical Outliers** - Detects yield values beyond acceptable standard deviations
7. **Spatial Outliers (Local Moran's I)** - Identifies spatially inconsistent readings
8. **Bearing/Heading Validation** - Removes points with unusual travel directions
9. **Time Lag Detection** - Flags points with abnormal distance measurements

### Advanced Features

- **Treatment-aware outlier detection** - Optionally groups data by treatment zones
- **Flexible flagging system** - Points are marked rather than deleted (configurable)
- **Comprehensive logging** - Detailed console output tracking cleaning progress
- **Safety checks** - Prevents over-aggressive filtering with validation stops

## Requirements

### R Packages

```r
library(sf)        # Spatial data handling
library(tidyverse) # Data manipulation
library(spdep)     # Spatial autocorrelation analysis
```

### Installation

```r
install.packages(c("sf", "tidyverse", "spdep"))
```

## Data Requirements

### Input Data Structure

**Yield Points** (`paddock_points`)
- Must be an `sf` spatial dataframe (point geometries)
- Required columns:
  - `Yld_Mass_D` (or custom yield column name) - Numeric yield values
  - Geometry column (spatial coordinates)
- Optional columns:
  - Speed/velocity column
  - Swath width column
  - Bearing/heading column
  - Distance column
  - Treatment zone identifier

**Boundary Polygon** (`paddock_boundary`)
- Must be an `sf` spatial dataframe (polygon geometry)
- Required columns:
  - `Field` - Field identifier for labeling points
  - Geometry column (polygon boundaries)

### Example Data Format

```r
# Yield points structure
paddock_points
#> Simple feature collection with 15000 features and 5 fields
#> Geometry type: POINT
#> Fields: Yld_Mass_D, Speed, Swath, Heading, Distance

# Boundary structure
paddock_boundary
#> Simple feature collection with 1 feature and 1 field
#> Geometry type: POLYGON
#> Fields: Field
```

## Usage

### Basic Usage

```r
# Source the script
source("Error_cleaning_code.R")

# Load your data
paddock_points <- your_yield_data
paddock_boundary <- your_boundary_polygon

# Run cleaning with default parameters
cleaned_data <- clean_paddock_yield(
  data = paddock_data,
  boundary = boundary_lines,
  yield_col = "Yld_Mass_D",
  yield_limits = c(0.1, 20),
  edge_buffer = 30
)
```

### Advanced Configuration

```r
cleaned_data <- clean_paddock_yield(
  data = paddock_data,
  boundary = boundary_lines,
  
  # Yield settings
  yield_col = "Yld_Mass_D",
  yield_limits = c(0.1, 20),           # Min/max acceptable yield (t/ha)
  
  # Spatial settings
  edge_buffer = 30,                     # Edge removal buffer (meters)
  
  # Equipment settings
  speed_col = "Speed",
  speed_limits = c(1, 10),              # Min/max speed (km/h)
  swath_col = "Swath",
  swath_threshold = 0.7,                # 70% of mode swath width
  
  # Direction settings
  bearing_col = "Heading",
  bearing_threshold = 0.05,             # 5% minimum frequency
  
  # Additional columns
  distance_col = "Distance",
  treatment_col = "Treatment",          # For experimental plots
  
  # Outlier detection toggles
  remove_global_outliers = TRUE,
  remove_spatial_outliers = TRUE,
  
  # Output behavior
  delete_flagged = FALSE                # Keep flagged points for review
)
```

## Parameters Reference

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | sf dataframe | Required | Spatial point data with yield readings |
| `boundary` | sf MULTILINESTRING | Required | Field boundary for edge detection |
| `yield_col` | character | `"Yld_Mass_D"` | Name of yield column |
| `yield_limits` | numeric vector | `c(0.1, 20)` | Min/max acceptable yield values |
| `edge_buffer` | numeric | `30` | Edge removal distance in meters |
| `speed_col` | character | `NULL` | Name of speed column (optional) |
| `speed_limits` | numeric vector | `c(1, 10)` | Min/max acceptable speed |
| `swath_col` | character | `NULL` | Name of swath width column |
| `swath_threshold` | numeric | `0.7` | Proportion of mode swath (0-1) |
| `bearing_col` | character | `NULL` | Name of bearing/heading column |
| `bearing_threshold` | numeric | `0.05` | Minimum proportion for valid bearing |
| `distance_col` | character | `NULL` | Name of distance column |
| `treatment_col` | character | `NULL` | Treatment zone grouping variable |
| `remove_global_outliers` | logical | `TRUE` | Enable statistical outlier removal |
| `remove_spatial_outliers` | logical | `TRUE` | Enable Local Moran's I filtering |
| `delete_flagged` | logical | `FALSE` | Remove flagged points vs. mark them |

## Output

### Returned Data Structure

The function returns the input dataframe with an additional `cleaned` column:

- `NA` - Point retained (clean data)
- `"yield_range"` - Yield outside acceptable limits
- `"edge"` - Within edge buffer zone
- `"speed"` - Speed outside acceptable range
- `"swath"` - Swath width too narrow
- `"yield_outlier"` - Statistical outlier (±2.5 SD)
- `"spatial_low_high"` - Low yield surrounded by high (spatial outlier)
- `"spatial_high_low"` - High yield surrounded by low (spatial outlier)
- `"bearing"` - Unusual travel direction
- `"time_lag"` - Abnormal distance measurement

### Console Output Example

```
==== Starting paddock cleaning ====
Initial point count: 15234 

Filtering: Yield range... 342 points
Filtering: Missing yields... 28 points
Filtering: Edge effects (30m buffer)... 1523 points
Filtering: Speed outliers... 87 points
Filtering: Swath width outliers... 156 points
Filtering: Global yield outliers... 234 points
Filtering: Spatial outliers (Local Moran's I)... 89 points

==== Cleaning Summary ====
Points retained: 12775
Points flagged: 2459

Flagging reasons:
           edge   yield_range         speed         swath  yield_outlier 
           1523           342            87           156            234 
spatial_low_high spatial_high_low 
             52            37 

==== Cleaning complete ====
```

## Workflow Integration

### Typical Workflow

```r
# 1. Load raw data
raw_data <- st_read("yield_monitor_data.shp")
boundary <- st_read("field_boundary.shp")

# 2. Run cleaning
cleaned <- clean_paddock_yield(
  data = raw_data,
  boundary = st_cast(boundary, "MULTILINESTRING"),
  delete_flagged = FALSE
)

# 3. Review flagged points
table(cleaned$cleaned, useNA = "ifany")

# 4. Visualize results
library(ggplot2)
ggplot(cleaned) +
  geom_sf(aes(color = is.na(cleaned)), size = 0.5) +
  scale_color_manual(
    values = c("TRUE" = "forestgreen", "FALSE" = "red"),
    labels = c("Retained", "Flagged")
  ) +
  theme_minimal() +
  labs(title = "Yield Data Quality Control Results")

# 5. Export clean data only
clean_data <- cleaned %>% filter(is.na(cleaned))
st_write(clean_data, "cleaned_yield_data.shp")
```

## Helper Functions

### `calculate_mode(x)`
Finds the most frequent value in a vector (used for swath width determination).

### `identify_outliers(x, cutoff = 2.5)`
Identifies statistical outliers beyond a specified number of standard deviations.

### `assign_field_labels(points, boundary)`
Spatially joins point data to field boundary labels using geometric intersection.

## Spatial Outlier Detection Details

The script uses **Local Moran's I** to identify spatial autocorrelation patterns:

- **High-High**: High yield surrounded by high yield (cluster)
- **Low-Low**: Low yield surrounded by low yield (cluster)
- **High-Low**: High yield surrounded by low yield (spatial outlier) ⚠️
- **Low-High**: Low yield surrounded by high yield (spatial outlier) ⚠️

Only High-Low and Low-High patterns are flagged as they represent spatially inconsistent readings likely caused by equipment errors.

## Customization Tips

### Adjusting Sensitivity

**More Conservative** (keep more data):
```r
remove_global_outliers = FALSE
remove_spatial_outliers = FALSE
edge_buffer = 15
```

**More Aggressive** (stricter cleaning):
```r
yield_limits = c(0.5, 15)
edge_buffer = 50
swath_threshold = 0.8
```

### Working with Multiple Paddocks

```r
# Loop through multiple paddocks
paddock_names <- c("North", "South", "East")
cleaned_list <- list()

for (paddock in paddock_names) {
  points <- raw_data %>% filter(Field == paddock)
  boundary <- boundaries %>% filter(Field == paddock)
  
  cleaned_list[[paddock]] <- clean_paddock_yield(
    data = points,
    boundary = st_cast(boundary, "MULTILINESTRING")
  )
}

# Combine results
all_cleaned <- bind_rows(cleaned_list)
```

## Troubleshooting

### Common Issues

**Error: "All data flagged. Stopping cleaning."**
- Your parameters are too strict
- Try widening `yield_limits` or reducing `edge_buffer`

**Warning: "Skipped (insufficient unflagged points)"**
- Normal when previous filters removed most data
- Spatial analysis requires minimum of 10 points

**Edge filtering removes too much data**
- Reduce `edge_buffer` value
- Check that boundary geometry is correct

**Spatial outlier detection fails**
- Ensure coordinate system is projected (not geographic)
- Check for sufficient point density (>40m average spacing)

### Performance Considerations

- Large datasets (>100,000 points): Spatial outlier detection may be slow
- Consider sub-sampling or disabling: `remove_spatial_outliers = FALSE`
- Boundary geometry complexity affects edge filtering speed

## Citation

If you use this code in your research, please cite accordingly and acknowledge the data cleaning methodology.

## License

[Specify your license here]

## Author

Refactored for clarity and single-paddock use

## Version History

- **Current**: Single-paddock optimized version with comprehensive documentation
- Supports flexible parameter configuration
- Enhanced console logging and safety checks

## Additional Resources

- [sf package documentation](https://r-spatial.github.io/sf/)
- [spdep spatial statistics](https://r-spatial.github.io/spdep/)
- Precision agriculture data cleaning best practices

---

**Note**: Always review flagged points visually before deleting them. Set `delete_flagged = FALSE` initially to preserve data for quality assurance checks.
