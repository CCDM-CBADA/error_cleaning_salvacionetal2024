library(sf)
library(tidyverse)
library(spdep)

# =============================================================================
# SINGLE PADDOCK YIELD DATA CLEANING PIPELINE
# =============================================================================
# Purpose: Clean and filter yield monitor data for a single paddock
# Author: Kai Bagley, Matthew Nguyen, Arnold Salvacion
# =============================================================================

# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

# Input data (replace with your actual data sources)
paddock_points <- raw[[1]][[1]]  # Point data from yield monitor
paddock_boundary <- bounds_poly[[1]]  # Paddock boundary polygon

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

#' Calculate the mode (most frequent value)
calculate_mode <- function(x) {
  unique_vals <- unique(x)
  unique_vals[which.max(tabulate(match(x, unique_vals)))]
}

#' Identify statistical outliers based on standard deviation
#' @param x Numeric vector
#' @param cutoff Number of standard deviations (default: 2.5)
#' @return Logical vector indicating outliers
identify_outliers <- function(x, cutoff = 2.5) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  (x < mean_val - cutoff * sd_val) | (x > mean_val + cutoff * sd_val)
}

#' Assign paddock field labels to points based on intersection
#' @param points SF point dataframe
#' @param boundary SF polygon with Field column
#' @return SF dataframe with Field labels assigned
assign_field_labels <- function(points, boundary) {
  # Find which boundary polygon each point intersects
  intersection_indices <- st_intersects(points, boundary)
  
  # Convert to numeric indices
  points_with_index <- points %>%
    mutate(index = as.numeric(intersection_indices)) %>%
    select(-any_of("Field"))
  
  # Create lookup table from boundary
  field_lookup <- boundary %>%
    mutate(index = row_number()) %>%
    st_drop_geometry() %>%
    select(index, Field)
  
  # Merge field labels
  merge(points_with_index, field_lookup, by = "index", all = FALSE)
}

# -----------------------------------------------------------------------------
# DATA PREPARATION
# -----------------------------------------------------------------------------

# Remove duplicate geometries
paddock_points <- distinct(paddock_points, geometry, .keep_all = TRUE)

# Assign field labels to points
paddock_data <- assign_field_labels(paddock_points, paddock_boundary)

# Create boundary linestring for edge buffer calculation
boundary_lines <- st_cast(paddock_boundary, "MULTILINESTRING")

# -----------------------------------------------------------------------------
# MAIN CLEANING FUNCTION
# -----------------------------------------------------------------------------

#' Clean yield monitor data for a single paddock
#'
#' @param data SF dataframe with yield point data
#' @param boundary SF boundary (MULTILINESTRING for edge detection)
#' @param yield_col Name of yield column (default: "Yld_Mass_D")
#' @param yield_limits Vector of [min, max] acceptable yield values
#' @param edge_buffer Buffer distance in meters for edge removal (default: 30)
#' @param speed_col Name of speed/velocity column (optional)
#' @param speed_limits Vector of [min, max] acceptable speed values
#' @param swath_col Name of swath width column (optional)
#' @param swath_threshold Proportion of mode for swath filtering (default: 0.7)
#' @param bearing_col Name of bearing/heading column (optional)
#' @param bearing_threshold Minimum proportion for valid bearing bin (default: 0.05)
#' @param distance_col Name of distance column for time lag detection (optional)
#' @param treatment_col Name of treatment column for grouped outlier detection (optional)
#' @param remove_global_outliers Remove statistical yield outliers (default: TRUE)
#' @param remove_spatial_outliers Remove spatial autocorrelation outliers (default: TRUE)
#' @param delete_flagged Remove flagged points instead of marking them (default: FALSE)
#'
#' @return SF dataframe with 'cleaned' column indicating removal reason (NA = kept)

clean_paddock_yield <- function(
    data,
    boundary,
    yield_col = "Yld_Mass_D",
    yield_limits = c(0.1, 20),
    edge_buffer = 30,
    speed_col = NULL,
    speed_limits = c(1, 10),
    swath_col = NULL,
    swath_threshold = 0.7,
    bearing_col = NULL,
    bearing_threshold = 0.05,
    distance_col = NULL,
    treatment_col = NULL,
    remove_global_outliers = TRUE,
    remove_spatial_outliers = TRUE,
    delete_flagged = FALSE
) {
  
  # Initialize cleaning flag column
  data$cleaned <- NA_character_
  
  cat("\n==== Starting paddock cleaning ====\n")
  cat("Initial point count:", nrow(data), "\n\n")
  
  # -------------------------------------------------------------------------
  # 1. YIELD RANGE FILTER
  # -------------------------------------------------------------------------
  cat("Filtering: Yield range... ")
  data$cleaned[data[[yield_col]] < yield_limits[1] | 
                 data[[yield_col]] > yield_limits[2]] <- "yield_range"
  cat(sum(data$cleaned == "yield_range", na.rm = TRUE), "points\n")
  
  # -------------------------------------------------------------------------
  # 2. MISSING YIELD VALUES
  # -------------------------------------------------------------------------
  cat("Filtering: Missing yields... ")
  na_count <- sum(is.na(data[[yield_col]]))
  data <- data[!is.na(data[[yield_col]]), ]
  cat(na_count, "points\n")
  
  if (nrow(data) == 0) {
    cat("\nWARNING: All data removed. Returning empty dataset.\n")
    return(data)
  }
  
  # -------------------------------------------------------------------------
  # 3. EDGE EFFECTS
  # -------------------------------------------------------------------------
  if (!is.null(edge_buffer) && edge_buffer > 0) {
    cat("Filtering: Edge effects (", edge_buffer, "m buffer)... ", sep = "")
    
    # Create buffer around boundary
    buffer_poly <- st_buffer(boundary, edge_buffer)
    
    # Find points within buffer zone
    edge_indices <- st_intersects(buffer_poly, data) %>% unlist()
    
    data$cleaned[edge_indices[is.na(data$cleaned[edge_indices])]] <- "edge"
    cat(sum(data$cleaned == "edge", na.rm = TRUE), "points\n")
    
    if (all(!is.na(data$cleaned))) {
      cat("\nAll data flagged. Stopping cleaning.\n")
      return(data)
    }
  }
  
  # -------------------------------------------------------------------------
  # 4. SPEED/VELOCITY OUTLIERS
  # -------------------------------------------------------------------------
  if (!is.null(speed_col)) {
    cat("Filtering: Speed outliers... ")
    
    data <- data %>%
      mutate(cleaned = if_else(
        is.na(cleaned) & 
          (!!sym(speed_col) < speed_limits[1] | 
             !!sym(speed_col) > speed_limits[2]),
        "speed",
        cleaned
      ))
    
    cat(sum(data$cleaned == "speed", na.rm = TRUE), "points\n")
    
    if (all(!is.na(data$cleaned))) {
      cat("\nAll data flagged. Stopping cleaning.\n")
      return(data)
    }
  }
  
  # -------------------------------------------------------------------------
  # 5. SWATH WIDTH OUTLIERS
  # -------------------------------------------------------------------------
  if (!is.null(swath_col)) {
    cat("Filtering: Swath width outliers... ")
    
    mode_swath <- calculate_mode(data[[swath_col]])
    min_acceptable <- mode_swath * swath_threshold
    
    data <- data %>%
      mutate(cleaned = if_else(
        is.na(cleaned) & !!sym(swath_col) < min_acceptable,
        "swath",
        cleaned
      ))
    
    cat(sum(data$cleaned == "swath", na.rm = TRUE), "points\n")
    
    if (all(!is.na(data$cleaned))) {
      cat("\nAll data flagged. Stopping cleaning.\n")
      return(data)
    }
  }
  
  # -------------------------------------------------------------------------
  # 6. GLOBAL YIELD OUTLIERS (Statistical)
  # -------------------------------------------------------------------------
  if (remove_global_outliers) {
    cat("Filtering: Global yield outliers... ")
    
    if (!is.null(treatment_col)) {
      # Group by treatment for outlier detection
      data <- data %>%
        group_by(!!sym(treatment_col)) %>%
        mutate(cleaned = if_else(
          is.na(cleaned) & identify_outliers(!!sym(yield_col)),
          "yield_outlier",
          cleaned
        )) %>%
        ungroup()
    } else {
      # Whole-paddock outlier detection
      data <- data %>%
        mutate(cleaned = if_else(
          is.na(cleaned) & identify_outliers(!!sym(yield_col)),
          "yield_outlier",
          cleaned
        ))
    }
    
    cat(sum(data$cleaned == "yield_outlier", na.rm = TRUE), "points\n")
    
    if (all(!is.na(data$cleaned))) {
      cat("\nAll data flagged. Stopping cleaning.\n")
      return(data)
    }
  }
  
  # -------------------------------------------------------------------------
  # 7. SPATIAL OUTLIERS (Local Moran's I)
  # -------------------------------------------------------------------------
  if (remove_spatial_outliers) {
    cat("Filtering: Spatial outliers (Local Moran's I)... ")
    
    # Only analyze points not already flagged
    unflagged_data <- data[is.na(data$cleaned), ]
    
    if (nrow(unflagged_data) > 10) {  # Need sufficient points for spatial analysis
      
      # Create distance-based neighborhood (40m threshold)
      neighbors <- st_coordinates(unflagged_data) %>%
        dnearneigh(0, 40)
      
      # Check for points with no neighbors
      neighbor_counts <- card(neighbors)
      isolated_points <- which(neighbor_counts == 0)
      
      if (length(isolated_points) > 0) {
        # Remove isolated points for Moran's I calculation
        analysis_data <- unflagged_data[-isolated_points, ]
        neighbors <- st_coordinates(analysis_data) %>%
          dnearneigh(0, 40)
      } else {
        analysis_data <- unflagged_data
      }
      
      # Calculate Local Moran's I
      moran_results <- localmoran(
        analysis_data[[yield_col]],
        nb2listw(neighbors, style = "W")
      )
      
      # Classify spatial patterns
      analysis_data <- analysis_data %>%
        mutate(
          yield_centered = !!sym(yield_col) - mean(!!sym(yield_col)),
          moran_i = moran_results[, 1],
          moran_p = p.adjust(moran_results[, 5], method = "bonferroni"),
          
          # Classify spatial autocorrelation pattern
          spatial_pattern = case_when(
            moran_p < 0.1 & yield_centered > 0 & moran_i > 0 ~ "High-High",
            moran_p < 0.1 & yield_centered < 0 & moran_i < 0 ~ "Low-High",
            moran_p < 0.1 & yield_centered < 0 & moran_i > 0 ~ "Low-Low",
            moran_p < 0.1 & yield_centered > 0 & moran_i < 0 ~ "High-Low",
            moran_p >= 0.1 ~ "Not Significant"
          ),
          
          # Flag spatial outliers
          cleaned = case_when(
            spatial_pattern == "Low-High" ~ "spatial_low_high",
            spatial_pattern == "High-Low" ~ "spatial_high_low",
            TRUE ~ cleaned
          )
        )
      
      # Map results back to original data
      match_indices <- st_intersects(data, analysis_data) %>% as.numeric()
      data$cleaned[which(!is.na(match_indices))] <- 
        analysis_data$cleaned[match_indices[!is.na(match_indices)]]
      
      cat(sum(data$cleaned %in% c("spatial_low_high", "spatial_high_low"), 
              na.rm = TRUE), "points\n")
    } else {
      cat("Skipped (insufficient unflagged points)\n")
    }
    
    if (all(!is.na(data$cleaned))) {
      cat("\nAll data flagged. Stopping cleaning.\n")
      return(data)
    }
  }
  
  # -------------------------------------------------------------------------
  # 8. BEARING/HEADING OUTLIERS
  # -------------------------------------------------------------------------
  if (!is.null(bearing_col)) {
    cat("Filtering: Bearing outliers... ")
    
    data <- data %>%
      mutate(
        # Normalize bearing to 0-360
        !!sym(bearing_col) := as.numeric(!!sym(bearing_col)) %% 360,
        
        # Round to 10-degree bins
        bearing_bin = round(!!sym(bearing_col), -1)
      ) %>%
      group_by(bearing_bin) %>%
      mutate(
        # Calculate proportion of unflagged points in this bin
        bin_proportion = n() / 
          sum(is.na(data$cleaned) & !is.na(data[[bearing_col]]))
      ) %>%
      ungroup() %>%
      mutate(
        # Flag bearings that appear too infrequently
        cleaned = if_else(
          is.na(cleaned) & bin_proportion <= bearing_threshold,
          "bearing",
          cleaned
        )
      ) %>%
      select(-bearing_bin, -bin_proportion)
    
    cat(sum(data$cleaned == "bearing", na.rm = TRUE), "points\n")
    
    if (all(!is.na(data$cleaned))) {
      cat("\nAll data flagged. Stopping cleaning.\n")
      return(data)
    }
  }
  
  # -------------------------------------------------------------------------
  # 9. TIME LAG / DISTANCE OUTLIERS
  # -------------------------------------------------------------------------
  if (!is.null(distance_col)) {
    cat("Filtering: Time lag outliers... ")
    
    data <- data %>%
      mutate(cleaned = if_else(
        is.na(cleaned) & identify_outliers(!!sym(distance_col)),
        "time_lag",
        cleaned
      ))
    
    cat(sum(data$cleaned == "time_lag", na.rm = TRUE), "points\n")
  }
  
  # -------------------------------------------------------------------------
  # SUMMARY
  # -------------------------------------------------------------------------
  cat("\n==== Cleaning Summary ====\n")
  cat("Points retained:", sum(is.na(data$cleaned)), "\n")
  cat("Points flagged:", sum(!is.na(data$cleaned)), "\n")
  
  if (sum(!is.na(data$cleaned)) > 0) {
    cat("\nFlagging reasons:\n")
    print(table(data$cleaned, useNA = "no"))
  }
  
  # -------------------------------------------------------------------------
  # OPTIONAL: DELETE FLAGGED POINTS
  # -------------------------------------------------------------------------
  if (delete_flagged) {
    cat("\nDeleting flagged points...\n")
    data <- data %>%
      filter(is.na(cleaned)) %>%
      select(-cleaned)
  }
  
  cat("\n==== Cleaning complete ====\n\n")
  return(data)
}

# -----------------------------------------------------------------------------
# EXECUTE CLEANING
# -----------------------------------------------------------------------------

cleaned_paddock <- clean_paddock_yield(
  data = paddock_data,
  boundary = boundary_lines,
  yield_col = "Yld_Mass_D",
  yield_limits = c(0.1, 20),
  edge_buffer = 30,
  remove_global_outliers = TRUE,
  remove_spatial_outliers = TRUE,
  delete_flagged = FALSE  # Set to TRUE to remove flagged points
)

# -----------------------------------------------------------------------------
# OPTIONAL: SAVE OR VISUALIZE RESULTS
# -----------------------------------------------------------------------------

# View cleaning results
# table(cleaned_paddock$cleaned, useNA = "ifany")

# Plot results (requires ggplot2)
# ggplot(cleaned_paddock) +
#     geom_sf(aes(color = is.na(cleaned)), size = 0.5) +
#     scale_color_manual(
#         values = c("TRUE" = "green", "FALSE" = "red"),
#         labels = c("Kept", "Flagged"),
#         name = "Status"
#     ) +
#     theme_minimal() +
#     labs(title = "Paddock Yield Cleaning Results")