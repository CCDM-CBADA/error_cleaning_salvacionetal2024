library(sf)
library(tidyverse)

## Insert paddock dataframe here
paddock_df <- raw[[1]][[1]]

## Insert paddock boundary here
boundaries <- bounds_poly[[1]]




## Remove geometry duplicates
paddock <- distinct(paddock_df, geometry, .keep_all = TRUE)

faster_intersection <- function(polygon, points) {
    ## Returns a list of indices, index of test_pol it intersects with
    inters <- st_intersects(points, polygon)

    temp_dat <- points %>%
        mutate(index = as.numeric(inters)) %>%
        select(-any_of(c("Field")))

    table <-
        polygon %>%
        mutate(index = row_number()) %>%
        st_drop_geometry() %>%
        select(index, Field)

    ## Now we can Inner Join
    return(merge(temp_dat, table, by = "index", all = FALSE))
}

paddock <- faster_intersection(boundaries, farm_list)

## Allows the buffer to be created from the boundaries
bounds_mls <- st_cast(boundaries, "MULTILINESTRING")


paddock_cleaning <- function(paddock,
                             boundary,
                             response = "Yld_Mass_D",
                             resp_limit = c(0.1, 20),
                             remove_outliers = FALSE,
                             Zero = TRUE,
                             Edge = TRUE,
                             Bearing = NULL,
                             bin_threshold = 0.05,
                             Speed = NULL,
                             vel_upperlimit = 100,
                             Swath_Width = NULL,
                             Distance = NULL,
                             TRT = NULL,
                             Global = TRUE,
                             Local = TRUE) {
    ## Converting Response to Symbol
    ## response_sym <- sym(response)

    mode_avg <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }

    cat("\n**** Marking ")
    outlier <- function(x, cutoff = 2.5) {
        mn <- mean(x)
        st <- stats::sd(x)

        ## TRUE if value is further than 2.5 sd from mean
        cond <- (x < mn - cutoff * st) | (x > mn + cutoff * st)
        return(cond)
    }

    ## Drop yield above threshold
    paddock$cleaned[
                paddock[[response]] < resp_limit[1] |
                paddock[[response]] > resp_limit[2]
            ] <- "yld"

    ## Drop NA yields
    cat("NAs, ")
    paddock <- paddock[!is.na(paddock[[response]]), ]

    ## Remove Edge effects (Requires Boundary)
    if (Edge) {
        cat("edges, ")
        radius <- 30

        ## union them so we can go back to one operation
        buffer <- st_buffer(boundary, radius)

        ## remove buffer points from the trial area (requires s2 to be off)
        ## (this should???? be fine as the area we are working with is close enough that
        ## there are no problems interpreting it to be a flat plane rather than points on
        ## a sphere)

        edge_idx <- st_intersects(buffer, paddock) %>% unlist()

        paddock <- paddock %>%
            mutate(cleaned = replace(
                       cleaned,
                       edge_idx,
                       "edge"
                   ))

        # Exit if all data has been cleaned
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }

    ## Remove Velocity Outliers
    if (!is.null(Speed)) {
        cat("velocities, ")
        paddock <-
            paddock %>%
            mutate(cleaned = if_else(
                       "&"(
                           is.na(cleaned),
                           (!!sym(Speed) < 1 | !!sym(Speed) > vel_upperlimit)
                       ),
                       "vel",
                       cleaned
                   ))
        # Exit if all data has been cleaned
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }

    ## Remove Swath Width outliers
    if (!is.null(Swath_Width)) {
        cat("swath widths, ")
        paddock <- paddock %>%
            mutate(
                cleaned =
                    if_else(
                        "&"(
                            is.na(cleaned),
                            (!!sym(Swath_Width) < mode_avg(!!sym(Swath_Width))*0.7)
                        ),
                        "swth",
                        cleaned
                    )
            )
        # Exit if all data has been cleaned
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }

    ## Global Yield Outliers
    if (Global) {
        cat("yields, ")
        if (!is.null(TRT)) {
            paddock <- paddock %>%
                group_by(TRT) %>%
                mutate(
                    cleaned =
                        if_else(
                            "&"(
                                is.na(cleaned),
                                outlier(!!sym(response))
                            ),
                            "yld",
                            cleaned
                        )
                ) %>%
                mutate(cleaned = if_else(is.na(!!sym(response)), "missing", cleaned)) %>%
                ungroup()
        } else {
            paddock <- paddock %>%
                mutate(
                    cleaned =
                        if_else(
                            "&"(
                                is.na(cleaned),
                                outlier(!!sym(response))
                            ),
                            "yld",
                            cleaned
                        )
                ) %>%
                mutate(cleaned = if_else(is.na(!!sym(response)), "missing", cleaned))
        }
        # Exit if all data has been cleaned
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }

    ## Spatial Outliers (Local)
    ## Only consider non-outliers
    if (Local) {
        cat("spatial yields, ")

        ## Select non-outliers
        moran_padd <-
            paddock[is.na(paddock$cleaned), ]

        ## Distance-based neighbourhood method:
        neighbours <-
            sf::st_coordinates(moran_padd) |>
            spdep::dnearneigh(0, 40)

        c.nb <- spdep::card(neighbours)
        bad <- which(c.nb == 0)

        moran_out <- NULL
        padd_temp <- moran_padd

        if (length(bad) != 0) {
            padd_temp <- padd_temp[-bad, ]

            neighbours <-
                sf::st_coordinates(padd_temp) |>
                spdep::dnearneigh(0, 40)

            moran_out <- localmoran(
                padd_temp[[response]],
                nb2listw(neighbours, style = "W")
            )
        } else {
            moran_out <-
                spdep::localmoran(
                           padd_temp[[response]],
                           spdep::nb2listw(neighbours, style = "W")
                       )
        }

        ## Local Moran
        paddock_spr <-
            padd_temp %>%
            mutate(
                yield_c = !!sym(response) - mean(!!sym(response)),
                locmoran = moran_out,
                morani = locmoran[, 1],
                pr = locmoran[, 5] %>% stats::p.adjust(method = "bonferroni"),
                ## moran_clf = attr(moran_out,"quadr")$pysal,

                moran_clf =
                    case_when(
                        pr < 0.1 & yield_c > 0 & morani > 0 ~ "High-High",
                        pr < 0.1 & yield_c < 0 & morani < 0 ~ "Low-High",
                        pr < 0.1 & yield_c < 0 & morani > 0 ~ "Low-Low",
                        pr < 0.1 & yield_c > 0 & morani < 0 ~ "High-Low",
                        pr >= 0.1 ~ "NS"
                    ),
                cleaned = ifelse(
                    "&"(is.na(cleaned), moran_clf == "Low-High"),
                    "spat_lh",
                    cleaned
                ),
                cleaned = ifelse(
                    "&"(is.na(cleaned), moran_clf == "High-Low"),
                    "spat_hl",
                    cleaned
                )
            ) %>%
select(-locmoran)

        ## Adding the Spatial Outlier label to each point
        idx <- st_intersects(paddock, paddock_spr) %>% as.numeric()

        paddock$cleaned[which(!is.na(idx))] <- paddock_spr$cleaned

        # Exit if all data has been cleaned
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }


    ## Remove "GPS errors" based on heading/bearing of header
    ## if (!is.null(Bearing)) {
    ##     cat("bearing, ")
    ##     paddock <-
    ##         paddock %>%
    ##         mutate(
    ##             cleaned = if_else(
    ##                 "&"(
    ##                     is.na(cleaned),
    ##                     round(!!sym(Bearing), -1) %% 180 != mode_avg(round(!!sym(Bearing), -1) %% 180)
    ##                 ),
    ##                 "gps",
    ##                 cleaned
    ##             )
    ##         )
    ##     # Exit if all data has been cleaned
    ##     if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    ## }


    ## Multimodal alternative
    if (!is.null(Bearing)) {
        cat("bearing, ")
        paddock <-
            paddock %>%
            mutate(!!sym(Bearing) := as.numeric(!!sym(Bearing))) %>%
            mutate(!!sym(Bearing) := !!sym(Bearing) %% 360) %>%
            ## Rounding method
            mutate(BearingBin = round(!!sym(Bearing), -1)) %>%
            ## Binning method (alternative)
            ## mutate(
            ##     BearingBin = cut(!!sym(Bearing), breaks = seq.int(0, 360, by = 10), right = FALSE)
            ## ) %>%
            group_by(BearingBin) %>%
            ## Only want the proportion of non NA Bearing/outlier points
            mutate(
                BinProp = n()/
                    nrow({{.}} %>% filter("&"(!is.na(BearingBin), is.na(cleaned))))
            ) %>%
            ungroup() %>%
            ## Legitimate direction if it comprises at least bin_threshold% of considered bearings
            mutate(
                cleaned = if_else(
                    "&"(
                        is.na(cleaned),
                        BinProp <= bin_threshold
                    ),
                    "gps",
                    cleaned
                )
            ) %>%
            select(-BearingBin, -BinProp)
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }


    ## Remove Time lag Outliers (For each Product)
    if (!is.null(Distance)) {
        cat("time lags, ")
        paddock <- paddock %>%
            mutate(
                cleaned =
                    if_else(
                        "&"(
                            is.na(cleaned),
                            outlier(!!sym(Distance))
                        ),
                        "tlag",
                        cleaned
                    )
            )
        # Exit if all data has been cleaned
        if (sum(is.na(paddock$cleaned)) == 0) return(paddock)
    }


    ## Remove outliers from the data if remove_outlier is true
    if (remove_outliers) {
        cat("deleting, ")
        paddock <- paddock %>%
            filter(!(is.na(cleaned))) %>%
            select(-cleaned)
    }

    cat("done!")
    return(paddock)
}

paddock_cl <-
    paddock_cleaning(
        paddock,
        bounds_mls
    )
