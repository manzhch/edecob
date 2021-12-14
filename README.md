# edecob

The *edecob* package can detect sustained change in digital biomarker data. We account for noise using an autoregressive model and use confidence bounds to detect the change.

### Example

```
library(edecob)

# Let us examine the example_data dataset
head(example_data)
#>     subject study_day jump_height detect_lower detect_upper
#> 1 Subject 1         1    55.60844         -Inf     54.41227
#> 2 Subject 1         4    57.77688         -Inf     54.41227
#> 3 Subject 1         7    57.59584         -Inf     54.41227
#> 4 Subject 1        10    59.92832         -Inf     54.41227
#> 5 Subject 1        13    53.33169         -Inf     54.41227
#> 6 Subject 1        16    60.17763         -Inf     54.41227

# We apply the main fuction of the package onto our example_data
example_event <- edecob(example_data, med_win = c(-21,21), bt_tot_rep = 50,
                        min_change_dur = 50)
#> Warning in edecob(example_data, med_win = c(-21, 21), bt_tot_rep = 50,
#> min_change_dur = 50) :
#>   Removing rows where value is NA
names(example_event)
#> [1] "Subject 1"  "Subject 2"  "Subject 3"  "event_info"

# example_event contains the event data for each source
plot(example_event$`Subject 1`)
plot(example_event$`Subject 2`)
plot(example_event$`Subject 3`)

# example_event also contains a data frame containing the event information for all patients
example_event$event_info
#>           event_detected event_onset event_duration event_stop
#> Subject 1           TRUE         169             87       TRUE
#> Subject 2           TRUE         205             51       TRUE
#> Subject 3          FALSE         306             38      FALSE

# Using this data frame, we can draw a survival plot
library("survival")
plot(survfit(Surv(time = event_onset, event = event_detected) ~ 1,
             data = example_event$event_info),
     conf.int = FALSE, xlim = c(0,350), ylim = c(0,1), mark.time = TRUE,
     xlab = "Study Day", ylab = "Survival Probability", main = "Survival plot")
```

Methodology based on Bühlmann, P. (1998). Sieve Bootstrap for Smoothing in Nonstationary Time Series. *The Annals of Statistics*, 26(1), 48-83.
