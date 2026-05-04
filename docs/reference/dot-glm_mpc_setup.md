# GLM Setup: Transport Keys + Standardization

Initializes transport keys on all servers and standardizes features.
Pure Ring63 MPC.

## Usage

``` r
.glm_mpc_setup(
  datasources,
  server_names,
  server_list,
  non_label_servers,
  y_server,
  y_var,
  x_vars,
  data_name,
  family,
  session_id,
  verbose,
  standardize_y_override = NULL,
  std_mode = "full"
)
```

## Value

List with transport_pks, x_means, x_sds, y_mean, y_sd, std_data,
standardize_y, .dsAgg, .sendBlob
