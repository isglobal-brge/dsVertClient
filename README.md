# dsVertClient

The `dsVertClient` package complements the `dsVert` package, facilitating the client-side operations for analyzing vertically partitioned data within the DataSHIELD environment. This package enables users to interact securely with data stored across multiple sources without compromising privacy.

## Features

- **Client-Side Interface**: Provides the necessary tools for clients to interact with vertically partitioned datasets.
- **Privacy Protection**: Ensures that the privacy of the data is maintained, following the strict protocols of DataSHIELD.

## Installation

Install `dsVertClient` from GitHub using R:

```R
devtools::install_github("isglobal-brge/dsVertClient")
```

## Usage

To use dsVertClient, load it into your R session:

```R
library(dsVertClient)
```
For more detailed instructions on how to use the package, refer to the vignettes:

```R
browseVignettes(package = "dsVertClient")
```

## License

dsVertClient is licensed under the MIT License. See the LICENSE file for more details.
