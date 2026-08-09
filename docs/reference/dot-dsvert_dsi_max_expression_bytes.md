# Portable upper bound for one deparsed DSI expression (internal)

DSI 1.8 with DSLite 1.4.1 accepts a 640 KiB Base64url text payload
inside a complete expression measured below the portable 768 KiB
ceiling, but rejects the larger candidates in the recorded data-free
sweep. The default stays strictly below that observed boundary. A
deployment raises it only after every target accepts the stateless
transport probe; the complete expression geometry is then frozen before
its first payload transmission.

## Usage

``` r
.dsvert_dsi_max_expression_bytes()
```
