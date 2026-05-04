# K\>=3 Ring63 DCF + Beaver Gradient Loop

Pure Ring63 protocol for K\>=3 binomial/Poisson GLM. Uses additive
secret sharing of features to 2 DCF parties, DCF wide spline for
sigmoid/exp, and Beaver matvec for gradient. Pure Ring63 additive secret
sharing.

## Details

Protocol:

1.  One-time: all servers share features+y with 2 DCF parties (Ring63
    FP)

2.  Per iteration: DCF parties compute eta from X_share\*beta

3.  DCF wide spline between DCF parties -\> mu shares

4.  Beaver matvec between DCF parties -\> gradient shares

5.  Client aggregates Ring63 shares -\> plaintext gradient

6.  Client L-BFGS update

Requires: transport keys only. No CPK, Galois, or RLK. Precision: ~1e-3
(DCF approximation).
