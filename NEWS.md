# ph2rand 0.1

* Support established for a variety of single-stage and two-stage randomized
comparative designs, assuming the primary outcome variable is Bernoulli
distributed (`des_one_stage()` and `des_two_stage()`). In addition, functions
created to return the probability mass functions (`pmf()`), terminal points
(`terminal()`), and operating characteristics (`opchar()` and `sim()`) of such
designs. A selection of S3 method functions are also available
(`plot.ph2rand_des()`, `plot.ph2rand_pmf()`, `plot.ph2rand_terminal()`, and
`summary.ph2rand_des()`).