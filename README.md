# Joint longitudinal-survival model for multivariate intensive longitudinal data

This repository constains example code for fitting a joint longitudinal-survival model.  This joint model, which summarizes multiple longitudinal outcomes as a smaller number of time-varying latent factors that capture the instantaneous risk of a time-to-event outcome, consists of three submodels: (i) a measurement submodel, (ii) a structural submodel, and (iii) a survival submodel.  We take a Bayesian approach to fitting this model and use Stan.  The portion of code corresponding to longitudinal submodels (i) and (ii) builds on code written by Trung Dung Tran; the original code is available on the author's Github: https://github.com/tdt01/LOUmodels.


This repository contains simulated data for demonstration purposes.  Settings 1 and 2 correspond to two different sets of true parameter values.  These parameters values are given in ``true_values_s1`` and ``true_values_s2.csv``.  For both settings, we consider four different measurement patterns for the longitudinal outcomes.  These measurement patterns vary in both the total number and intervals between the observations of the longitudinal outcomes.

* The longitudinal data consists of measurements of four outcomes and is provided in ``sim_dat/long_dat_s[setting]_meas_pattern_[measurement pattern]_g1.csv``
* The time-to-event outcome is provided in ``sim_dat/surv_dat_s[setting]_meas_pattern_[measurement pattern]_g1.csv``

To fit the model, run ``fit_model/fit_jm.R``.

For questions, please contact mrabbott@umich.edu.
