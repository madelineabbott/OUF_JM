# Joint longitudinal-survival model for intensive longitudinal data

This repository constains example code for fitting a joint longitudinal-survival model.  This joint model, which summarizes multiple longitudinal outcomes as a smaller number of time-varying latent factors that capture the instantaneous risk of a time-to-event outcome, consists of three submodels: (i) a measurement submodel, (ii) a structural submodel, and (iii) a survival submodel.  We take a Bayesian approach to fitting this model and use Stan.  The portion of code corresponding to longitudinal submodels (i) and (ii) uses code written by Trung Dung Tran; the original code is available on Github: https://github.com/tdt01/LOUmodels.


This repository contains simulated data for demonstration; this data consists of:

* Four measured longitudinal outcomes
* A time-to-event outcome

To fit the model...

For questions, please contact mrabbott@umich.edu.
