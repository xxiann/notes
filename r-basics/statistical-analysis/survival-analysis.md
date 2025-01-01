# survival analysis

[https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/](https://rviews.rstudio.com/2017/09/25/survival-analysis-with-r/)

[https://www.emilyzabor.com/tutorials/survival\_analysis\_in\_r\_tutorial.html](https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html)

[https://bookdown.org/mpfoley1973/survival/](https://bookdown.org/mpfoley1973/survival/)

video: [https://youtube.com/playlist?list=PLqzoL9-eJTNDdnKvep\_YHIwk2AMqHhuJ0\&si=rbeRM6CgH5fSarpx](https://youtube.com/playlist?list=PLqzoL9-eJTNDdnKvep_YHIwk2AMqHhuJ0\&si=rbeRM6CgH5fSarpx)

* HR = 1: No effect
* HR < 1: Reduction in the hazard
* HR > 1: Increase in Hazard
* censored data point refer to subjects who either die of causes other than the disease of interest or are lost to follow-up/ data that we do not know the exact time of event (eg death)
  * are assumed to be non informative for the below 2 models
  * but not true when for eg some subjects that are getting better/worse stopped going back for observations

#### **1. Kaplan Meier Analysis**

* non parametric - no equation and no HR

#### **2. Cox Proportional Hazards Model**

* [http://www.sthda.com/english/wiki/cox-proportional-hazards-model](http://www.sthda.com/english/wiki/cox-proportional-hazards-model)
* semi parametric - can calculate HR but no equation
* suitable for describing the effect of multiple covariates on hazard ratio, and not for predicting

Overall test of the model

* concordance - if close to 50% == more likely predicted at random
* wald, likelihood and log rank tests whether one of the covariate of the model is significant

Fit of model with multiple covariates

* p value of the coefficient
* difference between 2 models (using likelihood ratio test)
* significant covariates but not significant different model can still be added to increase predictive abilities of the model
* log-log plot - testing for proportional hazard assumption
* `cox.zph` - Test the proportional hazards assumption for a Cox regression model fit (goodness of fit)
* whether the factor is time dependent
