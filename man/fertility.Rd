\name{fertility}
\alias{fertility}
\docType{data}
\title{Fertility history of rural women in Shiraz: competing risks data}
\description{
Time from marriage to a live birth (the event of interest) or to a stillbirth or abortion (the
competing event), or to the end of follow-up without either (censored), in a cross-sectional
study of 652 women aged 15 to 49 years in rural areas of Shiraz, Iran. The data were analysed
by Shayan, Ayatollahi and Zare (2011) to model the cumulative incidence of live birth with
competing risks.
}
\usage{fertility}
\format{
A data frame with 652 rows and 5 variables:
\describe{
  \item{\code{id}}{subject identifier.}
  \item{\code{time}}{time from marriage to the event or to censoring, in years.}
  \item{\code{status}}{\code{0} = censored (45 women), \code{1} = live birth (574),
    \code{2} = stillbirth or abortion (33).}
  \item{\code{education}}{factor: \code{"illiterate"} (93), \code{"primary/secondary"} (468),
    \code{"higher"} (91).}
  \item{\code{age}}{age of the woman (years) at the time of the event or of censoring, as recorded
    in the source.}
}
}
\details{
Times are recorded to about a day and many are tied (390 of 652 values are repeats), which makes
the data a good test of the handling of tied times. Because \code{age} is measured at the end of
follow-up, it is not a baseline covariate: use it for description, not as a predictor of the
time to the event. The same data are distributed with the \pkg{cmpp} package of the same author
(as \code{fertility_data}, with \code{Event} and \code{Education} coded as numbers).
}
\source{
Shayan Z, Ayatollahi SMT, Zare N (2011). A parametric method for cumulative incidence modeling
with a new four-parameter log-logistic distribution. Theoretical Biology and Medical Modelling
8:43. \doi{10.1186/1742-4682-8-43}
}
\examples{
str(fertility)
table(fertility$status)

## cumulative incidence of live birth and of stillbirth/abortion by education
fit <- aj(Event(time, status) ~ education, fertility)
fit
plot(fit)

## effect of education on the cumulative incidence of live birth
fg(Event(time, status) ~ education, fertility, cause = 1)
}
\keyword{datasets}
