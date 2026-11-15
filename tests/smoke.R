library(SurvX)
set.seed(1)
n <- 300
x <- rbinom(n, 1, 0.5); z <- runif(n, -1, 1)
t1 <- rexp(n, 0.5 * exp(0.5 * x)); t2 <- rexp(n, 0.3 * exp(0.3 * z)); cens <- runif(n, 0, 4)
d <- data.frame(time = pmin(t1, t2, cens), x = x, z = z)
d$status <- ifelse(d$time == cens, 0L, ifelse(t1 < t2, 1L, 2L))
d$dead <- as.integer(d$status > 0)

k <- km(Event(time, dead) ~ x, d); print(k)
a <- aj(Event(time, status) ~ x, d); print(a)
f1 <- cox(Event(time, dead) ~ x + z, d); print(f1)
f2 <- fg(Event(time, status) ~ x + z, d); print(f2)
f3 <- csc(Event(time, status) ~ x + z, d); print(f3)
f4 <- directcif(Event(time, status) ~ x + z, d); print(f4)
nd <- data.frame(x = c(0, 1), z = 0)
p2 <- predict(f2, nd, times = c(1, 2)); p3 <- predict(f3, nd, times = c(1, 2)); p4 <- predict(f4, nd, times = c(1, 2))
stopifnot(all(p2 >= 0 & p2 <= 1), all(p3 >= 0 & p3 <= 1), all(p4 >= 0 & p4 <= 1), is.finite(AIC(f4)), length(coef(f4)) == 8L)
stopifnot(abs(p2 - p4) < 0.25)

# user-facing methods
for (obj in list(k, a, f1, f2, f3, f4)) { capture.output(print(obj)); capture.output(print(summary(obj))) }
stopifnot(is.data.frame(estimates(f1)), is.data.frame(estimates(f3)), is.data.frame(estimates(f4)),
          nrow(confint(f4)) == length(coef(f4)), is.list(confint(f3)))
f5 <- update(f4, link = "PO"); stopifnot(f5@link[1] == "PO", is.finite(AIC(f4, f5)$AIC[2]))
f6 <- directcif(Event(time, status) ~ x + z, d, estimation = "weighted"); capture.output(print(summary(f6)))
f7 <- directcif(Event(time, status) ~ x + z, d, estimation = "constrained", cond = "weibull"); capture.output(print(summary(f7)))
cm <- compare_models(Event(time, status) ~ x + z, d, times = c(1, 2), folds = 2); capture.output(print(cm))
grDevices::pdf(NULL)
plot(k); plot(a); plot(f1, newdata = nd); plot(f2, newdata = nd, B = 20); plot(f3, B = 20); plot(f4, newdata = nd); plot(f6); plot(f7, cause = 2); plot(cm)
cif_plot(list(direct = f4, "Fine-Gray" = f2, CSC = f3), B = 20, layout = "overlay"); cif_plot(f4, layout = "stacked"); cif_plot(f4, type = "unadjusted", reference = TRUE)
stopifnot(is.data.frame(cif_data(f4, type = "marginal")))
invisible(grDevices::dev.off())
