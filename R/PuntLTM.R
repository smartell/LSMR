# This a simple R-script to test the Punt et al (1997) paper.


# vonB growth parameters
dt		<- 1/4
da		<- seq(0, 30, dt)
linf	<- 150
k		<- 0.3

len 	<- linf*(1-exp(-k*da))
plot(da, len, type="o")


# equivalent parameters for the generalized Schnute model.
c	<- 1
cv	<- 0.2
l1	<- 80
l2	<- 100
lam1	<- l1 + (linf-l1)*(1-exp(-k*dt))
lam2	<- l2 + (linf-l2)*(1-exp(-k*dt))

a	<- (lam2^c-lam1^c)/(l2^c-l1^c)
b	<- (lam1^c*l2^c - l1^c*lam2^c) / (lam1^c-l1^c+l2^c-lam2^c)

# Size intervals
x	<- seq(0, linf/0.85, 5)
xa  <- x[1:(length(x)-1)] + 0.5*diff(x)

r	<- 1/cv^2
si	<- (a*xa^c+b*(1-a))^(1/c)-xa
gamma	<- 1/(si*cv^2)


