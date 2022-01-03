################################################################################
# Demographic stochasticity in a simple Malthusian model 
# for a homogeneous population (i.e., without age/size/spatial structure).
# by Giulio De Leo. First created in 2017, Updated in Dec. 2020
################################################################################

################################################################################
# Learning goals
################################################################################
# Here we will learn to:
# - simulate a small population whose abundance is described by integer numbers
# - run the code in a smart way by using the binomial probability distribution 
# - replicate simulations one at a time
# - replicate lot of simulations a large number of times
# - speed up the computation by "vectorizing" the for loop for the replicates
# - assess probability of extinction and quasi-extinction
# - assess how these probability changes with time and initial conditions
# - derive the probability distribution of the population at any given time step
# - compute the average time to reach a quasi-extinction threshold  
################################################################################

# Let's start by removing all previous variables from R's memory
rm(list=ls(all=TRUE)) 
################################################################################

# Let's define model parameters for a semelparus, asexual, Malthusian population
fec   <- 2    # this is per-capita fecundity
sigma <- 0.51 # this is probability of surviving to next time step

# compute the finite growth rate:
(lambda<- sigma*fec) 

# note that lambda is slightly larger than 1 so, we expect the population to
# increase

# For instance, we can compute the doubling time, from:
#    N(x) = 2*N(t=0) = N(t=0)*lambda^x
# where "x" is the (unknown) time at which the population doubles in size, so is
# twice as large as the initial population size at time t=0.
# Dividing both sides by N(t=0), we get:
#    2 = lambda^x
# We then take the natural logarithm of each side   
#    log(2) = x * log(lambda)
# and solve for x, namely:
#    x = log(2)/log(lambda)

log(2)/log(lambda)

# Let's now set initial condition and the length of the simulation window
Tmax <- 100         # simulation time
No <- 20            # initial conditions
N <- numeric(Tmax)  # set memory apart for the vector of pop. density
N[1] <- No          # assign the initial condition

# now run the deterministic model...
for(t in (1:(Tmax-1))) { 
   N[t+1] = lambda * N[t] }

# and plot the result
plot(N, type = "l", col = "grey", ylim=c(0,max(N)*1.5), lwd=5)

################################################################################
# Now, let's see how a  stochastic simulation looks like
################################################################################
# For simplicity, we will first simulate the stochastic dynamics of a 
# semelparous population where the individuals reproduces and then die.
# We will track only females in the population, and assume that "fec" is the 
# (exact) number of female offspring/female, so, there is no stochasticity 
# in reproduction (We will later relax this hypothesis).
# Here we will assume that stochasticity is only in the survival process: 
# basically, we will flip a (slightly biased) coin to check whether it will be
# able to survive to the next time step   

Ns = numeric(Tmax)  # initialize the vector for pop. density
Ns[1] = No          # set the same initial conditions

# Let set the seed for the random sequence to a common value, so we can all
# have the same sequence:
set.seed(3)

for (t in (1:(Tmax-1))) {
	NG = fec * Ns[t] 	# compute the number of recruits at this time step
	j = 0				      # initialize the variable
	
	for(k in 1:NG) {  # here we flip a coin for each new recruit in the population
	  
	  if (runif(1)<=sigma) j <- j+1  # in runif() is smaller than sigma, then this 
	                                 # individual survives, otherwise it does not 
	}
	Ns[t+1] <- j  		# this is the overall number of individual that survived
	if (j==0) break
}

# now plots the result, along with that of the deterministic model

points(Ns, type = "l", col = sample(colors(),1))

# Repeat the exercise by running again the for cycle (but do NOT run the
# set.seed(3) line again!) and plotting the results. Because you get a different
# sample of random numbers, you should get each time a  different
# demographic trajectory

################################################################################
# There is a much faster way than flipping a coins for each individual to
# determine whether it dies or survives: we can use the *binomial*
# distribution and the random generation rbinom() function in R.
# For instance, if I want to flip an unbiased coin just one time, 
# this is how I do it:
rbinom(n = 1, size = 1, prob = 0.5)

# If I want to do it 4 times:
rbinom(n = 4, size = 1, prob = 0.5)  

# If I have 10 coins, I flip them all and I want *to count* how many heads, 
# I can do it in this way:
rbinom(n = 1, size = 10, prob = 0.5)  

# this is equivalent to say that I have 10 individuals in the population and
# I want to check how many of them survive, by using their specific 
# survival to the next time step

# so, let's now use "rbinom" to re-run our simulations
Ns    <- Ns*0  # re-set to zero all the elements on Ns
Ns[1] <- No
# plot first the deterministic simulations
plot(N, type = "l", col = "grey", ylim = c(0,max(N, Ns)), lwd=5)

# in the case of our Malthusian model, this is how we change the code:
set.seed(3) 


for(t in (1:(Tmax-1))) {
	NG      <- fec * Ns[t]
	Ns[t+1] <- rbinom(1,NG,sigma)
	print(c(t,NG, Ns[t+1]))
	if (Ns[t+1]==0) {break}
}
points(Ns, type = "l", col = sample(colors(),1))

# so, each time we generate a different realization of the stochastic process.

################################################################################
# now, let's repeat this exercise a number of times, for instance 10
Nrep   <- 10  # number of replicates
Nr     <-   matrix(0, nrow = Tmax, ncol=Nrep) # initialize the matrix
Nr[1,] <- No  # set the the initial conditions 
ng     <- 0   # this is a dummy variable to temporary store the number 
              # of new recruits

for (r in 1:Nrep) {  # two nested for loops
	for (t in (1:(Tmax-1))) {
		ng  <- fec * Nr[t,r]
		Nr[t+1,r] <- rbinom(1,ng,sigma)
	}
}

# print the results, 10 years for the first 5 replicates
Nr[1:10, ]


# and plot these replicates
matplot(Nr, type = "l")
points(N, type = "l", col = "grey", lwd=4)

################################################################################
# Note that in some replicates the population grows, in others it vanishes. 
# We can compute
# the number of times the population goes extinct and divide for the total
# number of times we replicated this exercise. In order to get a robust
# statistics not overly affected by the stochastic nature of the process,
# we have to repeat this exercise a very large number of times, at least
# 1,000 or more...

Nrep   <- 10000 # number of replicates
Nr     <- matrix(0, nrow = Tmax, ncol=Nrep) #  set memory aside
Nr[1,] <- No    #initial conditions

# in this case, let's track precisely the time required to compute 
# the 10,000 simulations of 100 years each. In order to do so, we use the 
# R function system.time. The reasons will be clear in a moment

system.time(
for (r in 1:Nrep) {
	for (t in (1:(Tmax-1))) {
		ng <- fec * Nr[t,r]
		Nr[t+1,r] <- rbinom(1,ng,sigma)
	}
}
)

# and write down the computation time (it was about 17 seconds 
# on my old laptop and now it is only 2.29 !)
# Note also that above we have two nested for-loops: 
#   - the internal one is for the recursive equation to cycle through time
#   - the external one is to replicate the exercise Nrep times

################################################################################
# there is a way to speed up these computations a lot, i.e.
# by "vectorizing" the replicates. First we define the dummy variable: 

# first we reset the matrix to store the results.
Nr  <- matrix(0, nrow = Tmax, ncol=Nrep) 

# and assign again the initial conditions
Nr[1,]  <-No 

# we also create a dummy vector to temporary store the number of recruits for
# each of the NRep replicates at a specific time step
NG <- matrix(0, nrow = 1, ncol=Nrep) 

# and run the model again, but this time with only *one* (1) for-loop!
system.time(
for(t in (1:(Tmax-1))) {
  # these is how simple is to run all the replicates altogether: we select a row
  # of Nr, multiply for fecundity and assign the result to the dummy vector
	NG  <- fec * Nr[t,] 
	
	# then, we extract for each of the Nrep replicate how many individuals survive 
	Nr[t+1,] <- rbinom(NG,NG,sigma) 
}
)
# check the time: the vectorize computation should take substantially less than
# the two nested cycles
################################################################################

################################################################################
# In this exercise we assumed that demographic stochasticity is only about the
# discrete number of individuals that survive to the next generation. In
# reality, the number of individual generated by each parent is also discrete.
# If you have precise information on the distribution probability of number of
# offspring generated by a mother, then it is worth using it. For grey reef
# shark, for instance, the number of offspring ranges between 2 and 6, with
# probability 0.3,0.4,0.2,0.1 of generating 3, 4, 5 and 6 baby shark
# respectively. For Leopard shark, the mean is about 20 +/- 7 (SD). 
# Anyway, when information on the exact range and probability distribution is
# lacking, ecologists often use the discrete Poisson (as a side note, for the
# majority of marine organisms that are broadcast spawner and produce a large
# number of eggs at each reproductive event, environmental stochasticity is
# probably more relevant, or the number of individuals that are actually
# reproducing). We can simulate this by using, for instance a Poisson
# distribution "rpois(n,lambda)" with lambda=fec and n = Nr, i.e., the number of
# individual in the population at a given time t.  Accordingly, we can change
# the r line:
#   NG <-  fec * Nr[t,]
#
# with: 
#   NG  <-  sum(rpois(Nr[t,],fec))
#
# Note that the larger Nr[t,], the better is the mean approximation NG = fec *
# Nr[t,] As usual, it is a small population size that demographic stochasticity
# plays a very relevant role.


# and assign again the initial conditions
Nr[1,]  <-No 

# and run the model again, but this time with only *one* "for-loop"!
system.time(
  for(t in (1:(Tmax-1))) {
    # these is how simple is to run all the replicates altogether: we select a
    # row of Nr, multiply for fecundity and assign the result to the dummy
    # vector
    NG  <- sapply(Nr[t,], function(z) {sum(rpois(z,fec))}) 
    
    # then, we extract for each of the Nrep replicate how many individuals
    # survive
    Nr[t+1,] <- rbinom(NG,NG,sigma) 
  }
)

matplot(Nr, type = "l")
points(N, type = "l", col = "grey", lwd=4)


################################################################################
# now, let's finally compute the probability of extinction in 100(=Tmax) years
Nend <- Nr[Tmax,]# extract the last row

#  this tells us in which specific replicate the population got extinct 
which(Nend==0)  

# well... you see, they are way too many...
length(which(Nend==0))  # let's just track how many got extinct

#this should be simply equal to Nrep, in fact I could have used Nrep directly
#whatever, now we can compute what  the fraction of replicates in which the
#population got extinct - which is our probability of extinction after Tmax
#years
length(Nend) 

ext.prob <- length(which(Nend==0))/length(Nend); ext.prob*100 

################################################################################
# A comment: extinction is just extinction, no doubt about this, zero
# individual, period. But if the population drops below a very low number, say
# 10 individuals, it is not going to be good anyway, because we do not want to
# wait for the population to actually go extinct to take action, but we shall
# move way before.... So, it might be important to compute not just the
# probability to drop to zero, but the probability to drop below a given
# population threshold. This is not true extinction, and that's why we call
# these thresholds>0 "quasi-extinction thresholds" So, let's compute the
# quasi-extinction risk first, let's define a wide range of quasi extinction
# thresholds
# 
qet <- seq(0,max(Nend), by = 1) 

# then define the function to compute the quasi extinction risk 
qetf <- function (qet) length(which(Nend<=qet))/length(Nend)

# then "apply" this function for each "qet", namely: 
ext.prob  <-  sapply(qet, qetf) # NB sapply uses a vector and returns a vector
plot(ext.prob~qet, type = "l", ylim = c(0,1), 
          xlab = "quasi-extinction threshold",
          ylab = "probability") # natural scale

# we are actually more interested to derive probability of extinction at low
# abundances so, let's plot this graph in a semilogarithic scale
plot(ext.prob~qet, type = "l", log = "x") 

# compute mean population size at each time step (i.e. the mean of each row). 
# The "apply" function applies the function "mean" to each row (margin = 1) of
# the matrix Nr
Nmean  <- apply(Nr, MARGIN=1, mean) 
plot(Nmean) # this is the mean at time t for the Nrep replicates
points(N, type = "l", col = "red") # this was the deterministic solution!
# you can see that with so many replicates, the  mean of the stochastic process
# overlaps closely to the deterministic solution, as it should be

################################################################################
# So, expected population size increases in a Malthusian way 

# Question:
# Does it mean that the extinction risk decreases with time?
# well, let's see it...

# let's define a function to compute, at each time step (each row of the matrix
# Nr) the fraction of replicates in which the population got extinct (zero
# individuals):

roe.f <-function (Nvect) length(which(Nvect==0))/Nrep

# let's apply this function to each row of Nr (MARGIN = 1)
roe  <-  apply(Nr, MARGIN = 1, roe.f) 
plot(roe, xlab = "time", ylab = "risk of extinction")


# Multiple choice test (select the right answer)
# while the mean population size increases, the risk of extinction : 
#   a) also increases with time
#   b) decreases with time
#   c) remains constant with time
# 

################################################################################
#another interesting statistic is in fact the population size distribution at
#any given time step
gf <- par(mfrow=c(3,2), mai = c(0.3,0.35,0.25,0.25))  
hist(Nr[5,])  # here it is population size distribution at year 5 
hist(Nr[10,]) # here it is population size distribution at year 10 
hist(Nr[20,]) # ...and so on and so for
hist(Nr[30,])
hist(Nr[40,])
hist(Nr[Tmax,])

par(gf)
# we might be interested to compute some statistical property of this
# distribution for instance the quantiles at each time step
q0.975 = apply(Nr, 1, quantile, probs = 0.975)
q0.025 = apply(Nr, 1, quantile, probs = 0.025)
q0.750 = apply(Nr, 1, quantile, probs = 0.75)
q0.250 = apply(Nr, 1, quantile, probs = 0.25)

# and then make a plot of it 
plot(q0.975,   type = "l", col = "grey", lwd=2, xlab="time")
points(q0.025, type = "l", col = "grey", lwd=2)
points(q0.750, type = "l", col = "green", lwd=2)
points(q0.250, type = "l", col = "green", lwd=2)
points(Nmean,  type = "l", col = "red", lwd = 4)

# this is a nicer way to plot it (although, with ggplot2 now you can do so much
# better!)

xt = seq(1,Tmax, by=1)
plot(Nmean, type = "l", col = "red", lwd = 4, ylim= c(0,max(q0.975)))
polygon(c(xt, rev(xt)),c(q0.025, rev(q0.975)), col="light gray", border = NA)
polygon(c(xt, rev(xt)),c(q0.250, rev(q0.750)), col="dark gray", border = NA)
points(Nmean, type = "l", col = "red", lwd = 4) # redraw the data
axis(1, xlim=c(1,Tmax+1)) # redraw the x axis
 
################################################################################
# Finally, let's compute the time to reach a quasi extinction threshold
################################################################################
# we have to find for each replicate the time at which the pop.size drops below
# a "quasi-extinction threshol" (qet) and pick up the first one in the list,
# namely the first time the population reached the quasi extinction threshold:
ttqetf1 <- function(Nsim, qet) which(Nsim<=qet)[1] 

# here we set the quasi extinction threshold for instance to 1, instead of zero 
ttqet   <- apply(Nr, MARGIN=2,ttqetf1, qet=1) 

# "ttqet" is a vector whose elements report for each replicate the time when the
# population reached the quasi-extinction threshold, if ever!: if, within a
# specific replicate, the population never reached the quasi-extinction
# threshold, the function returns NA

# we can now compute the statistical properties of the distribution of the time
# taken to drop at or below "qet". Note that na.rm is to remove "NA" values.
quantile(ttqet, na.rm = TRUE) 

# manually change qet to 1, 2, etc., to see how time to extinction change 

# now let's do it in one shot for all the "qet" of interest
# Let's look just at the time taken for the population to drop at or below a 
# threshold between 0 and 9:
ttqetf2   <- function(qet) apply(Nr, 2,ttqetf1, qet=qet)
qet       <- seq(0,9, by=1)  
ttqet.all <- sapply(qet, ttqetf2) # apply this function for each threshold
time.quantile  <- apply(ttqet.all, 2, quantile, na.rm = TRUE); time.quantile

# and plot the result
plot(time.quantile[5,]~qet, type = "l", ylim=c(0,Tmax+1), col="black", 
  ylab="Time to quasi extinction threshold", xlab="quasi-exctintion threshold")
points(time.quantile[4,]~qet, type = "l", col="green")
points(time.quantile[3,]~qet, type = "l", col="red", lwd=4)
points(time.quantile[2,]~qet, type = "l", col="green")
points(time.quantile[1,]~qet, type = "l", col="black")

################################################################################
# check whether the extinction probability depends upon the initial conditions,
# i.e. the initial number of individual in the population
################################################################################

# here below is the function to generate random deviates of number of offspring
# generated by a grey reef shark
fv <- c(3,4,5,6)
fp <- c(0.3, 0.4, 0.2, 0.1)
fv[sample.int(n=4, size =1, prob = fp)]
xv <- sapply(1:1000, function(z){fv[sample.int(n=4, size =1, prob =fp)]})
mean(xv)


 

