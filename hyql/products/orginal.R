# numerical implementation in R for replicating the values
# in tables 1 to 4 in Hull/White 2012: CVA and wrong-way risk
# code for candidate project "Counterparty Credit Risk" at the University of Copenhagen
# Christian Klausen

#======================
# simulation parameters
#======================
discretization.steps <- 10	# N in the text
simulations <- 500000		# n in the text
iterations <- 2

#=================
# model parameters 
#=================
maturity <- 1
x0 <- 1	
sigma <- 0.15	# volatility of FX rate (no drift)
strike <- 1 
s <- 0.0125 	# 125 bps spread
R <- 0.4 		# 60 % recovery rate
r <- 0.05 		# 5% interest rate
mm <- 100 		# notional of FX swap in millions 

#=================
# table parameters
#=================
b <- 0.03		# b parameter
long = TRUE		# TRUE when buying contract, FALSE when selling
c <- 15		# cure period in days (or margin period of risk)
K <- 0		# collateral treshold. ignored if collateralized=FALSE
collateralized=TRUE # set to FALSE for no collateral, else TRUE

#==================
#==================

epsilon <- 0.000000015  # epsilon for spread sensitivities
epsilon_x <- 0.002 	# epsilon for fx sensitivities


factor <- ifelse(long,1,-1)
dt <- maturity/discretization.steps
curetime <- c/365
curetime_fac <- round(curetime/dt)
tolerance <-  10^-30 #.Machine$double.eps

datavec <- matrix(NA,iterations,5)
system.time(for(m in 1:iterations){

# CVA 
q2 <- a_hzrd   <- rep(0,discretization.steps+1)
hzrd <- q <- exposure <- matrix(0, discretization.steps+1,simulations)
portfolio_value <- rep(0, simulations)
X  <- rep(x0, simulations)



#delta s / gamma s values
s_plus <- s+epsilon
s_minus <- s-epsilon
hzrd_deltas_plus <- q_deltas_plus <- hzrd_deltas_minus <- q_deltas_minus <-  matrix(0, discretization.steps+1,simulations)
a_hzrd_deltas_plus <- a_hzrd_deltas_minus <- rep(0,discretization.steps+1)

#delta x / gamma x values
x_plus <- x0+epsilon_x
x_minus <- x0-epsilon_x
X_plus  <- rep(x_plus, simulations)
X_minus <- rep(x_minus,simulations)
hzrd_deltax_plus <- q_deltax_plus <- exposure_plus <- exposure_minus <- matrix(0, discretization.steps+1,simulations)
hzrd_deltax_minus <- q_deltax_minus <- q_deltax_b0 <- matrix(0, discretization.steps+1,simulations)
a_hzrd_deltax_minus <- a_hzrd_deltax_plus  <- rep(0,discretization.steps+1)
collateral <-collateral_minus <- collateral_plus <- matrix(max(0,-K), discretization.steps+1,simulations)

#hazard functions
hzrd_a_quick  <- function(x,w,sk,tk){
	return(mean(exp(-dt*(exp(x+b*w))))-exp(-sk*tk/(1-R)))
}

hzrd_a_long  <- function(x,w,sk,tk,counter,hzrd_matrix){
	return(mean(exp(-dt*(colSums(hzrd_matrix[1:(counter-1),])+exp(x+b*w))))-exp(-sk*tk/(1-R)))
}


#collateral function
collateral_fun <- function(w,k){
	return(pmax(w-k,0))
}


#=========================
# Simulation of CVA values
#=========================

portfolio_value <- mm*factor*(X-strike)
portfolio_value_x_plus <- mm*factor*(X_plus-strike)
portfolio_value_x_minus <- mm*factor*(X_minus-strike)
collateral[1+curetime_fac, ] <- collateral_fun(portfolio_value,K)
collateral_plus[1+curetime_fac, ] <- collateral_fun(portfolio_value_x_plus,K)
collateral_minus[1+curetime_fac, ] <- collateral_fun(portfolio_value_x_minus,K)



# Calculation of the first step. Seperate treatment, as it is only half of the usual step length
vec <- rnorm(simulations)
X <- X*exp(-0.5*0.5*dt*sigma^2+sqrt(dt*0.5)*sigma*vec)
X_plus <- X_plus*exp(-0.5*0.5*dt*sigma^2+sqrt(dt*0.5)*sigma*vec)
X_minus <- X_minus*exp(-0.5*0.5*dt*sigma^2+sqrt(dt*0.5)*sigma*vec)

portfolio_value <- exp(-r*(maturity-0.5*dt))*mm*factor*(X-strike)
portfolio_value_x_plus <- exp(-r*(maturity-0.5*dt))*mm*factor*(X_plus-strike)
portfolio_value_x_minus <- exp(-r*(maturity-0.5*dt))*mm*factor*(X_minus-strike)


# Hazard rates for calculation of CVA^b
collateral[2+curetime_fac, ] <- collateral_fun(portfolio_value,K)
a_hzrd[2] <- uniroot( hzrd_a_quick, c(-20,10),w=portfolio_value , sk=s, tk=dt, tol = tolerance)$root
hzrd[2,] <- exp(a_hzrd[2]+b*portfolio_value)

# Hazard rates for changes in spread for calculation of CVA^b_s
a_hzrd_deltas_plus[2] <- uniroot( hzrd_a_quick, c(-20,10),w=portfolio_value , sk=s_plus, tk=dt, tol = tolerance)$root
hzrd_deltas_plus[2,] <- exp(a_hzrd_deltas_plus[2]+b*portfolio_value)
a_hzrd_deltas_minus[2] <- uniroot( hzrd_a_quick, c(-20,10),w=portfolio_value , sk=s_minus, tk=dt, tol = tolerance)$root
hzrd_deltas_minus[2,] <- exp(a_hzrd_deltas_minus[2]+b*portfolio_value)

# Hazard rates for changes in exchange rate for calculation of CVA^b_x
collateral_plus[2+curetime_fac, ] <- collateral_fun(portfolio_value_x_plus,K)
a_hzrd_deltax_plus[2] <- uniroot( hzrd_a_quick, c(-20,10),w=portfolio_value_x_plus , sk=s, tk=dt, tol = tolerance)$root
hzrd_deltax_plus[2,] <- exp(a_hzrd_deltax_plus[2]+b*portfolio_value_x_plus)


collateral_minus[2+curetime_fac, ] <- collateral_fun(portfolio_value_x_minus,K)
a_hzrd_deltax_minus[2] <- uniroot( hzrd_a_quick, c(-20,10),w=portfolio_value_x_minus , sk=s, tk=dt, tol = tolerance)$root
hzrd_deltax_minus[2,] <- exp(a_hzrd_deltax_minus[2]+b*portfolio_value_x_minus)


# q_1 as in HW12 for the standard wrong-way/right-way risk case
q[2,] <- 1-exp(-dt*hzrd[2,])

# q_1 as in HW12 for the independent case
q2[2] <- 1-exp(-dt*s/(1-R))

# q_1 as in HW12 for small changes in spread and fx rate
q_deltas_plus[2,] <- 1-exp(-dt*hzrd_deltas_plus[2,])
q_deltas_minus[2,] <- 1-exp(-dt*hzrd_deltas_minus[2,])
q_deltax_plus[2,] <- 1-exp(-dt*hzrd_deltax_plus[2,])
q_deltax_minus[2,] <- 1-exp(-dt*hzrd_deltax_minus[2,])

# calculation of discounted exposure
if(collateralized){
	exposure[2,] <- exp(-r*dt*0.5)*(pmax(portfolio_value-collateral[2,],0))
	exposure_plus[2,] <- exp(-r*dt*0.5)*(pmax(portfolio_value_x_plus-collateral_plus[2,],0))
	exposure_minus[2,] <- exp(-r*dt*0.5)*(pmax(portfolio_value_x_minus-collateral_minus[2,],0))
} else {
	exposure[2,] <- exp(-r*dt*0.5)*(pmax(portfolio_value,0))
	exposure_plus[2,] <- exp(-r*dt*0.5)*(pmax(portfolio_value_x_plus,0))
	exposure_minus[2,] <- exp(-r*dt*0.5)*(pmax(portfolio_value_x_minus,0))

}


for(i in 3:(discretization.steps+1)){
	vec <- rnorm(simulations)
	
	X <- X*exp(-0.5*dt*sigma^2+sqrt(dt)*sigma*vec)
	X_plus <- X_plus*exp(-0.5*dt*sigma^2+sqrt(dt)*sigma*vec)
	X_minus <- X_minus*exp(-0.5*dt*sigma^2+sqrt(dt)*sigma*vec)

	portfolio_value <- exp(-r*(maturity-(i-1.5)*dt))*mm*factor*(X-strike)
	portfolio_value_plus <- exp(-r*(maturity-(i-1.5)*dt))*mm*factor*(X_plus-strike)
	portfolio_value_minus <- exp(-r*(maturity-(i-1.5)*dt))*mm*factor*(X_minus-strike)


	if(collateralized){
		tryCatch({collateral[i+curetime_fac,] <- collateral_fun(portfolio_value,K)
		collateral_plus[i+curetime_fac,] <- collateral_fun(portfolio_value_plus,K)
		collateral_minus[i+curetime_fac,] <- collateral_fun(portfolio_value_minus,K)
		}, error = function(e){})
		exposure[i,] <- exp(-r*(dt*(i-1)-0.5*dt))*(pmax(portfolio_value-collateral[i,],0))
		exposure_plus[i,] <- exp(-r*(dt*(i-1)-0.5*dt))*(pmax(portfolio_value_plus-collateral_plus[i,],0))
		exposure_minus[i,] <- exp(-r*(dt*(i-1)-0.5*dt))*(pmax(portfolio_value_minus-collateral_minus[i,],0))

	} else {
		exposure[i,] <- exp(-r*(dt*(i-1)-0.5*dt))*(pmax(portfolio_value,0))
		exposure_plus[i,] <- exp(-r*(dt*(i-1)-0.5*dt))*(pmax(portfolio_value_plus,0))
		exposure_minus[i,] <- exp(-r*(dt*(i-1)-0.5*dt))*(pmax(portfolio_value_minus,0))
	}

	a_hzrd[i] <- uniroot( hzrd_a_long, c(-10,5),w=portfolio_value , sk=s, tk=(i-1)*dt,counter=i, hzrd_matrix=hzrd, tol = tolerance)$root
	hzrd[i,] <- exp(a_hzrd[i]+b*portfolio_value)
	q[i,] <- exp(-dt*colSums(hzrd[1:(i-1),]))-exp(-dt*colSums(hzrd[1:i,]))

	
	#delta s hzrd matrices
	a_hzrd_deltas_plus[i] <- uniroot( hzrd_a_long, c(-10,5),w=portfolio_value , sk=s_plus, tk=(i-1)*dt,counter=i,hzrd_matrix=hzrd_deltas_plus, tol = tolerance)$root
	hzrd_deltas_plus[i,] <- exp(a_hzrd_deltas_plus[i]+b*portfolio_value)
	q_deltas_plus[i,] <- exp(-dt*colSums(hzrd_deltas_plus[1:(i-1),]))-exp(-dt*colSums(hzrd_deltas_plus[1:i,]))

	a_hzrd_deltas_minus[i] <- uniroot( hzrd_a_long, c(-10,5),w=portfolio_value , sk=s_minus, tk=(i-1)*dt,counter=i,hzrd_matrix=hzrd_deltas_minus, tol = tolerance)$root
	hzrd_deltas_minus[i,] <- exp(a_hzrd_deltas_minus[i]+b*portfolio_value)
	q_deltas_minus[i,] <- exp(-dt*colSums(hzrd_deltas_minus[1:(i-1),]))-exp(-dt*colSums(hzrd_deltas_minus[1:i,]))

	#delta x hzrd matrices
	a_hzrd_deltax_plus[i] <- uniroot( hzrd_a_long, c(-10,5),w=portfolio_value_plus , sk=s, tk=(i-1)*dt,counter=i, hzrd_matrix=hzrd_deltax_plus, tol = tolerance)$root
	hzrd_deltax_plus[i,] <- exp(a_hzrd_deltax_plus[i]+b*portfolio_value_plus)
	q_deltax_plus[i,] <- exp(-dt*colSums(hzrd_deltax_plus[1:(i-1),]))-exp(-dt*colSums(hzrd_deltax_plus[1:i,]))

	a_hzrd_deltax_minus[i] <- uniroot( hzrd_a_long, c(-10,5),w=portfolio_value_minus , sk=s, tk=(i-1)*dt,counter=i, hzrd_matrix=hzrd_deltax_minus, tol = tolerance)$root
	hzrd_deltax_minus[i,] <- exp(a_hzrd_deltax_minus[i]+b*portfolio_value_minus)
	q_deltax_minus[i,] <- exp(-dt*colSums(hzrd_deltax_minus[1:(i-1),]))-exp(-dt*colSums(hzrd_deltax_minus[1:i,]))
	q2[i] <- exp(-(i-2)*dt*s/(1-R))-exp(-(i-1)*dt*s/(1-R))
	

}
#=====================================
# Calculation of CVA and sensitivities
#==================================== 

#=====
# CVA
#=====
CVA_dummy <- expo <- rep(0,discretization.steps+1)
for (i in 2:(discretization.steps+1)){
	CVA_dummy[i] <- mean(q[i,]*exposure[i,])
	expo[i] <- mean(exposure[i,])
}


CVA0 <- sum(q2*expo)*(1-R)
CVA <- sum(CVA_dummy)*(1-R)

delta0 <- deltas_plus <- deltas_minus <- deltas3 <- rep(0,discretization.steps+1)
gammas0 <-  gammas <- rep(0,(discretization.steps+1))
deltax0 <- deltax <- rep(0,discretization.steps+1)
gammax0 <- gammax <- rep(0,discretization.steps+1)

for(i in 2:(discretization.steps+1)){

#=============
# delta spread
#=============	
	
	delta0[i] <- ( (i-1)*dt*exp(-s*(i-1)*dt/(1-R))-(i-2)*dt*exp(-s*(i-2)*dt/(1-R))   )*expo[i]
	deltas_plus[i] <- mean(q_deltas_plus[i,]*exposure[i,])
	deltas_minus[i] <- mean(q_deltas_minus[i,]*exposure[i,])

#=============
# gamma spread
#=============
	gammas0[i] <- ( ((i-2)*dt)^2*exp(-s*(i-2)*dt/(1-R))-((i-1)*dt)^2*exp(-s*(i-1)*dt/(1-R))   )*expo[i]
	gammas[i] <- mean(q_deltas_plus[i,]*exposure[i,]-2*q[i,]*exposure[i,]+q_deltas_minus[i,]*exposure[i,])/(epsilon^2)
	
#===========
# fx  delta
#===========

	deltax0[i] <- q2[i]*(mean(exposure_plus[i,])-mean(exposure_minus[i,]))
	deltax[i] <- mean(q_deltax_plus[i,]*exposure_plus[i,]-q_deltax_minus[i,]*exposure_minus[i,])


#===========
# fx  gamma
#===========
	gammax0[i] <- (q2[i]*(mean(exposure_plus[i,])+mean(exposure_minus[i,])-2*expo[i]))/(epsilon_x^2)
	gammax[i] <- (mean(q_deltax_plus[i,]*exposure_plus[i,]+ q_deltax_minus[i,]*exposure_minus[i,])-2*CVA_dummy[i])/(epsilon_x^2)
}
deltas <- (1-R)*(sum(deltas_plus)-sum(CVA_dummy))/epsilon
deltas0 <- sum(delta0)

gammas0 <- (1/(1-R))*sum(gammas0)
gammas <- (1-R)*sum(gammas)

deltax0 <- (1-R)*sum(deltax0)/(2*epsilon_x)
deltax <- (1-R)*sum(deltax)/(2*epsilon_x)

gammax0 <- (1-R)*sum(gammax0)
gammax <- (1-R)*sum(gammax)

datavec[m,] <- c((CVA/CVA0-1)*100,(deltax/deltax0 - 1)*100,(sum(gammax)/sum(gammax0)  - 1)*100,(deltas/deltas0 - 1)*100,(gammas/gammas0 - 1)*100)
	print(m)
})


data <- matrix(NA,5,3)

dummy <- sort(datavec[,1], TRUE)
CVA_mean <- mean(datavec[,1])
CVA_left <- dummy[95]
CVA_right <- dummy[5]
data[1,]  <- c(CVA_left, CVA_mean, CVA_right)

dummy <- sort(datavec[,2], TRUE)
delta_fx_mean <- mean(datavec[,2])
delta_fx_mean_left <- dummy[95]
delta_fx_mean_right <- dummy[5]
data[2,] <- c(delta_fx_mean_left,delta_fx_mean, delta_fx_mean_right)

dummy <- sort(datavec[,3], TRUE)
gamma_fx_mean <- mean(datavec[,3])
gamma_fx_mean_left <- dummy[95]
gamma_fx_mean_right <- dummy[5]
data[3,] <- c(gamma_fx_mean_left,gamma_fx_mean,gamma_fx_mean_right)

dummy <- sort(datavec[,4], TRUE)
delta_s_mean <- mean(datavec[,4])
delta_s_mean_left <- dummy[95]
delta_s_mean_right <- dummy[5]
data[4,] <- c(delta_s_mean_left, delta_s_mean,delta_s_mean_right)

dummy <- sort(datavec[,5], TRUE)
gamma_s_mean <- mean(datavec[,5])
gamma_s_mean_left <- dummy[95]
gamma_s_mean_right <- dummy[5]
data[5,] <- c(gamma_s_mean_left,gamma_s_mean,gamma_s_mean_right)

#datavec
data
