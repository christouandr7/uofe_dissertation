patients_d1 <- vector("list", length = 36)
library(readr)
library(data.table)
library(fastICA)

#load dataset - mixtures
for (i in c(1:36)){
	path <- paste("./EMG_data_for_gestures-master", i, "d1.txt", sep = "/")
	print(path)
	patients_d1[[i]] <- read.delim(path)
}

train <- patients_d1[c(1:36)]

train_ics_multiplied <- vector("list", length = 36)
unmix_matrices <- vector("list", length = 36)

#run fastica only on one random subject
mixture_signal <- train[[1]][,2:9]
unmix <- fastICA(mixture_signal, 8, alg.typ = "parallel", fun = "logcosh", alpha = 1,
								 method = "C", row.norm = FALSE, maxit = 200,
								 tol = 0.000001, verbose = TRUE)
unmix_matrix <- unmix$W	#save the unmixing matrix to be used on all other subjects

#multiply all subjects' signals with the above unmixing matrix
for (i in c(1:36)){
	
	new_record <- data.table(as.matrix(train[[i]][,2:9] ) %*% as.matrix(unmix_matrix))
	colnames(new_record) <- c("IC1", "IC2", "IC3", "IC4", "IC5", "IC6", "IC7", "IC8")
	train_ics_multiplied[[i]] <- new_record
}

#initialize 2 lists - mixtures and ics will be splited into gestures
#lists will contain the gestures
#see when gesture label changes during the signal and define it as the end of the gesture
mixtures_gestures <- vector("list",200)
ics_gestures <- vector("list",200)
s <- 0

for (i in c(1:length(train))){
	pat <- train[[i]]
	ica_pat <- train_ics_multiplied[[i]]
	signal_end <- 0
	for (j in c(3:dim(pat)[1])){
		signal_start <- signal_end + 1
		if (pat$class[j-1] != pat$class[j-2]){
			s <- s + 1
			signal_end <- j-1
			new_mix_gesture <- data.frame(pat[signal_start:signal_end,c(1:10)])
			new_ica_gesture <- data.frame(ica_pat[signal_start:signal_end,c(1:8)])
			new_ica_gesture <- cbind(new_ica_gesture, pat[signal_start:signal_end, 9])
			new_ica_gesture <- cbind(new_ica_gesture, pat[signal_start:signal_end, 10])
			
			colnames(new_ica_gesture)[9] <- "Time"
			colnames(new_ica_gesture)[10] <- "Class"
			mixtures_gestures[[s]] <- data.table(new_mix_gesture)
			ics_gestures[[s]] <- data.table(new_ica_gesture)
		}
	}
}

#Windows method - use only the middle 50% of the signal
length(ics_gestures)

for (i in c(1:length(ics_gestures))){
	b <- trunc((dim(ics_gestures[[i]])[1]*.25))
	c <- trunc((dim(ics_gestures[[i]])[1]*.75))
	ics_gestures[[i]] <- ics_gestures[[i]][b:c]
}


#datatable with the rms values - no ordering ambiguity since they have the same order 
ics_rms_dt <- data.table(rms1 = numeric(), rms2 = numeric(), rms3 = numeric(), 
												 rms4 = numeric(), rms5 = numeric(), rms6 = numeric(), 
												 rms7 = numeric(), rms8 = numeric(), class = numeric())

#calculate the rms values for each gesture
for (i in c(1:length(ics_gestures))){
	
	ics_record <-ics_gestures[[i]]
	new_record <- data.frame()
	
	for (j in c(1:8)){
		a <- paste("IC", j, sep = "")
		rms <- sqrt((1/dim(ics_gestures[[i]])[1]) * sum(ics_gestures[[i]][,get(a)]^2, na.rm = TRUE))
		new_record <- rbind(new_record, rms)
		
	}
	
	new_record <- rbind(new_record, ics_gestures[[i]]$Class[1])
	new_record <- data.frame(t(new_record))
	row.names(new_record) <- i
	ics_rms_dt <- rbind(ics_rms_dt, new_record, use.names = FALSE)
	
}


#Run the 3 classification problems with the same train and test set
set.seed(1)
idx.tr <- sample(1:length(ics_gestures), round(length(ics_gestures)*.75), replace=F)

#1st classification problem - original problem
tr <- ics_rms_dt[idx.tr]
te <- ics_rms_dt[!idx.tr]

for (i in c(1:dim(tr)[1])){
	if (tr[i]$class == 0){
		tr[i]$class = 1
	}
}

for (i in c(1:dim(te)[1])){
	if (te[i]$class == 0){
		te[i]$class = 1
	}
}

library(nnet)
glmfit <- multinom(class ~ rms1 + rms2 + rms3 + rms4 + rms5 + rms6 
									 + rms7 + rms8
									 , data = tr, maxit = 500)

a <- predict(glmfit, te[,1:8])
a <- list(a)
mean((a) == te[,9])
table(factor(a[[1]], levels=min(0):max(7)), 
			factor(te$class, levels=min(0):max(7)))




#2nd classification problem (gestures only)
tr <- ics_rms_dt[idx.tr]
te <- ics_rms_dt[!idx.tr]

for (i in c(1:dim(tr)[1])){
	if (tr[i]$class == 0){
		tr[i]$class = 1
	}
}

for (i in c(1:dim(te)[1])){
	if (te[i]$class == 0){
		te[i]$class = 1
	}
}

tr <- tr[tr$class != 1]
te <- te[te$class != 1]

library(nnet)
glmfit <- multinom(class ~ rms1 + rms2 + rms3 + rms4 + rms5 + rms6 
									 + rms7 + rms8
									 , data = tr, maxit = 500)

a <- predict(glmfit, te[,1:8])
a <- list(a)
mean((a) == te[,9])
table(factor(a[[1]], levels=min(0):max(7)), 
			factor(te$class, levels=min(0):max(7)))






#3rd classification problem (binary)
tr <- ics_rms_dt[idx.tr]
te <- ics_rms_dt[!idx.tr]


for (i in c(1:dim(tr)[1])){
	if (tr[i]$class == 0){
		tr[i]$class = 1
	}
	if (tr[i]$class != 0 & tr[i]$class != 1){
		tr[i]$class = 2
	}
	
}

for (i in c(1:dim(te)[1])){
	if (te[i]$class == 0){
		te[i]$class = 1
	}
	
	if (te[i]$class != 0 & te[i]$class != 1){
		te[i]$class = 2
	}
}

library(nnet)
glmfit <- multinom(class ~ rms1 + rms2 + rms3 + rms4 + rms5 + rms6 
									 + rms7 + rms8
									 , data = tr, maxit = 500)

a <- predict(glmfit, te[,1:8])
a <- list(a)
mean((a) == te[,9])
table(factor(a[[1]], levels=min(1):max(2)), 
			factor(te$class, levels=min(1):max(2)))











