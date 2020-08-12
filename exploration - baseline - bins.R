patients_d1 <- vector("list", length = 36)
library(readr)
library(data.table)
library(fastICA)

for (i in c(1:36)){
	path <- paste("./EMG_data_for_gestures-master", i, "d1.txt", sep = "/")
	print(path)
	patients_d1[[i]] <- read.delim(path)
}

train <- patients_d1[c(1:36)]

train_ics <- vector("list", length = 36)
unmix_matrices <- vector("list", length = 36)


for (i in c(1:36)){
	
	mixture_signal <- train[[i]][,2:9]
	unmix <- fastICA(mixture_signal, 8, alg.typ = "parallel", fun = "logcosh", alpha = 1,
									 method = "C", row.norm = FALSE, maxit = 200,
									 tol = 0.000001, verbose = TRUE)
	new_record <- data.table(unmix$S)
	colnames(new_record) <- c("IC1", "IC2", "IC3", "IC4", "IC5", "IC6", "IC7", "IC8")
	train_ics[[i]] <- new_record
	unmix_matrices[[i]] <- data.table(unmix$W)
	
}




mixtures_gestures <- vector("list",200)
ics_gestures <- vector("list",200)
s <- 0

for (i in c(1:length(train))){
	pat <- train[[i]]
	ica_pat <- train_ics[[i]]
	signal_end <- 0
	for (j in c(3:dim(pat)[1])){
		signal_start <- signal_end + 1
		if (pat$class[j-1] != pat$class[j-2]){
			s <- s + 1
			signal_end <- j-1
			new_mix_gesture <- data.frame(pat[signal_start:signal_end,c(1:10)])
			new_ica_gesture <- data.frame(ica_pat[signal_start:signal_end,c(1:8)])
			new_ica_gesture <- cbind(new_ica_gesture, pat[signal_start:signal_end, 1])
			new_ica_gesture <- cbind(new_ica_gesture, pat[signal_start:signal_end, 10])
			
			colnames(new_ica_gesture)[9] <- "Time"
			colnames(new_ica_gesture)[10] <- "Class"
			mixtures_gestures[[s]] <- data.table(new_mix_gesture)
			ics_gestures[[s]] <- data.table(new_ica_gesture)
		}
	}
}


length(ics_gestures)

for (i in c(1:length(ics_gestures))){
		b <- trunc((dim(ics_gestures[[i]])[1]*.25))
		c <- trunc((dim(ics_gestures[[i]])[1]*.75))
		ics_gestures[[i]] <- ics_gestures[[i]][b:c]
}



##########################################################################
#CLASS DISTRIBUTION

class_instances <- c(0,0,0,0,0,0,0)
for (j in c(1:length(mixtures_gestures))){
	for (i in c(1:7)){
		if (mixtures_gestures[[j]]$class[2] == i){
			class_instances[i] = class_instances[i] + 1
		}
	}
}

class_instances[1] <- class_instances[1] + 432
class_instances_binary <- c(0,0)
class_instances_binary[1] <- class_instances[1]
class_instances_binary[2] <- length(mixtures_gestures) - class_instances_binary[1]


counts <- class_instances

df=data.frame(counts)
library(data.table)
df <- setDT(df, keep.rownames = TRUE)[]
df$rn <- as.integer(df$rn)
df$rn <- as.character(df$rn)
colnames(df) <- c("class", "instances")

library(ggplot2)

par(mfrow = c(1,2))
plot1 <- ggplot(data = df, aes(x=class, y=instances)) +
	geom_bar(stat="identity", fill="steelblue")+
	geom_text(aes(label=instances), vjust=1.6, color="white", size=3.5)+
	theme_minimal()


counts <- class_instances_binary

df=data.frame(counts)
library(data.table)
df <- setDT(df, keep.rownames = TRUE)[]
df$rn <- as.integer(df$rn)
df$rn <- as.character(df$rn)
colnames(df) <- c("class", "instances")
plot2 <- ggplot(data = df, aes(x=c("No Gesture", "Gesture"), y=instances)) +
	geom_bar(stat="identity", fill="steelblue", width = 0.6)+
	geom_text(aes(label=instances), vjust=1.6, color="white", size=3.5)+
	labs(x = "class") +
	theme_minimal()

#install.packages("gridExtra")
library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)


##########################################################################
##########################################################################



#### BASELINE #######

mix_rms_dt_new_baseline <- data.table(mean_rms = numeric(), class = numeric()) 

for (i in c(1:length(mixtures_gestures))){
	new_record_mix <- 0
	for (j in c(1:8)){
		a <- paste("channel", j, sep = "")
		rms_mix <- sqrt((1/dim(mixtures_gestures[[i]])[1]) * sum(mixtures_gestures[[i]][,get(a)]^2, na.rm = TRUE))
		new_record_mix = new_record_mix + rms_mix
	}

	mean_rms <- new_record_mix/8
	class <- mixtures_gestures[[i]]$class[1]
	df <- data.frame(mean_rms, class)
	mix_rms_dt_new_baseline <- rbind(mix_rms_dt_new_baseline, df)


}


set.seed(1)
idx.tr <- sample(1:dim(mix_rms_dt_new_baseline)[1], round(dim(mix_rms_dt_new_baseline)[1]*.75), replace=F)
#idx.tr
tr <- mix_rms_dt_new_baseline[idx.tr]
te <- mix_rms_dt_new_baseline[-idx.tr]

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
glmfit <- multinom(class ~ mean_rms, data = tr, maxit = 500)

a <- predict(glmfit, te[,1])
a <- list(a)
mean((a) == te[,2])
table(factor(a[[1]], levels=min(0):max(7)), 
			factor(te$class, levels=min(0):max(7)))




tr <- mix_rms_dt_new_baseline[idx.tr]
te <- mix_rms_dt_new_baseline[-idx.tr]

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
glmfit <- multinom(class ~ mean_rms, data = tr, maxit = 500)

a <- predict(glmfit, te[,1])
a <- list(a)
mean((a) == te[,2])
table(factor(a[[1]], levels=min(0):max(7)), 
			factor(te$class, levels=min(0):max(7)))




tr <- mix_rms_dt_new_baseline[idx.tr]
te <- mix_rms_dt_new_baseline[-idx.tr]

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
glmfit <- multinom(class ~ mean_rms, data = tr, maxit = 500)

a <- predict(glmfit, te[,1])
a <- list(a)
mean((a) == te[,2])
table(factor(a[[1]], levels=min(1):max(2)), 
			factor(te$class, levels=min(1):max(2)))

##################################################################################

ics_rms_dt <- data.table(rms1 = numeric(), rms2 = numeric(), rms3 = numeric(), 
												 rms4 = numeric(), rms5 = numeric(), rms6 = numeric(), 
												 rms7 = numeric(), rms8 = numeric(), class = numeric())

mixtures_rms_dt <- data.table(rms1 = numeric(), rms2 = numeric(), rms3 = numeric(), 
												 rms4 = numeric(), rms5 = numeric(), rms6 = numeric(), 
												 rms7 = numeric(), rms8 = numeric(), class = numeric())


for (i in c(1:length(ics_gestures))){
	
	ics_record <-ics_gestures[[i]]
	mixture_record <- mixtures_gestures[[i]]
	new_record <- data.frame()
	new_record_mix <- data.frame()

	for (j in c(1:8)){
		a <- paste("IC", j, sep = "")
		rms <- sqrt((1/dim(ics_gestures[[i]])[1]) * sum(ics_gestures[[i]][,get(a)]^2, na.rm = TRUE))
		new_record <- rbind(new_record, rms)
		
		
		a <- paste("channel", j, sep = "")
		rms_mix <- sqrt((1/dim(mixtures_gestures[[i]])[1]) * sum(mixtures_gestures[[i]][,get(a)]^2, na.rm = TRUE))
		new_record_mix<- rbind(new_record_mix, rms_mix)
	}
	
	new_record <- rbind(new_record, ics_gestures[[i]]$Class[1])
	new_record <- data.frame(t(new_record))
	row.names(new_record) <- i
	ics_rms_dt <- rbind(ics_rms_dt, new_record, use.names = FALSE)
	
	
	new_record_mix <- rbind(new_record_mix, mixtures_gestures[[i]]$class[1])
	new_record_mix <- data.frame(t(new_record_mix))
	rownames(new_record_mix) <- i
	mixtures_rms_dt <- rbind(mixtures_rms_dt, new_record_mix, use.names = FALSE)
	
}



#install.packages("OneR")
library(OneR)
new_df <- bin(ics_rms_dt, nbins = 10, labels = c("bin1", "bin2", "bin3", "bin4", "bin5", "bin6"
																											, "bin7", "bin8", "bin9", "bin10"))

bins.dt <- data.table("bin1" = numeric(), "bin2" = numeric(), "bin3" = numeric(), 
											"bin4" = numeric(), "bin5" = numeric(), "bin6" = numeric(), 
											"bin7" = numeric(), "bin8" = numeric(), "bin9" = numeric(), 
											"bin10" = numeric() , 
											"class" = numeric() )

dim(bins.dt)[2]

for (k in c(1:dim(new_df)[1])){
	vec <- vector("integer", dim(bins.dt)[2])
	for (j in c(1:8)){
		a <- paste("rms", j, sep = "")
		
		for (i in c(1:dim(bins.dt)[2] - 1)){
			b <- paste("bin", i, sep = "")
			if (new_df[k,a, with = FALSE] == b){
				vec[i] <- vec[i] + 1
			}
		}
		
	}
	
	vec[dim(bins.dt)[2]] <- as.numeric(new_df[k]$class) - 1
	df <- data.frame(t(vec))
	for (m in c(1:dim(bins.dt)[2] - 1)){
	 		b <- paste("bin", m, sep = "")
	 		colnames(df)[m] <- b
	}
	
	colnames(df)[dim(bins.dt)[2]] <- "class"
	bins.dt <- rbind(bins.dt, df)
}


bins_.dt <- bins.dt


set.seed(1)
idx.tr <- sample(1:length(ics_gestures), round(length(ics_gestures)*.75), replace=F)
tr <- bins_.dt[idx.tr]
te <- bins_.dt[!idx.tr]

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
glmfit <- multinom(class ~ bin1 + bin2 + bin3 + bin4 + bin5 + bin6 + bin7 + bin8 + bin9 + bin10
									 	, data = tr, maxit = 500)

a <- predict(glmfit, te[,1:10])
a <- list(a)
mean((a) == te[,11])
table(factor(a[[1]], levels=min(0):max(7)), 
			factor(te$class, levels=min(0):max(7)))





tr <- bins_.dt[idx.tr]
te <- bins_.dt[!idx.tr]

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
glmfit <- multinom(class ~ bin1 + bin2 + bin3 + bin4 + bin5 + bin6 + bin7 + bin8 + bin9 + bin10
									 	, data = tr, maxit = 500)

a <- predict(glmfit, te[,1:10])
a <- list(a)
mean((a) == te[,11])
table(factor(a[[1]], levels=min(0):max(7)), 
			factor(te$class, levels=min(0):max(7)))





tr <- bins_.dt[idx.tr]
te <- bins_.dt[!idx.tr]

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
glmfit <- multinom(class ~ bin1 + bin2 + bin3 + bin4 + bin5 + bin6 + bin7 + bin8 + bin9 + bin10
									 	, data = tr, maxit = 500)

a <- predict(glmfit, te[,1:10])
a <- list(a)
mean((a) == te[,11])
table(factor(a[[1]], levels=min(1):max(2)), 
			factor(te$class, levels=min(1):max(2)))








