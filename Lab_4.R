install.packages("jsonlite")
library("jsonlite")
library("ggplot2")
library("readr")
setwd("C:/Users/Kiruthika Velusamy/Documents/Sem-2/ASM")
df_bus_TO <- stream_in(file("Business_Toronto_Restaurant.json"))
names(df_bus_TO)
?dim
head(df_bus_TO)
attach(df_bus_TO)
dim(df_bus_TO)
df_bus_TO

sel_nbhd <- df_bus_TO$neighborhood == "Chinatown" | df_bus_TO$neighborhood == "Koreatown" 
class(df_bus_TO$categories)
head(df_bus_TO$categories)
cat_total <- unlist(df_bus_TO$categories)
class(cat_total)
cat_total <- factor(cat_total)
class(cat_total)
cat_names_sort <- sort(table(cat_total), decreasing = TRUE)
head(cat_names_sort, n = 25)
tail(cat_names_sort, n = 25)
cat_names <- names(cat_names_sort)[2:18]
?sapply()
cat_bus_ind_mat <- sapply(df_bus_TO$categories, function(y) as.numeric(cat_names %in% y))
?t()
cat_bus_ind_mat <- t(cat_bus_ind_mat)
?colnames()
colnames(cat_bus_ind_mat) <- cat_names
df_TO_tidy_cat <- cbind(df_bus_TO, cat_bus_ind_mat)
names(df_TO_tidy_cat)
?apply
head(df_TO_tidy_cat)
attach(df_TO_tidy_cat)




# Neighborhood data frame for  Scarborough and  Etobicoke
df_ck <- subset(df_TO_tidy_cat, (neighborhood == "Scarborough" | neighborhood == "Etobicoke")&Indian==1&is_open==1, select =  c("stars", "neighborhood","Indian","is_open"))
e<-df_ck$neighborhood=="Etobicoke"
length(e)
df_ck$index<-ifelse(df_ck$neighborhood=='Scarborough',1,2)
df_ck$index<-as.factor(df_ck$index)
dim(df_ck)
head(df_ck)
attach(df_ck)
library(ggplot2)
ggplot(df_ck) + geom_boxplot(aes(index, stars, fill = index)) + geom_jitter(aes(index, stars, shape = df_ck$index))

tapply(df_ck$stars,df_ck$neighborhood, mean)
tapply(df_ck$stars, df_ck$neighborhood, median)
tapply(df_ck$stars,df_ck$neighborhood, sd)
t.test(stars ~ neighborhood, data=df_ck, var.equal = TRUE)









compare_2_gibbs <- function(y, ind, mu0 = 2.5, tau0 = 1/1.265625, del0 = 0, gamma0 = 1/1.265625, a0 = 2, b0 = 2.53125, maxiter = 5000)
{
  y1 <- y[ind == 1]
  y2 <- y[ind == 2]
  
  n1 <- length(y1) 
  n2 <- length(y2)
  
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  
  for(s in 1 : maxiter) 
  {
    
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    
    ##update mu
    taun <-  tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    
    ##update del
    gamman <-  gamma0 + tau*(n1 + n2)
    deln <- ( del0 * gamma0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}


library("MCMCpack")
fit <- compare_2_gibbs(df_ck$stars, as.factor(df_ck$index))
plot(as.mcmc(fit))
apply(fit$params, 2, mean)

raftery.diag(as.mcmc(fit))

y1_sim <- rnorm(5000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
y2_sim <- rnorm(5000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))
ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))
mean(y1_sim > y2_sim)
1/sd(y1_sim)
mean(y2_sim)
mean((y1_sim+y2_sim)/2)
ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) + geom_abline(slope = 1, intercept = 0)




# question 1 part 2 comparing multiple groups
df_open<-subset(df_bus_TO, is_open==1, select =c("stars", "neighborhood","is_open"))
names(df_open)
?apply
head(df_open)
attach(df_open)
df_open
dim(df_open)

drop1<- df_open$neighborhood==""
df_open2 <- df_open[!drop1, ]
dim(df_open2)

df_open2$number <- as.numeric(as.factor(df_open2$neighborhood))
df_open2
newdata<-df_open2[order(df_open2$number),]
newdata <- newdata[!newdata$neighborhood=="Meadowvale Village",]
newdata <- newdata[!newdata$neighborhood=="Cooksville",]
library("ggplot2")

ggplot(df_open2) + geom_boxplot(aes(x = reorder(number, stars,median), stars, fill = reorder(number, stars, median)), show.legend=FALSE)

ggplot(df_open2, aes(x = reorder(number, neighborhood, length))) + stat_count()
ggplot(df_open2, aes(stars)) + stat_bin()




ggplot(data.frame(size = tapply(df_open2$stars, df_open2$number, length), mean_score = tapply(df_open2$stars, df_open2$number, mean)), aes(size, mean_score)) + geom_point()


compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
  ### weakly informative priors
  a0 <- 2 ; b0 <- 2.53125 ## tau_w hyperparameters
  eta0 <- 2 ; t0 <- 2.53125 ## tau_b hyperparameters
  mu0<- 2.5  ; gamma0 <- 1/(1.265625)
  ###
  
  ### starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
      
    }
   
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
      
    }
    
    bn <- b0 + ss/2
    
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta-mu)^2)/2
    tau_b <- rgamma(1, etam, tm)
    
    
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(theta_mat) <- levels(ind)
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}


newdata$number<-as.factor(newdata$number)
fit2 <- compare_m_gibbs(newdata$stars, newdata$number)
apply(fit2$params, 2, mean)
apply(fit2$params, 2, sd)
mean(1/sqrt(fit2$params[, 3]))
levels(newdata$number)
theta_hat <- apply(fit2$theta, 2,mean)
which(theta_hat==min(theta_hat))
newdata[newdata$number==41,]
high<- newdata$number==41
newdata[high, ]
mean(fit2$theta[,41]>mean(theta_hat))

y_bar= tapply(df_open2$stars, df_open2$number, mean)-theta_hat
ggplot(data.frame(size = tapply(df_open2$stars, df_open2$number, length), y_bar = y_bar), aes(size, y_bar)) + geom_point()
ggplot(, aes(x = reorder(neighborhood, neighborhood, length))) + stat_count()


# data analysis for review count  and rating
drop1<- df_bus_TO$neighborhood==""
df_open2 <- df_bus_TO[!drop1, ]
dim(df_open2)


#inference from data 
#  How class is related to review count 
# hypothesis 1: higher number of review count is contributing to good class

df_Review<-subset(df_open2, is_open==1 , select =c("stars", "neighborhood","is_open","review_count"))
df_review_reserve<-cbind(df_Review, df_bus_TO[df_bus_TO$is_open==1 & !(df_bus_TO$neighborhood==""),]$attributes$RestaurantsReservations)
df_review_reserve

median<-median(df_review_reserve$stars)
df_review_reserve$class<-ifelse(df_review_reserve$stars>=median,"Good","Bad")
colnames(df_review_reserve)[5] <- "Reservation"
df_review_reserve<-df_review_reserve[!is.na(df_review_reserve$Reservation),]
attach(df_review_reserve)
dim(df_review_reserve)

#library("ggmosaic")
#ggplot(data = df_review_reserve) +
#  geom_mosaic(aes(x = product(class), fill=Reservation), na.rm=TRUE) +
 # labs(x="class ", title='class vs reservation')


ggplot(df_review_reserve) + geom_boxplot(aes(class, review_count, fill = class)) 

 # feature selection:
#taken original dataframe with categories converted from df_bus_TO
class(df_bus_TO$categories)
head(df_bus_TO$categories)
cat_total <- unlist(df_bus_TO$categories)
class(cat_total)
cat_total <- factor(cat_total)
class(cat_total)
cat_names_sort <- sort(table(cat_total), decreasing = TRUE)
head(cat_names_sort, n = 25)
tail(cat_names_sort, n = 25)
cat_names <- names(cat_names_sort)[2:60]


cat_bus_ind_mat <- sapply(df_bus_TO$categories, function(y) as.numeric(cat_names %in% y))

cat_bus_ind_mat <- t(cat_bus_ind_mat)

colnames(cat_bus_ind_mat) <- cat_names
df_TO_tidy_cat <- cbind(df_bus_TO, cat_bus_ind_mat)
names(df_TO_tidy_cat)
?apply
head(df_TO_tidy_cat)
attach(df_TO_tidy_cat)
head(df_TO_tidy_cat$state)

# given new variable for second question part 2(a) for factors contributing to closed resta..
df_closed_features<-df_TO_tidy_cat
dim(df_closed_features)

# data imputation -50 % data are missing values

# removing the name, address, business id ,categories  from the dataframe.
# separate opening and closing hrs inorder to determine whether the restaruant is being closed or not.
ggplot(Tidy_df, aes(y=df_closed_features$postal_code, x=factor(df_closed_features$postal_code))) +
  geom_boxplot()
head(df_bus_TO$postal_code)
sum(is.na(df_closed_features$categories))
sum(df_closed_features$neighborhood=="")
names(df_closed_features)
attach(df_closed_features)
df_closed_features[,c("name","business_id","address","categories")] <- list(NULL)
colnames(df_closed_features)

names(df_closed_features)
df_closed_features$hours

# splitting hrs as per day basis
df_closed_features$state<-as.numeric(df_closed_features$state)
plot(df_closed_features$is_open,df_closed_features$state)

head(df_closed_features$state)
df_Monday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Monday),'-',fixed=TRUE)))
df_closed_features$df_Monday_open<-df_Monday$X1
df_closed_features$df_Monday_closed<-df_Monday$X2


df_Tuesday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Tuesday),'-',fixed=TRUE)))
df_closed_features$df_Tuesday_open<-df_Tuesday$X1
df_closed_features$df_Tuesday_closed<-df_Tuesday$X2

df_Wednesday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Wednesday),'-',fixed=TRUE)))
df_closed_features$df_Wednesday_open<-df_Wednesday$X1
df_closed_features$df_Wednesday_closed<-df_Wednesday$X2

df_Thursday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Thursday),'-',fixed=TRUE)))
df_closed_features$df_Thursday_open<-df_Thursday$X1
df_closed_features$df_Thursday_closed<-df_Thursday$X2

df_Friday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Friday),'-',fixed=TRUE)))
df_closed_features$Friday_open<-df_Friday$X1
df_closed_features$Friday_closed<-df_Friday$X2

df_Saturday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Saturday),'-',fixed=TRUE)))
df_closed_features$df_Saturday_open<-df_Saturday$X1
df_closed_features$df_Saturday_closed<-df_Saturday$X2

df_Sunday <- data.frame(do.call('rbind', strsplit(as.character(df_closed_features$hours$Sunday),'-',fixed=TRUE)))
df_closed_features$df_Sunday_open<-df_Sunday$X1
df_closed_features$df_Sunday_closed<-df_Sunday$X2





df_closed_features

names(df_closed_features$attributes)
df_attributes<-df_closed_features$attributes
df_BusinessParking<-df_closed_features$attributes$BusinessParking
df_HairSpecializesIn<-df_closed_features$attributes$HairSpecializesIn
#names(df_closed_features$attributes$BusinessParking)
#names(df_closed_features$attributes$HairSpecializesIn)
#names(df_closed_features$attributes$Music)
#names(df_closed_features$attributes$Ambience)
df_Music<-df_closed_features$attributes$Music
df_Ambience<-df_closed_features$attributes$Ambience
df_BestNights<-df_closed_features$attributes$BestNights
df_GoodForMeal<-df_closed_features$attributes$GoodForMeal
df_DietaryRestrictions<-df_closed_features$attributes$DietaryRestrictions
df_Feature_Selection<-cbind(df_closed_features,df_attributes,df_BusinessParking,df_HairSpecializesIn,df_Music,df_Ambience,df_BestNights,df_GoodForMeal,df_DietaryRestrictions)
# dropped hours and attributes (attributes dataframe combined with main data frame)
df_Feature_Selection[,"hours"] = NULL;
df_Feature_Selection[,"attributes"] = NULL;
df_Feature_Selection[,"BusinessParking"] = NULL;
df_Feature_Selection[,"HairSpecializesIn"] = NULL;
df_Feature_Selection[,"Music"] = NULL;
df_Feature_Selection[,"Ambience"] = NULL;
df_Feature_Selection[,"BestNights"] = NULL;
df_Feature_Selection[,"GoodForMeal"] = NULL;
df_Feature_Selection[,"DietaryRestrictions"] = NULL;

names(df_Feature_Selection)

colnames(df_Feature_Selection)[111]<-"BusinessParking_garage"
colnames(df_Feature_Selection)[112]<-"BusinessParking_street"
colnames(df_Feature_Selection)[113]<-"BusinessParking_validated"
colnames(df_Feature_Selection)[114]<-"BusinessParking_lot"
colnames(df_Feature_Selection)[115]<-"BusinessParking_valet"
colnames(df_Feature_Selection)[116]<-"HairSpecializesIn_coloring"
colnames(df_Feature_Selection)[117]<-"HairSpecializesIn_curly"
colnames(df_Feature_Selection)[118]<-"HairSpecializesIn_perms"
colnames(df_Feature_Selection)[119]<-"HairSpecializesIn_kids"
colnames(df_Feature_Selection)[120]<-"HairSpecializesIn_extensions"

names(df_closed_features$attributes$Music)

colnames(df_Feature_Selection)[121]<-"Music_dj"
colnames(df_Feature_Selection)[122]<-"Music_background_music"
colnames(df_Feature_Selection)[123]<-"Music_no_music"
colnames(df_Feature_Selection)[124]<-"Music_karaoke"
colnames(df_Feature_Selection)[125]<-"Music_live"
colnames(df_Feature_Selection)[126]<-"Music_video"
colnames(df_Feature_Selection)[127]<-"Music_jukebox"

names(df_closed_features$attributes$BestNights)



colnames(df_Feature_Selection)[128]<-"Ambience_romantic"
colnames(df_Feature_Selection)[129]<-"Ambience_intimate"
colnames(df_Feature_Selection)[130]<-"Ambience_classy"
colnames(df_Feature_Selection)[131]<-"Ambience_hipster"
colnames(df_Feature_Selection)[132]<-"Ambience_touristy"
colnames(df_Feature_Selection)[133]<-"Ambience_trendy"
colnames(df_Feature_Selection)[134]<-"Ambience_upscale"

colnames(df_Feature_Selection)[135]<-"Ambience_casual"


names(df_closed_features$attributes$BestNights)
colnames(df_Feature_Selection)[136]<-"BestNights_monday"
colnames(df_Feature_Selection)[137]<-"BestNights_tuesday"
colnames(df_Feature_Selection)[138]<-"BestNights_friday"
colnames(df_Feature_Selection)[139]<-"BestNights_wednesday"
colnames(df_Feature_Selection)[140]<-"BestNights_thursday"
colnames(df_Feature_Selection)[141]<-"BestNights_sunday"
colnames(df_Feature_Selection)[142]<-"BestNights_saturday"
names(df_closed_features$attributes$GoodForMeal)
colnames(df_Feature_Selection)[143]<-"GoodForMeal_dessert"
colnames(df_Feature_Selection)[144]<-"GoodForMeal_latenight"
colnames(df_Feature_Selection)[145]<-"GoodForMeal_lunch"
colnames(df_Feature_Selection)[146]<-"GoodForMeal_dinner"
colnames(df_Feature_Selection)[147]<-"GoodForMeal_breakfast"
colnames(df_Feature_Selection)[148]<-"GoodForMeal_brunch"
names(df_closed_features$attributes$DietaryRestrictions)
colnames(df_Feature_Selection)[149]<-"DietaryRestrictions_dairy-free"
colnames(df_Feature_Selection)[150]<-"DietaryRestrictions_gluten-free"
colnames(df_Feature_Selection)[151]<-"DietaryRestrictions_vegan"
colnames(df_Feature_Selection)[152]<-"DietaryRestrictions_kosher"
colnames(df_Feature_Selection)[153]<-"DietaryRestrictions_halal"
colnames(df_Feature_Selection)[154]<-"DietaryRestrictions_soy-free"
colnames(df_Feature_Selection)[155]<-"DietaryRestrictions_vegetarian"

names(df_Feature_Selection)
Missing<-sort(sapply(df_Feature_Selection,function(a) sum(is.na(a), decreasing=TRUE)))
write.csv(Missing,"missing.csv")

# check all null values 
sort(sapply(df_Feature_Selection,function(a) sum(is.na(a), decreasing=TRUE)))
cond1 <- sapply(df_Feature_Selection, function(col) sum(is.na(col))<3000)
df_feature2<-df_Feature_Selection[, cond1,drop=FALSE]


names(df_feature2)

# state is not necessary since state ontario
# city is not necessary since it contains only toronto
# Postal code is removed since it contains lesser variation in values & less impact on is_open
# removing food since it is type of categories 

names(df_feature2)

df_feature2[,"city"] = NULL;
df_feature2[,"state"] = NULL;
df_feature2[,"postal_code"] = NULL;
df_feature2[,"Food"] = NULL;

head(df_feature2)
names(df_feature2)
# hypothesis: A rest is closed if a bad food, service or nightlife
# co-relation between bars and night life is 0.979
# removing highly co-related variable bars
library("ggpubr")
cor(df_feature2$Nightlife,df_feature2$Bars , method = c("pearson", "kendall", "spearman"))
df_feature2[,"Bars"] = NULL;

names(df_feature2)
Tidy_df<-df_feature2
Tidy_df[30:63] = NULL;
head(Tidy_df)
names(Tidy_df)
cor(df_feature2$Nightlife,df_feature2$Pubs , method = c("pearson", "kendall", "spearman"))

#tidy plot after removing categories feom ferom pubs to tepas_bars
names(Tidy_df)

#plot to infer nightlife, pubs, is_open
# if it is sushibar and rating between 3 to 3.5 , rest is closed.
# if it is chinese and rating between app.2.8 to 3.5 , rest is closed.
## if it is burgers and rating between 3 to 3.5 , rest is closed.
## if it is asiaan_fusion and rating between 3 to 3.5 , rest is closed.
# for pizza, there is diffe. in the mean of being open and closed
attach(Tidy_df)
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Nightlife), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$`Sushi Bars`), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Chinese), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$`Asian Fusion`), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Pizza), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Indian), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Thai), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Mexican), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Seafood), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Friday_closed), fill=factor(is_open))) +
  geom_boxplot()
names(Regul)

# removing  variables that doesnt have significance on rest being closed
#Nightlife,Canadian(new),Sanwidches,Breakfast &brunch, Italien,Cafes,Coffe &tea
# Fast food, Japanese , Middle Eastern.Mediterranean, Korean,pubs
sum(is.na(Tidy_df$Pubs))
Remove_list<-cbind("Nightlife","Canadian (New)","Sandwiches","Breakfast & Brunch","Italian","Cafes","Coffee & Tea","Fast Food","Japanese","Middle Eastern","Mediterranean","Korean","Pubs")
Tidy_df[,Remove_list]=NULL

class(Tidy_df$df_Saturday_closed)

names(Tidy_df)

cor(as.numeric(Tidy_df$df_Monday_open),as.numeric(Tidy_df$df_Tuesday_open) , method = c("pearson", "kendall", "spearman"))
Tidy_df$df_Monday_open<-as.numeric(df_Monday_open)
class(df_Monday_open)

data <- Tidy_df[,cbind("df_Monday_open","df_Tuesday_open","df_Wednesday_open","df_Thursday_open","df_Saturday_open","Friday_open","df_Sunday_open")]
data$df_Tuesday_open<-as.numeric(df_Tuesday_open)
data$df_Wednesday_open<-as.numeric(df_Wednesday_open)
data$df_Thursday_open<-as.numeric(df_Thursday_open)
data$Friday_open<-as.numeric(Friday_open)
data$df_Sunday_open<-as.numeric(df_Sunday_open)
data$df_Saturday_open<-as.numeric(df_Saturday_open)
class(data$df_Saturday_open)
data$df_Monday_open<-as.numeric(data$df_Monday_open)
class(data)

Heatmap_open <- round(cor(data,use="complete.obs"),2)
head(Heatmap_open)
library(reshape2)
melted_cormat <- melt(Heatmap_open)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# heatmap for closed hrs 

data1 <- Tidy_df[,cbind("df_Monday_closed","df_Tuesday_closed","df_Wednesday_closed","df_Thursday_closed","df_Saturday_closed","Friday_closed","df_Sunday_closed")]
names(data1)
data1$df_Tuesday_closed<-as.numeric(df_Tuesday_closed)
data1$df_Wednesday_closed<-as.numeric(df_Wednesday_closed)
data1$df_Thursday_closed<-as.numeric(df_Thursday_closed)
data1$Friday_closed<-as.numeric(Friday_closed)
data1$df_Sunday_closed<-as.numeric(df_Sunday_closed)
data1$df_Saturday_closed<-as.numeric(df_Saturday_closed)
class(data$df_Saturday_open)
data1$df_Monday_closed<-as.numeric(data1$df_Monday_closed)


Heatmap_closed <- round(cor(data1,use="complete.obs"),2)
head(Heatmap_closed)
library(reshape2)
melted_cormat1 <- melt(Heatmap_closed)
library(ggplot2)
ggplot(data = melted_cormat1, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

sum(is.na(as.numeric(Tidy_df$Friday_closed)))
sum(is.na(as.numeric(Tidy_df$df_Thursday_closed)))

sum(is.na(as.numeric(Tidy_df$Friday_open)))
sum(is.na(as.numeric(Tidy_df$df_Saturday_open)))
#  since it is highly correlated  with all of the weekday opening and closing hrs.
#except  remove other things . Since weekend and weekdays differs , one variable from weekdays and one from weekend are chosen.

Remove_open_hrs<-cbind("df_Monday_open","df_Tuesday_open","df_Wednesday_open","df_Thursday_open","df_Sunday_open")
Tidy_df[,Remove_open_hrs]=NULL

Remove_closed_hrs<-cbind("df_Monday_closed","df_Tuesday_closed","df_Wednesday_closed","df_Saturday_closed","df_Sunday_closed")
Tidy_df[,Remove_closed_hrs]=NULL

names(Tidy_df)

# removing BusinessAcceptsCreditCards because it doesnt have any imapct at all
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BusinessAcceptsCreditCards), fill=factor(is_open))) +
  geom_boxplot()

#  there is impact  RestaurantsPriceRange2

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsPriceRange2), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForKids), fill=factor(is_open))) +
  geom_boxplot()
# there is some difference in being closed 
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Alcohol), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$HasTV), fill=factor(is_open))) +
  geom_boxplot()
# there is some difference in noise level in being closed (loud)
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$NoiseLevel), fill=factor(is_open))) +
  geom_boxplot()
#there is some diff
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsAttire), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsGoodForGroups), fill=factor(is_open))) +
  geom_boxplot()


ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsGoodForGroups), fill=factor(is_open))) +
  geom_boxplot()

# there is some diff in wifi
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$WiFi), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsReservations), fill=factor(is_open))) +
  geom_boxplot()
#ther is some diff

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsTakeOut), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$RestaurantsTableService), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BikeParking), fill=factor(is_open))) +
  geom_boxplot()
# there is some difference
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BusinessParking_garage), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BusinessParking_street), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BusinessParking_validated), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BusinessParking_lot), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$BusinessParking_valet), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_romantic), fill=factor(is_open))) +
  geom_boxplot()

# there is some variation in mean of open
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_intimate), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_classy), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_hipster), fill=factor(is_open))) +
  geom_boxplot()
#there is differ
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_touristy), fill=factor(is_open))) +
  geom_boxplot()
#there is some variation in mean of open
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_trendy), fill=factor(is_open))) +
  geom_boxplot()
#there is some variation in mean of open
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_upscale), fill=factor(is_open))) +
  geom_boxplot()
#there is some variation in mean of open
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_casual), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForMeal_dessert), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForMeal_latenight), fill=factor(is_open))) +
  geom_boxplot()

ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForMeal_lunch), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForMeal_dinner), fill=factor(is_open))) +
  geom_boxplot()
#there is some variation in mean of open
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForMeal_breakfast), fill=factor(is_open))) +
  geom_boxplot()
ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$GoodForMeal_brunch), fill=factor(is_open))) +
  geom_boxplot()


remove1<-c("GoodForMeal_dessert","GoodForMeal_lunch","BusinessParking_Valet","BusinessParking_validated","GoodForMeal_dinner","GoodForMeal_latenight","BikeParking","GoodForMeal_brunch","Ambience_classy","Ambience_hipster","lot","OutdoorSeating","RestaurantsTableService","RestaurantsReservations","RestaurantsGoodForGroups","NoiseLevel","HasTV","GoodForKids")
Tidy_df[,remove1]=NULL
sort(sapply(Tidy_df,function(a) sum(is.na(a), decreasing=FALSE)))

# business accept credit cards have more values yes so it is more biased .
sum(Tidy_df$WiFi=='no',na.rm=TRUE)
sum(df_feature2$BusinessAcceptsCreditCards==0,na.rm=TRUE)
Tidy_df[,"BusinessParking_garage"]=NULL
Tidy_df[,"BusinessParking_street"]=NULL
# more na values more than 1800 WIFI 2093 BusinessAcceptsCreditCards 1811
# wifi paid values are less than 34 (mean variations) for free and no, there is no mean variation is observed.
Tidy_df[,"WiFi"]=NULL
Tidy_df[,"BusinessAcceptsCreditCards"]=NULL

sort(sapply(Tidy_df,function(a) sum(is.na(a), decreasing=FALSE)))

# removing other variables except ambience becuase other varibales have very less value compared to ambience casual.

sum(Tidy_df$Ambience_romantic==1,na.rm=TRUE)
sum(Tidy_df$Ambience_intimate==1,na.rm=TRUE)
sum(Tidy_df$Ambience_touristy==1,na.rm=TRUE)
sum(Tidy_df$Ambience_trendy==1,na.rm=TRUE)
sum(Tidy_df$Ambience_upscale==1,na.rm=TRUE)
sum(Tidy_df$Ambience_casual==1,na.rm=TRUE)
sum(Tidy_df$BusinessParking_lot==1,na.rm=TRUE)
sum(Tidy_df$RestaurantsAttire==1,na.rm=TRUE)
sum(Tidy_df$`Asian Fusion`==1,na.rm=TRUE)
as.factor(df_feature2$WiFi)
sort(sapply(Tidy_df,function(a) sum(is.na(a), decreasing=FALSE)))
Ambience_removal<-cbind("Ambience_romantic","Ambience_intimate","Ambience_touristy","Ambience_trendy","Ambience_upscale")
Tidy_df[,Ambience_removal]=NULL

# removing restaurant attire variable because box plot 
#of dressy and formal category shows diff in ditribution 
#but the data of categrory dressy and formal are very small in count doesnt have an impact on rest on being closed

Tidy_df[,"RestaurantsAttire"]=NULL
Tidy_df[,"BusinessParking_valet"]=NULL
Tidy_df[,"GoodForMeal_breakfast"]=NULL
Tidy_df[,"BusinessParking_lot"]=NULL

Tidy_df$BusinessAcceptsCreditCards<-df_feature2$BusinessAcceptsCreditCards




names(Tidy_df)
Regul<-Tidy_df
names(Regul)
library("MASS")

tbl<-table(Regul$Chinese, Regul$`Asian Fusion`)
chi2 = chisq.test(tbl, correct=F)
sqrt(chi2$statistic / sum(tbl))

tbl<-table(Regul$Thai, Regul$`Asian Fusion`)
chi2 = chisq.test(tbl, correct=F)
sqrt(chi2$statistic / sum(tbl))

tbl<-table(Regul$Pizza, Regul$Burgers)
chi2 = chisq.test(tbl, correct=F)
sqrt(chi2$statistic / sum(tbl))
# chinese , indian, thai ,chinese are independent of is _open . removing them wouldnt have much impact
tbl<-table(Regul$is_open, Regul$Chinese)
chi2 = chisq.test(tbl, correct=F)
sqrt(chi2$statistic / sum(tbl))

tbl<-table(Regul$is_open, Regul$Indian)
chi2 = chisq.test(tbl, correct=F)
sqrt(chi2$statistic / sum(tbl))

tbl<-table(Regul$is_open, Regul$`Asian Fusion`)
chi2 = chisq.test(tbl, correct=F)
sqrt(chi2$statistic / sum(tbl))

tbl1<-table(Regul$is_open, Regul$Thai)
chi2 = chisq.test(tbl1, correct=F)
sqrt(chi2$statistic / sum(tbl))
#ggplot(Tidy_df, aes(y=Tidy_df$stars, x=factor(Tidy_df$Ambience_touristy), fill=factor(is_open))) +
  geom_boxplot()


tbl = matrix(data=c(55, 45, 20, 30), nrow=2, ncol=2, byrow=T)
dimnames(tbl) = list(City=c('B', 'T'), Gender=c('M', 'F'))
chi2 = chisq.test(tbl, correct=F)
c(chi2$statistic, chi2$p.value)
  
# Since pizza and burgers are fastfood category removing it doesnt not affect . Keep pizza 
Regul[,"Burgers"]=NULL
Regul[,"Thai"]=NULL
Regul[,"Indian"]=NULL
Regul[,"Asian Fusion"]=NULL
  
# asian fusion, chinese, thai,indian can be clustered into one category .
  

# p value is more than 0.74 we can remove it
tbl1<-table(Regul$is_open, Regul$Mexican)
chi2 = chisq.test(tbl1, correct=F)
sqrt(chi2$statistic / sum(tbl1))

am<-table(Regul$is_open, Regul$`American (Traditional)`)
chiam = chisq.test(am, correct=F)
sqrt(chiam$statistic / sum(am))

# p value is more than 0.74 we can remove it
sea<-table(Regul$is_open, Regul$Seafood)
chi2 = chisq.test(tbl1, correct=F)
sqrt(chi2$statistic / sum(sea))
# p value is more than 0.74 we can remove it
tbl1<-table(Regul$is_open, Regul$`Sushi Bars`)
chi2 = chisq.test(tbl1, correct=F)
sqrt(chi2$statistic / sum(tbl))
names(Regul)
#ggplot(data = df_feature2) +geom_mosaic(aes(x = product(), fill=Reservation), na.rm=TRUE)

Regul[,"Sushi Bars"]=NULL
Regul[,"Seafood"]=NULL
Regul[,"Mexican"]=NULL

names(Regul)



b<-table(Regul$is_open, Regul$BusinessAcceptsCreditCards)
bi = chisq.test(b, correct=F)
sqrt(bi$statistic / sum(b))
# Rest delivery can be removed since it is independant of is_open p value is more than 0.05, therefore accepting the null hypo

deliv<-table(Regul$is_open, Regul$RestaurantsDelivery)
deli = chisq.test(deliv, correct=F)
sqrt(deli$statistic / sum(deliv))

restakeout<-table(Regul$is_open, Regul$RestaurantsTakeOut)
chitakeout = chisq.test(restakeout, correct=F)
sqrt(chitakeout$statistic / sum(restakeout))


Regul[,"RestaurantsDelivery"]=NULL
names(Regul)
sort(sapply(Regul,function(a) sum(is.na(a), decreasing=TRUE)))

sort(sapply(df_bus_TO$hours,function(a) sum(is.na(a), decreasing=TRUE)))


###Imputataion KNN

Data<-cbind(df_feature2)
Data
str(Data)
attach(Data)
head(Data)
summary(Data)





library(VIM)
# friday imputing using KNN
impute1 <- kNN(Data, variable="Friday_open",k=8)
summary(impute1)
sum(is.na(impute1$Friday_open))
impute1$Friday_open
Regul$Friday_open<-impute1$Friday_open

impute2 <- kNN(Data, variable="Friday_closed",k=12)
summary(impute2)
Regul$Friday_closed<-impute2$Friday_closed

impute3 <- kNN(Data, variable="df_Saturday_open",k=8)
summary(impute3)

Regul$df_Saturday_open<-impute3$df_Saturday_open

impute4 <- kNN(Data, variable="df_Thursday_closed",k=10)
summary(impute4)
Regul$df_Thursday_closed<-impute4$df_Thursday_closed
impute_al <- kNN(Data, variable="Alcohol",k=3)
summary(impute_al)
Regul$Alcohol<-impute_al$Alcohol

impute_6 <- kNN(Data, variable="RestaurantsPriceRange2",k=4)
summary(impute_al)
Regul$RestaurantsPriceRange2<-impute_6$RestaurantsPriceRange2


impute_7 <- kNN(Data, variable="RestaurantsTakeOut",k=2)
summary(impute_al)
Regul$RestaurantsTakeOut<-impute_7$RestaurantsTakeOut


Final <- kNN(Data, variable="BusinessAcceptsCreditCards",k=2)
summary(Final)
Regul$BusinessAcceptsCreditCards<-Final$BusinessAcceptsCreditCards

sort(sapply(Regul,function(a) sum(is.na(a), decreasing=TRUE)))


Regul$Ambience_casual[(is.na(Regul$Ambience_casual))] = "missing"

# modelling :







require(caTools)
table(Regul)
set.seed(101) 
sample = sample.split(Regul$is_open, SplitRatio = .70)
train_d = subset(Regul, sample == TRUE)
test_d  = subset(Regul, sample == FALSE)




sort(sapply(Regul, function(y) sum(length(which(is.na(y))))))

logitMod <- glm(is_open~neighborhood+longitude+latitude+stars+review_count+Chinese
                +Pizza+`American (Traditional)`+df_Thursday_closed+Friday_open+Friday_closed
                +df_Saturday_open+RestaurantsPriceRange2+Alcohol+RestaurantsTakeOut+Ambience_casual+BusinessAcceptsCreditCards
                ,data=train_d)


logitMod

library(caret)
# Use your model to make predictions, in this example newdata = training set, but replace with your test set   
logitMod$xlevels$neighborhood <- union(logitMod$xlevels$description, levels(train_d$neighborhood))
logitMod$xlevels$df_Thursday_closed <- union(logitMod$xlevels$description, levels(test_d$df_Thursday_closed))
logitMod$xlevels$Friday_open <- union(logitMod$xlevels$description, levels(test_d$Friday_open))
logitMod$xlevels$df_Saturday_open <- union(logitMod$xlevels$description, levels(test_d$df_Saturday_open))
logitMod$xlevels$Friday_closed <- union(logitMod$xlevels$description, levels(test_d$Friday_closed))

pdata <- predict(logitMod, newdata = test_d, type = "response",se.fit=FALSE)


# use caret and compute a confusion matrix
confusionMatrix(data = factor(as.numeric(pdata>0.5)), reference = factor(test_d$is_open))

library('glmnet')

library("glmnetUtils")

regularized <- glmnet(is_open~neighborhood+longitude+latitude+stars+review_count+Chinese
                +Pizza+`American (Traditional)`+df_Thursday_closed+Friday_open+Friday_closed
                +df_Saturday_open+RestaurantsPriceRange2+Alcohol+RestaurantsTakeOut+Ambience_casual+BusinessAcceptsCreditCards
                ,data=train_d,family = "binomial",lambda=0.004)




regularized

library(caret)
# Use your model to make predictions, in this example newdata = training set, but replace with your test set   
regularized$xlevels$neighborhood <- union(regularized$xlevels$description, levels(train_d$neighborhood))
pdata <- predict(regularized, newdata = test_d, type = "response",se.fit=FALSE)


# use caret and compute a confusion matrix
confusionMatrix(data = factor(as.numeric(pdata>0.5)), reference = factor(test_d$is_open))

coef(regularized)
plot(regularized, label = TRUE)

regularized_cv <- glmnet(is_open~neighborhood+longitude+latitude+stars+review_count+Chinese
                      +Pizza+`American (Traditional)`+df_Thursday_closed+Friday_open+Friday_closed
                      +df_Saturday_open+RestaurantsPriceRange2+Alcohol+RestaurantsTakeOut+Ambience_casual+BusinessAcceptsCreditCards
                      ,data=train_d,family = "binomial")

plot(regularized_cv)

# 3 question Mixture models 

library(mclust)
jpeg("Test.jpg")
Mixture<- df_bus_TO[,c("latitude","longitude")]
names(Mixture)
fit <- Mclust(Mixture)
fit$BIC
plot(fit, what = "BIC")
dev.off()

#initial default fit showed that 9 cluster is best component . there might be possiblity that other clusters as well .
# Fitting for 20 groups 
jpeg("fit2.jpg")
fit2 <- Mclust(Mixture, G=1:20)
plot(fit2, what = "BIC")
dev.off()
fit2$BIC
max(fit2$BIC)

# plot for VVE model
jpeg("fit2VVV.jpg")
fit2 <- Mclust(Mixture, G = 2, modelNames = "VVV")
plot(fit2, what = "classification")
dev.off()

# After fit  of 9 clusters showed that 20 cluster is best clusters . there might be possiblity that other clusters as well .
# Fitting for 40 groups 
jpeg("fit3.jpg")
fit3 <- Mclust(Mixture, G=20:30)
fit3$BIC
max(fit3$BIC)
plot(fit3, what = "BIC")
dev.off()

jpeg("fit30to40.jpg")
fit4 <- Mclust(Mixture, G=30:40)
fit4$BIC

plot(fit4, what = "BIC")
dev.off()
max(fit4$BIC)

fit2 <- Mclust(Mixture, G = 35, modelNames = "VVV")
plot(fit2, what = "classification",col=c(68:118))


plot(fit2, what = "uncertainty",col=c(68:118))
table(fit$classification, fit2$classification)



#since high number clusters are obtsained for the fit , its not able visuliza properly.
# so lets take 19 as upper bound


VVV <- Mclust(Mixture, G = 35, modelNames = "VVV")
plot(VVV, what = "classification",col=c(68:110))
df_bus_TO$cluster_fit<-as.factor(VVV$classification)


inter<-tapply(df_bus_TO$stars, df_bus_TO$cluster_fit, mean)
write.csv(inter,"3b.csv")
class(df_bus_TO$attributes$RestaurantsPriceRange2)

price<-tapply(Regul$RestaurantsPriceRange2, df_bus_TO$cluster_fit, mean)
write.csv(price,"3b2.csv")
#  

plot(inter,price)

datacluster<-data.frame(VVV$classification,Regul$is_open)


colnames(datacluster)[1]<-"cluster_fit"
colnames(datacluster)[2]<-"is_open"
names(datacluster)
datacluster$cluster_fit<-as.factor(datacluster$cluster_fit)
datacluster$is_open<-as.factor(datacluster$is_open)
class(datacluster$is_open)
glm1 <- glm(is_open ~ cluster_fit, data =datacluster , family = binomial())
head(predict(glm1))
pred_glm <- plogis(predict(glm1)) ## this is logistic function, maps to [0,1]
boxplot(pred_glm  ~ datacluster$is_open)
table(pred_glm > 0.5, datacluster$is_open)


library(caret)
# confusionMatrix(data = predict(model,json_data_flattened), reference = json_data_flattened$clusters)
tab = table(predict(glm1),datacluster_neigh$cluster_fit)

cm = as.matrix(tab)
n = sum(cm) # number of instances
nc = nrow(cm) # number of classes
diag = diag(cm) # number of correctly classified instances per class 
rowsums = apply(cm, 1, sum) # number of instances per class
colsums = apply(cm, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n

precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1) 

datacluster_neigh<-data.frame(VVV$classification,Regul$neighborhood)

require(foreign)
require(nnet)
colnames(datacluster_neigh)[1]<-"cluster_fit"
colnames(datacluster_neigh)[2]<-"neighborhood"
names(datacluster_neigh)
datacluster_neigh$cluster_fit<-as.factor(datacluster_neigh$cluster_fit)
datacluster_neigh$neighborhood<-as.factor(datacluster_neigh$neighborhood)
class(datacluster$is_open)
model <- multinom(neighborhood ~ cluster_fit, data =datacluster_neigh,MaxNWts = 845544)
head(predict(model))






library(caret)
# confusionMatrix(data = predict(model,json_data_flattened), reference = json_data_flattened$clusters)
tab = table(predict(model,datacluster_neigh),datacluster_neigh$cluster_fit)

cm = as.matrix(tab)
n = sum(cm) # number of instances
nc = nrow(cm) # number of classes
diag = diag(cm) # number of correctly classified instances per class 
rowsums = apply(cm, 1, sum) # number of instances per class
colsums = apply(cm, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n

precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 

data.frame(precision, recall, f1) 
