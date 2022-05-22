library(dominanceanalysis)

var.mod <- c(1:11)
#var.mod <- c(12:20)

alltog <- alltog.in2[,c(6,7,(var.mod+7))]

names.var <- paste0("absolute BC", var.mod[1:6])
names.var[7] <- paste("absolute", "elevation")
names.var[8:11] <- paste0("SH BC", (var.mod[8:11]-20))
names(alltog) <- c("countss","species.t", names.var)


f <- cbind(alltog$countss,alltog$species.t)~ . + ac
model.all <- glm(f, data = alltog, 
                 family = binomial)

summary(model.all)

da.model <- dominanceAnalysis(model.all)
dominanceMatrix(da.model, type = "general", ordered = T)


dominanceMatrix(da.model, type = "conditional", ordered = T)

contributionByLevel(da.model, fit.functions = "r2.m")


plot(da.model, which.graph ="conditional",fit.function = "r2.m")
