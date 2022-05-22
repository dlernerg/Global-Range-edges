library(dominanceanalysis)

var.mod <- c(2,3,8,12,13,16,20,26,32,33,36,40) 
#var.mod <- c(21:31) 

alltog <- alltog.in2[,c(6,7,(var.mod+7))]

names.var <- paste0("absolute BC", var.mod[1:6])
names.var[7] <- paste("absolute", "elevation")
names.var[8:11] <- paste0("SH BC", (var.mod[8:11]-20))
names.var[12] <- "SH elevation"
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

########prepare da.model figure
da2 <- (da.model$contribution.average$r2.m) %>% as.data.frame()
da2$variable <- c(names.var,"ac") %>% factor(., levels = unique(. ))

ggplot(da2, aes(x = variable, y = .)) + 
  geom_col(fill = "grey", color = "black", width = 0.7, alpha = 0.5)+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab(bquote(mean ~ R^2))

ggsave(filename = "dominance(succ_fail) autocor.png", width = 4.76, height = 3.4)
  
a + theme_classic()
#1,2,3,8,12,13,20,2,3,6,7,12,13,20,