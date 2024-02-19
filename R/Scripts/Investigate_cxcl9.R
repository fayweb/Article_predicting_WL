##############later

lab %>%
    ggplot(aes(x = CXCL9, y = WL_max)) + 
    geom_jitter() + 
    geom_smooth(method = "lm")

model <- lm(WL_max ~ CXCL9, data = lab)
summary(model)

model <- lm(WL_max ~ IFNy + CXCR3 + IL.6 + IL.13 +
                IL1RN + CASP1 + CXCL9 + IDO1 + IRGM1 + MPO + 
                MUC2 + MUC5AC + MYD88 + NCR1 + PRF1 + RETNLB + SOCS1 + 
                TICAM1 + TNF , data = lab)
summary(model)



ggplot(Field, aes(x = HI, y = CXCL9, color = Sex)) +
    geom_jitter() +
    geom_smooth(aes(fill = Sex), 
                method = "loess", se = TRUE, alpha = 0.2)
model <- lm(CXCL9 ~ Sex*HI, data = Field)    
summary(model)    