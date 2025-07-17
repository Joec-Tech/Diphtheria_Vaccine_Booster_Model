##CHOLERA DATA CLEANING
library(openxlsx)
matilda_cholera_data <- read.xlsx("cholera_outbreak.xlsx")
matilda_cholera_data

##Cleanin categorical variables

names(matilda_cholera_data)

unique(matilda_cholera_data$Sex)
unique(matilda_cholera_data$Education)
unique(matilda_cholera_data$Marital_Status)
unique(matilda_cholera_data$Locality)
unique(matilda_cholera_data$Status)
unique(matilda_cholera_data$Test_Result)

library(gtsummary)
matilda_cholera_data%>% 
  tbl_summary()

cholera_clean <- matilda_cholera_data %>% 
  mutate(Gender = recode(Sex, 
                         "femil" = "Female",
                         "femle" = "Female",
                         "feminist" = "Female",
                         "mal" = "Male",
                         "masculine" = "Male"))

# To view the unique values after renaming
unique(cholera_clean$Gender)




