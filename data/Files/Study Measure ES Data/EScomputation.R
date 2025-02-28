
#Raw and computed ES 

library(dplyr)

#Load dataset--select ESdata.cvs when window opens
ESdata<-read.csv(file.choose())

ESdata1 <- mutate (ESdata,
#1. Compute pooled SD and then Hedges' g from ns, means, and SDs
  Spooled1 = sqrt((SD_t*SD_t*(n_t-1) + SD_c*SD_c*(n_c-1))/(n_t+n_c-2)), 
  g1 = (mean_t - mean_c)/Spooled1, 
  
#2. Apply small sample correction
  ga = g1 - (3*g1 / (4*(n_t + n_c) - 9) ),
  
#3. Compute ns for rs that are converted from means and SDs and ns
  n_r1 = ifelse((is.na(Spooled1)), n_r, (n_t+n_c)),

#4. Convert corrected Hedges' g into rs, if ga is missing, then input entered r, otherwise convert ga to r                                       
   rga = ifelse((is.na(ga)), r, ga/(sqrt((ga*ga) + ((n_r*n_r)/(n_t*n_c)) ) )   ), 

#5. Reverse score rs that are comprised of 1 variable that is in the opposite direction, otherwise same direction
  rgaopp = ifelse((NumOppDir==1), rga*-1, rga),

#6. Reverse score gs that are comprised of 1 variable that is in the opposite direction, otherwise same direction
  gaopp = ifelse((NumOppDir==1), ga*-1, ga),

#7. Compute variance of corrected Hedges' g
  Vga  =       ((n_t + n_c) / (n_t*n_c) ) + ((ga*ga) / (2*(n_t + n_c)) ),

#8. Compute SE of un/corrected Hedges' g
  SEga = sqrt( ((n_t + n_c) / (n_t*n_c) ) + ((ga*ga) / (2*(n_t + n_c)) ) )
  )


# Collapse multiple measures of the same category at the same timepoint within study/tx-control comparison

ESdata2 <- ESdata1 %>% 
  group_by(studyid, TCid, Var1time, var1type, var1mul, Var2time, var2type, var2mul) %>% 
  summarise(rgaopp=mean(rgaopp), gaopp=mean(gaopp), SEga=mean(SEga), Vga=mean(Vga), n_r1=mean(n_r1), n_t=mean(n_t), n_c=mean(n_c)) 

ESdata3 <- ESdata2 %>% 
  group_by(studyid, TCid, Var1time, var1type, Var2time, var2type) %>% 
  summarise(rgaopp=mean(rgaopp), gaopp=mean(gaopp), SEga=mean(SEga), Vga=mean(Vga), n_r1=mean(n_r1), n_t=mean(n_t), n_c=mean(n_c)) 


#Apply reverse rules for r
#Reverse r
#If var1type=1 AND var2type=Group A: O, MA, ML2, ML29  
#If var1type=O, MA, ML2, or ML29 AND var2type=MB, MC, MD, ME
#If var1type=MB, MC, MD, ME AND var2type=O, MA, ML2, or ML29

#Same r
#If var1type=1 AND var2type=Group B: MB, MC, MD, ME, 1 
#If var1type=var2type 
#If both var group A var1type=O, MA, ML2, or ML29 AND var2type=O, MA, ML2, or ML29
#If both var group B var1type=MB, MC, MD, ME, 1 AND var2type=MB, MC, MD, ME, 1

ESdata4 <- mutate (ESdata3,
                     var1grp = ifelse((var1type=="O"|var1type=="MA"|var1type=="ML2"|var1type=="ML29"), "A", "B"),
                     var2grp = ifelse((var2type=="O"|var2type=="MA"|var2type=="ML2"|var2type=="ML29"), "A", "B"),
                     rga = ifelse((var1grp==var2grp), rgaopp, rgaopp*-1),
                     ga = ifelse((var1grp==var2grp), gaopp, gaopp*-1)
)

#Use to create data files for univariate meta-analysis and MASEM.
