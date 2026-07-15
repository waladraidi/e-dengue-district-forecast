source('./Model/R/99_load.R')


args <- commandArgs(trailingOnly = TRUE)
j <- as.numeric(args[1])  
k <- as.numeric(args[2])

all.mods <- list('mod1' = mod1, 'mod2' = mod2, 'mod3' = mod3,'mod4' = mod4,'mod5' = mod5)


#j=59

#k=5


 #modN_extract <- as.numeric(str_match(names(all.mods)[k], "mod(\\d+)")[1, 2])
 
  #mod1 <- inla_spacetime_mod(vintage_date = date.test2[j], formula1 = all.mods[[k]], modN=modN_extract ) 

# 
# 
for (j in 1:length(date.test2)) {
  for (k in 1:5) {
    tryCatch({
      modN_extract <- as.numeric(stringr::str_match(names(all.mods)[k], "mod(\\d+)")[1, 2])

      mod1 <- inla_spacetime_mod(
        vintage_date = date.test2[j],
        formula1     = all.mods[[k]],
        modN         = modN_extract
      )

    }, error = function(e) {
      message(paste("Skipping iteration j =", j, "k =", k, "due to error:", e$message))
    })
  }
}
