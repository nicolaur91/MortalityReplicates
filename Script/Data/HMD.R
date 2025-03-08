#' Download death and exposure data from HMD on England & Wales data
#' 

# Brugernavn ("nsl.fi@cbs.dk") og kode ("laursen03")
us <- "nsl.fi@cbs.dk"
pw <- "laursen03"

# Which countries available?
# getHMDcountries()

#######------- Retrieve Deaths data from England & Wales Total population from HMD ###########---------
# EW 
dat_EW_d <- readHMDweb(CNTRY ="GBRTENW",item = "Deaths_1x1",username = us,
                        password = pw)

#######------- Retrieve Exposure data from 6 low-mortality countries###########---------
# EW 
dat_EW_expc <- readHMDweb(CNTRY ="GBRTENW",item = "Exposures_1x1",username = us,
                           password = pw)