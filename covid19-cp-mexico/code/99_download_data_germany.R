# This file can be used to download the current RKI data
# We used it on 2021-01-05 to download the current dataset and safed it in the data/raw folder
file_path <- file.path("..","data/raw/", str_c("RKI-",Sys.Date(),  ".csv"))
if (!file.exists(file_path)) {   
  url_RKI_dashboard_data <- "https://www.arcgis.com/sharing/rest/content/items/f10774f1c63e40168479a1feb6c7ca74/data"
  download.file(url=url_RKI_dashboard_data, destfile=file_path, mode="wb")
}
