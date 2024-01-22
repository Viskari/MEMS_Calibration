cd /eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_merge

gdalbuildvrt maom.vrt /eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_out/MAOM*.tif
gdal_translate maom.vrt MAOM_MEMS1k.tif -co COMPRESS=LZW
gdalbuildvrt  pom.vrt /eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas1/_out/POM*.tif
gdal_translate pom.vrt POM_MEMS1k.tif -co COMPRESS=LZW



