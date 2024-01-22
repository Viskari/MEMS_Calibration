#!/usr/bin/bash
#===============================================================================
#
#   FILE:  				_sum_geotiff_by_band.sh
#   USAGE:  			./_sum_geotiff_by_band.sh
#
#   DESCRIPTION: 		sum two raster by band number
#
#	INPUT PARAMETERS:
#		start_year:		first year of the dataset (band = 1)
#		bands:			number of bands (use gdalinfo to the this info)
#		raster1:		name of the 1st raster file
#		raster2:		name of the 2nd raster file
#		no_data:		no data value (use gdalinfo to the this info)
#		output_file:	name of the output file
#
#	FOLDERS: (can be relative or absolute paths)
#		geotiff_folder:	folder of the input (raster) files  
#		output_folder:	folder of the output file 
#		tmp_folder:		temporary folder for calculations (don`t need to be changed)
#
#	OTHER PARAMETERS:
# 		version:		suggest to use the running date in the following format: YY-MM-DD
#		metadata:		add additional metadata info (calculation and band info) ### don`t add blank spaces #
#
#   CREATED:  25/09/2023
#===============================================================================


# __________________ INPUT PARAMETERS
start_year=1980
bands=41
raster1="MAOM_MEMS1k.tif"
raster2="POM_MEMS1k.tif"
no_data=-3.39999995214436425e+38

output_file="MEMS1K_SOC_FullStack.tif"

# __________________ FOLDERS
geotiff_folder="../_merge"
output_folder="../_sum"
tmp_folder="${output_folder}/_tmp"

mkdir -p ${tmp_folder}

# __________________ Metadata
version="2023-09-25"
metadata="-mo CALCULATION=${raster1}+${raster2} -mo VERSION=${version}"

# __________________ Processing (sum bands)

for band in $(seq 1 $bands)
do
	
	year=$((start_year -1 + $band))
	tmp_file="${tmp_folder}/${year}.tif"	

	echo "... year: ${year}   band: ${band} tmp: ${tmp_file}"		
	
	# _________ add raster by band number	
	gdal_calc.py -A "${geotiff_folder}/${raster1}" --A_band=$band -B "${geotiff_folder}/${raster2}" --B_band=$band --type Float32 --outfile=${tmp_file} --NoDataValue=${no_data} --calc="A+B" --overwrite
	
	# _________ create a metadata for each band with the year (printf to add leading zero to band name)
	metadata="${metadata} -mo BAND_$(printf %02d ${band})=${year}"

done

# __________________ Create output file
echo " ... create virtual layer"
gdalbuildvrt -separate "${output_file}.vrt" ${tmp_folder}/*.tif

echo " ... create geotifff file: ${output_file}"
gdal_translate "${output_file}.vrt" "${output_folder}/${output_file}" -co COMPRESS=LZW

echo " ... add custom metadata"
gdal_edit.py ${metadata} "${output_folder}/${output_file}" 


# __________________ Clean up
rm "${output_file}.vrt"
#~ rm -Rf ${tmp_folder}
