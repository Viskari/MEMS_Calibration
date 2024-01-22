#!/usr/bin/bash
#===============================================================================
#
#   FILE:  				_subtract_geotiff_by_year.sh
#   USAGE:  			./_subtract_geotiff_by_year.sh
#
#   DESCRIPTION: 		subtract two rasters by band number
#
#	INPUT PARAMETERS:
#		start_year:			first year of the dataset (band = 1)
#		bands:				number of bands (use gdalinfo to the this info)
#		reference_band:		number of the band related to the reference year
# 		reference_year:		reference year. The rasters will be subtracted to this layer
#		calc_year_start:	first year of the calculation
#		calc_year_end: 		last year of the calculation
#		no_data:			no data value (use gdalinfo to the this info)

#
#	FOLDERS: (can be relative or absolute paths)
#		input_folder:		folder of the input raster file  
#		input_raster:		name of the input raster file  
#		output_folder:		folder of the output file 
#		output_file:		name of the output file
#		tmp_folder:			temporary folder for calculations (don`t need to be changed)
#
#	OTHER PARAMETERS:
# 		version:			suggest to use the running date in the following format: YY-MM-DD
#		metadata:			add additional metadata info (calculation and band info) ### don`t add blank spaces #
#
#   CREATED:  02/10/2023
#===============================================================================


# __________________ INPUT PARAMETERS
start_year=1980
bands=41
reference_band=30 # band of the year 2009
reference_year=2009
calc_year_start=2010 # first year of the calculation
calc_year_end=2018 	# last year of the calculation
no_data=-3.39999995214436425e+38


# __________________ FOLDERS
input_folder="../_sum"
input_raster="MEMS1K_SOC_FullStack.tif"

output_folder="../_subtract"
output_file="MEMS1K_SOC_Subtract.tif"

tmp_folder="${output_folder}/_tmp"

mkdir -p ${tmp_folder}

# __________________ Metadata
version="2023-10-02"
metadata="-mo VERSION=${version}"

# __________________ Processing (sum bands)

COUNTER=1
for band in $(seq 1 $bands)
do
	
	year=$((start_year -1 + $band))
	tmp_file="${tmp_folder}/${year}.tif"	
		
	
	if ((year >= calc_year_start && year <= calc_year_end))
	then
	
		echo "... year: ${year}   calc: ${band} - ${reference_band}    band: ${COUNTER}"	
		
		# _________ add raster by band number	
		gdal_calc.py -A "${input_folder}/${input_raster}" --A_band=$band -B "${input_folder}/${input_raster}" --B_band=${reference_band} --type Float32 --outfile=${tmp_file} --NoDataValue=${no_data} --calc="A-B" --overwrite
	
		# _________ create a metadata for each band with the year (printf to add leading zero to band name)
		metadata="${metadata} -mo BAND_$(printf %02d ${COUNTER})=${year}-${reference_year}"	
		let COUNTER++
	
	fi	


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
rm -Rf ${tmp_folder}
