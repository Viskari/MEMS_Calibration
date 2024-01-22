Files:

	1. MEMS1K_SOC_FullStack.tif
		processing: 	sum each band of the geotiff files: MAOM_MEMS1k.tif and POM_MEMS1k.tif
		script: 		/_postp/_sum_geotiff_by_band.sh
		n bands:		41
		start year:	1980


	2. MEMS1K_SOC_2018-2009.tif
		
		processing: 	subtraction of the layers: 2018 (Band 39) and 2009 (Band 30)
		tool:			QGIS (raster calculator).
		formula: 		MEMS1K_SOC_2018-2009.tif = "MEMS1K_SOC_FullStack@39" - "MEMS1K_SOC_FullStack@30"