
import os
import openeo

import geopandas as gpd

from pathlib import Path

def BBOX_to_geojson(BBOX):

    feature_collection = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[
                        [BBOX["west"], BBOX["south"]],
                        [BBOX["west"], BBOX["north"]],
                        [BBOX["east"], BBOX["north"]],
                        [BBOX["east"], BBOX["south"]],
                        [BBOX["west"], BBOX["south"]]
                    ]]
                }   
            }
        ]
    }

    return(feature_collection)

def add_fapar(con,BBOX):
    fapar = con.load_collection("TERRASCOPE_S2_FAPAR_V2",
                                spatial_extent=BBOX,
                                temporal_extent=["2020-01-01","2024-12-31"],
                                bands=["FAPAR_10M","SCENECLASSIFICATION_20M"])
    fapar_cube_masked = fapar.process("mask_scl_dilation",data=fapar, scl_band_name="SCENECLASSIFICATION_20M")
    fapar_cube_masked = fapar_cube_masked.filter_bands("FAPAR_10M")

    return(fapar_cube_masked)

def add_fapar_CROPSAR(con,BBOX,year):

    fapar = con.datacube_from_process(
        "CropSAR_px",
        namespace="vito",
        geometry = BBOX_to_geojson(BBOX),
        startdate = f"{year}-01-01",
        enddate = f"{year}-12-31",
        output = "FAPAR"
    )

    return(fapar)

def add_varmap(con,BBOX):

    varmap = con.datacube_from_process(
        "variability_map",
        namespace="vito",
        date=["2020-01-01","2024-12-31"],
        polygon = BBOX_to_geojson(BBOX)
        )
    
    return(varmap)

def add_yield_potential(con,BBOX,year):

    yield_potential = con.datacube_from_process(
        "yieldpotentialmap_shub",
        namespace="vito",
        check=False,
        date=[f"{str(year)}-01-01",f"{str(year)}-12-31"],
        polygon = BBOX_to_geojson(BBOX),
        )
    
    return(yield_potential)

def extract_openeo(shp,var,outfile,CROPSAR=False,overwrite=False,YPM=False):

    if not os.path.exists(outfile) or overwrite:

        con = openeo.connect("openeo.vito.be").authenticate_oidc()

        shp = shp.to_crs("EPSG:4326")

        west,south,east,north = shp.geometry.total_bounds

        BBOX = {
            'east': east,
            'south': south,
            'west': west,
            'north': north,
            'crs': "EPSG:4326"}
        
        if CROPSAR|YPM:
            years = [2020,2021,2022,2023,2024]
            for year in years:
                outfile_y = outfile.replace(".nc",f"_{year}.nc")
                if var == "fapar":
                    datacube = add_fapar_CROPSAR(con,BBOX,year)
                elif var=="yield_potential":
                    datacube = add_yield_potential(con,BBOX,year)

                job = datacube.execute_batch(
                    title=f"Extracting {var} for field",
                    description=f"Extracting {var} for field",
                    out_format = "netcdf"
                )

                results = job.get_results()
                for asset in results.get_assets():
                    asset.download(outfile_y)

        else:
            if var == "fapar":
                datacube = add_fapar(con,BBOX)

            elif var == "varmap":
                datacube = add_varmap(con,BBOX)

            job = datacube.execute_batch(
                    title=f"Extracting {var} for field",
                    description=f"Extracting {var} for field",
                    out_format = "netcdf"
            )

            results = job.get_results()

            for asset in results.get_assets():
                asset.download(outfile)

            
    
if __name__ == "__main__":

    field_folder = "/data/sigma/AdamPriscilla/varmaps_check/field_data/"
    field_name = "veld.shp"

    field = gpd.read_file(os.path.join(field_folder, field_name))

    #Extracting Fapar from OpenEO
    outfile = os.path.join(Path(field_folder).parent,"openeo_extractions",field_name.replace(".shp","")+"_fapar.nc")
    os.makedirs(Path(outfile).parent,exist_ok=True)
    extract_openeo(field,"fapar",outfile)

    #Also extracting CROPSAR FAPAR for comparison                       
    outfile = os.path.join(Path(field_folder).parent,"openeo_extractions",field_name.replace(".shp","")+"_fapar_cropsar.nc")
    extract_openeo(field,"fapar",outfile,CROPSAR=True)

    #Extracting Variability Map from OpenEO
    outfile = os.path.join(Path(field_folder).parent,"openeo_extractions",field_name.replace(".shp","")+"_varmap.nc")
    extract_openeo(field,"varmap",outfile)

    #Extracting Yield Potential from OpenEO
    outfile = os.path.join(Path(field_folder).parent,"openeo_extractions",field_name.replace(".shp","")+"_yield_potential.nc")
    extract_openeo(field,"yield_potential",outfile,YPM=True)

    


