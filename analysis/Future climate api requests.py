import cdsapi

# CMIP5 ------------------------------------------------------------------------------------------------

# https://cds.climate.copernicus.eu/datasets/projections-cmip5-monthly-single-levels?tab=overview
dataset = "projections-cmip5-monthly-single-levels"
client = cdsapi.Client()

#Create request first, then run:
#client.retrieve(dataset, request).download()


request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "access1_0",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "access1_3",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "bcc_csm1_1",
    "ensemble_member": "r1i1p1",
    "period": ["200601-209912"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "bcc_csm1_1_m",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "bnu_esm",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": ["mean_precipitation_flux"],
    "model": "ccsm4",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "cesm1_bgc",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "cesm1_cam5",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "cmcc_cms",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-200912",
        "201001-201912",
        "202001-202912",
        "203001-203912",
        "204001-204912",
        "205001-205912",
        "206001-206912",
        "207001-207912",
        "208001-208912",
        "209001-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "cnrm_cm5",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-205512",
        "205601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "csiro_mk3_6_0",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "fgoals_g2",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-201512",
        "201601-202512",
        "202601-203512",
        "203601-204512",
        "204601-205512",
        "205601-206512",
        "206601-207512",
        "207601-208512",
        "208601-209512",
        "209601-210512"
        ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "fio_esm",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "gfdl_cm3",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-201012",
        "201101-201512",
        "201601-202012",
        "202101-202512",
        "202601-203012",
        "203101-203512",
        "203601-204012",
        "204101-204512",
        "204601-205012",
        "205101-205512",
        "205601-206012",
        "206101-206512",
        "206601-207012",
        "207101-207512",
        "207601-208012",
        "208101-208512",
        "208601-209012",
        "209101-209512",
        "209601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "gfdl_esm2g",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-201012",
        "201101-201512",
        "201601-202012",
        "202101-202512",
        "202601-203012",
        "203101-203512",
        "203601-204012",
        "204101-204512",
        "204601-205012",
        "205101-205512",
        "205601-206012",
        "206101-206512",
        "206601-207012",
        "207101-207512",
        "207601-208012",
        "208101-208512",
        "208601-209012",
        "209101-209512",
        "209601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "gfdl_esm2m",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-201012",
        "201101-201512",
        "201601-202012",
        "202101-202512",
        "202601-203012",
        "203101-203512",
        "203601-204012",
        "204101-204512",
        "204601-205012",
        "205101-205512",
        "205601-206012",
        "206101-206512",
        "206601-207012",
        "207101-207512",
        "207601-208012",
        "208101-208512",
        "208601-209012",
        "209101-209512",
        "209601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "giss_e2_h",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-205012",
        "205101-210012",
        "210101-215012",
        "215101-220012",
        "220101-225012",
        "225101-230012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "giss_e2_h_cc",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-205012",
        "205101-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "giss_e2_r",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-202512",
        "202601-205012",
        "205101-207512",
        "207601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "giss_e2_r_cc",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-202512",
        "202601-205012",
        "205101-207512",
        "207601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "hadgem2_cc",
    "ensemble_member": "r1i1p1",
    "period": [
        "200512-203011",
        "203012-205511",
        "205512-208011",
        "208012-209912",
        "210001-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "hadgem2_es",
    "ensemble_member": "r1i1p1",
    "period": [
        "200512-203011",
        "203012-205511",
        "205512-208011",
        "208012-209911",
        "209912-212411"
        ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "inmcm4",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "ipsl_cm5a_lr",
    "ensemble_member": "r1i1p1",
    "period": ["200601-230012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "ipsl_cm5a_mr",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "ipsl_cm5b_lr",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "mpi_esm_lr",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "mpi_esm_mr",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210012"]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "noresm1_m",
    "ensemble_member": "r1i1p1",
    "period": [
        "200601-210012"
    ]
}
client.retrieve(dataset, request).download()

request = {
    "experiment": ["rcp_4_5", "rcp_8_5"],
    "variable": [
        "2m_temperature",
        "mean_precipitation_flux"
    ],
    "model": "noresm1_me",
    "ensemble_member": "r1i1p1",
    "period": ["200601-210212"]
}
client.retrieve(dataset, request).download()



# CMIP6 ------------------------------------------------------------------------------------------------


dataset = "projections-cmip6"
request = {
    "temporal_resolution": "monthly",
    "experiment": ["ssp5_8_5"],              #"ssp2_4_5", 
    "variable": "air_temperature",
    "level": ["1"],
    "model": ["access_cm2"], #,"awi_cm_1_1_mr","bcc_csm2_mr", "canesm5", "canesm5_canoe", "cesm2", 
             # "ciesm", "cmcc_cm2_sr5", "cmcc_esm2", "cnrm_cm6_1", "cnrm_cm6_1_hr", "cnrm_esm2_1",
             # "ec_earth3_cc", "ec_earth3_veg_lr", "fgoals_f3_l", "fgoals_g3", "fio_esm_2_0", 
             # "gfdl_esm4", "hadgem3_gc31_ll", "hadgem3_gc31_mm", "iitm_esm", "inm_cm4_8", "inm_cm5_0",
             # "ipsl_cm6a_lr", "kiost_esm", "mcm_ua_1_0", "miroc6", "miroc_es2l", "mpi_esm1_2_lr",
             # "mri_esm2_0", "nesm3", "noresm2_lm", "noresm2_mm", "taiesm1", "ukesm1_0_ll"],
    "month": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12"
    ],
    "year": [
        "2015", "2016", "2017",
        "2018", "2019", "2020",
        "2021", "2022", "2023",
        "2024", "2025", "2026",
        "2027", "2028", "2029",
        "2030", "2031", "2032",
        "2033", "2034", "2035",
        "2036", "2037", "2038",
        "2039", "2040", "2041",
        "2042", "2043", "2044",
        "2045", "2046", "2047",
        "2048", "2049", "2050",
        "2051", "2052", "2053",
        "2054", "2055", "2056",
        "2057", "2058", "2059",
        "2060", "2061", "2062",
        "2063", "2064", "2065",
        "2066", "2067", "2068",
        "2069", "2070", "2071",
        "2072", "2073", "2074",
        "2075", "2076", "2077",
        "2078", "2079", "2080",
        "2081", "2082", "2083",
        "2084", "2085", "2086",
        "2087", "2088", "2089",
        "2090", "2091", "2092",
        "2093", "2094", "2095",
        "2096", "2097", "2098",
        "2099", "2100"
    ],
    "area": [15, 49.8, 14, 50.05]
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()


