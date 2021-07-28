using PyCall, JSON

py"""
import sys
sys.path.append($pwd())
"""

create_meteogram = pyimport("create_meteogram")["create_meteogram"]

include("ftp_download.jl")
include("new_data_check.jl")

localpath = "/network/aopp/chaos/pred/kloewer/tokyo2021"
ecmwf = JSON.parsefile("ecmwf.json")
stopdate = DateTime(2021,08,09)

while Dates.now() < stopdate    # repeat till Olympic Games are over
    
    # check for new data on ftp server
    download_success = synchronize_ftpdata(ecmwf["hostname"],
                                            ecmwf["username"],
                                            ecmwf["password"],
                                            ecmwf["ftppath"],
                                            localpath)
    
    # check whether data folders are complete
    for location in ["sapporo","tokyo"]
        new_data,timestr = new_data_check(localpath,location)
    
        # initiate meteogram creation and tweeting
        new_data ? create_meteogram(location,timestr) : nothing
    end

    # sync every hour unless unsucessful/incomplete then try again 20min later
    download_success ? sleep(3600) : sleep(1200)
end
