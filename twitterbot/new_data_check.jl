pngismeteogram(filename::String,loc::String) = occursin(loc*r"_\d+.png",filename)
meteograms_created(path::String,loc::String) = [pngfile for pngfile in readdir(path) if pngismeteogram(pngfile,loc)]
isdatafolder(folder::String) = occursin(r"2021\d+",folder)

function data_feedback(new_data::Bool,timestr::String,loc::String)
    if new_data
        println("New data found for $loc at $timestr.")
    else
        println("No new data found for $loc.")
    end
end
    
function new_data_check(path::String,location::String)
    
    all = readdir(path)
    alldatafolders = [folder for folder in all if isdatafolder(folder)]
    allmeteograms = meteograms_created(joinpath(path,"meteograms"),location)
    
    mmdd_meteograms = [split(split(meteogram,"_")[2],".")[1][1:4] for meteogram in allmeteograms]
    
    new_data = false
    timestr = "00000000"
    
    for datafolder in alldatafolders
        
        # check that at least 270 files are in the folder, otherwise ignore
        nfiles = length(readdir(joinpath(path,datafolder)))
        
        if nfiles >= 270
        
            # extract month & day
            # yyyy = datafolder[1:4]
            mmdd = datafolder[5:8]
            
            if ~(mmdd in mmdd_meteograms)
                
                new_data = true
                timestr = mmdd*"0000"
                
                data_feedback(new_data,timestr,location)
                return new_data,timestr
            end
        end
    end
    
    data_feedback(new_data,timestr,location)
    return new_data,timestr
end
        
    