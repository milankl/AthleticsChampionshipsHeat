using FTPClient, Dates

now_string() = Dates.format(now(),Dates.RFC1123Format)

function synchronize_ftpdata(   hostname::String,
                                username::String,
                                password::String,
                                ftppath::String,
                                localpath::String)

    ftp = FTP(;hostname,username,password)
    success = false         # set to true once syncing completed

    try
        cd(ftp,ftppath)
        println("----")
        println("FTP Connection to $hostname establised $(now_string())")
        allfolders_on_ftp = readdir(ftp)
        allfolders_local = readdir(localpath)
        
        n_newfiles = 0

        for ftp_folder in allfolders_on_ftp[end:-1:1]

            println("Synchronizing $ftp_folder:")
            n_newfiles = 0

            # create local folder if not existent
            ftp_folder in allfolders_local ? nothing : mkdir(joinpath(localpath,ftp_folder))

            # compare number of files on ftp vs locally
            cd(ftp,ftp_folder)

            allfiles_on_ftp = readdir(ftp)
            allfiles_local = readdir(joinpath(localpath,ftp_folder))

            for file in allfiles_on_ftp
                if ~(file in allfiles_local)
                    n_newfiles += 1
                    download(ftp,file,joinpath(localpath,ftp_folder,file))
                end
            end

            println(" $n_newfiles downloaded.")

            # move back one folder up
            cd(ftp,"..")
        end

        # if some files were downloaded check again in shorter time
        success = n_newfiles == 0

    finally
        close(ftp)
        println("Connection closed: $(now_string())")
    end

    return success
end