## uses lots of hacks to try to estimate a reasonable amount of memory to use, in the most common systems!
## currently only linux is actually supported
getMemLim <- function(factor=0.7) {
    M <- NA # to know if we've succeeded or not...
    
    if ( .Platform$OS.type == 'unix' ) {
        if (file.exists('/proc/meminfo')) {
            ## expanded from:
            ## https://stackoverflow.com/questions/6457290/how-to-check-the-amount-of-ram-in-r
            
            ## this reads file as a table with a single column
            meminfo <- read.table('/proc/meminfo', sep="\n", stringsAsFactors=FALSE)
            meminfo <- meminfo[[1]] ## reduce data frame to vector
            ## these annoying commands extract the two numbers of interest (MemTotal and MemFree)
            ## M2 <- as.numeric(strsplit(meminfo[grepl('MemTotal', meminfo)], split=' +')[[1]][2])
            M <- as.numeric(strsplit(meminfo[grepl('MemFree', meminfo)], split=' +')[[1]][2])
            ## M <- min(M2/2, M) # use the least of the free memory or half the total memory
            M <- M*1024 # previous units were KB, convert to bytes
        } else {
            ## ... MAC OSX behavior goes here
            ## https://stackoverflow.com/questions/14150626/understanding-vm-stat-in-mac-os-how-to-convert-those-numbers-to-something-simil
        }
    } else {
        ## only other option is WINDOWS
        ## expanded from:
        ## http://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
        
        ## get string with values from system
        M <- system('wmic OS get FreePhysicalMemory /Value', intern=TRUE)[3]
        ## remove carriage return
        M <- gsub("\r",'', M)
        ## split value and keep number only
        M <- as.numeric(strsplit(M, '=')[[1]][2])
        M <- M*1024 # previous units were KB, convert to bytes
    }
    if (is.na(M)) {
        warning("Could not infer available memory, will default to 2GB!\nPlease specify a memory limit if your run exceeds memory or if default is too low!")
        M <- 2*1024*1024*1024 # when we can't determine free or available memory from system, default to using 2GB!
    } else {
        M <- M*factor # shrink memory by a factor to leave some more memory
    }
    message('Will limit mem to ', round( M/(1024*1024*1024), 2 ), ' GB') # DEBUGGING
    return(M) # return default or better value if it was available
}
