# uses lots of hacks to try to estimate a reasonable amount of memory to use, in the most common systems!
# currently only linux and windows are actually supported
get_mem_lim <- function(factor = 0.7, verbose = FALSE) {
    # make sure the factor makes sense
    if (factor <= 0)
        stop('memory `factor` must be strictly positive!  Passed: ', factor)
    if (factor > 1)
        stop('memory `factor` must be <= 1!  Passed: ', factor)
    
    mem <- NA # to know if we've succeeded or not...
    
    if ( .Platform$OS.type == 'unix' ) {
        if (file.exists('/proc/meminfo')) {
            # expanded from:
            # https://stackoverflow.com/questions/6457290/how-to-check-the-amount-of-ram-in-r
            
            # this reads file as a table with a single column
            meminfo <- utils::read.table('/proc/meminfo', sep = "\n", stringsAsFactors = FALSE)
            meminfo <- meminfo[[1]] # reduce data frame to vector
            # these annoying commands extract the two numbers of interest (MemTotal and MemFree)
            # mem2 <- as.numeric(strsplit(meminfo[grepl('MemTotal', meminfo)], split=' +')[[1]][2])
            mem <- as.numeric(strsplit(meminfo[grepl('MemFree', meminfo)], split = ' +')[[1]][2])
            mem <- mem*1024 # previous units were KB, convert to bytes
        } else {
            # ... MAC OSX behavior goes here
            # https://stackoverflow.com/questions/14150626/understanding-vm-stat-in-mac-os-how-to-convert-those-numbers-to-something-simil
        }
    } else {
        # only other option is WINDOWS
        # expanded from:
        # http://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine
        
        # get string with values from system
        mem <- system('wmic OS get FreePhysicalMemory /Value', intern = TRUE)[3]
        # remove carriage return
        mem <- gsub("\r",'', mem)
        # split value and keep number only
        mem <- as.numeric(strsplit(mem, '=')[[1]][2])
        mem <- mem*1024 # previous units were KB, convert to bytes
    }
    if (is.na(mem)) {
        warning("Could not infer available memory, will default to 1GB!\nPlease specify a memory limit if your run exceeds memory or if default is too low!")
        mem <- GB # when we can't determine free or available memory from system, default to using 1GB!
    } else {
        mem <- mem * factor # shrink memory by a factor to leave some more memory
    }
    
    if (verbose)
        message('Will limit mem to about ', round( mem/GB, 2 ), ' GB') # DEBUGGING
    
    return(mem) # return default or better value if it was available
}
