# uses lots of hacks to try to estimate a reasonable amount of memory to use, in the most common systems!
# currently only linux is actually supported
get_mem_lim <- function(factor = 0.7, verbose = FALSE) {
    # make sure the factor makes sense
    if (factor <= 0)
        stop('memory `factor` must be strictly positive!  Passed: ', factor)
    if (factor > 1)
        stop('memory `factor` must be <= 1!  Passed: ', factor)
    
    mem <- NA # to know if we've succeeded or not...
    
    if ( .Platform$OS.type == 'unix' ) {
        # try this, will be NA if not linux
        mem <- parse_meminfo_linux()
        ## if ( is.na( mem ) ) {
        ##     # ... MAC OSX behavior goes here
        ##     # https://stackoverflow.com/questions/14150626/understanding-vm-stat-in-mac-os-how-to-convert-those-numbers-to-something-simil
        ## }
    } else {
        # only other option is WINDOWS
        # expanded from:
        # http://stackoverflow.com/questions/27788968/how-would-one-check-the-system-memory-available-using-r-on-a-windows-machine

        # started failing recently (2023-01-06), let's just use the default here too (as in MacOS, set further below)
        ## # get string with values from system
        ## mem <- system('wmic OS get FreePhysicalMemory /Value', intern = TRUE)[3]
        ## # remove carriage return
        ## mem <- gsub("\r",'', mem)
        ## # split value and keep number only
        ## mem <- as.numeric(strsplit(mem, '=')[[1]][2])
        ## mem <- mem*1024 # previous units were KB, convert to bytes
    }
    if (is.na(mem)) {
        # I used to have a warning, but not anymore because:
        # 1) It's too annoying (tests fail when they really shouldn't)
        # 2) This issue (of not really knowing the available memory) is not limiting (except on the largest numbers of individuals) now that popkin defaults to using the least memory possible anyway
        #warning("Could not infer available memory, will default to 1GB!\nPlease specify a memory limit if your run exceeds memory or if default is too low!")
        mem <- GB # when we can't determine free or available memory from system, default to using 1GB!
    } else {
        mem <- mem * factor # shrink memory by a factor to leave some more memory
    }
    
    if (verbose)
        message('Will limit mem to about ', round( mem/GB, 2 ), ' GB') # DEBUGGING
    
    return(mem) # return default or better value if it was available
}

# though it actually fully parses the memory table, it only returns available memory for now
# originally based on this, but since adapted to use MemAvailable, or include as fallback use free+buff+cache
# https://stackoverflow.com/questions/6457290/how-to-check-the-amount-of-ram-in-r
parse_meminfo_linux <- function() {
    # if this file doesn't exist, return NA
    # fails in non-linuxes, but shouldn't fail in any linux
    file_meminfo <- '/proc/meminfo'
    if ( !file.exists( file_meminfo ) )
        return( NA )
    # else continue
    
    # parse into a nice data frame
    data <- utils::read.table(
                       file_meminfo,
                       fill = TRUE, # needed since some rows are missing the unit
                       col.names = c('type', 'value', 'unit')
                   )
    
    # clean names by removing final semicolons
    data$type <- sub( ':$', '', data$type )

    # just for niceness, autocomplete trivial missing units
    # NOTE: in my computer all non-missing units are kB, but didn't want to assume that is always true...
    indexes <- data$unit == ''
    # let's leave as is if missing units don't all have values of zero
    if ( any( indexes ) && all( data$value[ indexes ] == 0 ) ) { # true on my machine
        # extract non-missing units
        units <- data$unit[ !indexes ]
        # in my machine there's only one value
        units <- table( units )
        # this recovers the most common unit (in my case "kB")
        units <- names( units[ which.max( units ) ] )
        # assign that unit to the missing cases
        data$unit[ indexes ] <- units
    }

    # the goal is usually to get available memory: MemAvailable
    # (this is oddly not memFree, which can be much smaller in some cases and its meaning is not what we desire.)
    # https://superuser.com/questions/980820/what-is-the-difference-between-memfree-and-memavailable-in-proc-meminfo
    # sadly, MemAvailable doesn't exist in all systems (added to kernel 3.14 in 2014), so we've got to have a fallback
    # https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git/commit/?id=34e431b0ae398fc54ea69ff85ec700722c9da773
    mem_avail <- get_mem_by_type( data, 'MemAvailable', fatal = FALSE )
    if ( !is.na( mem_avail ) )
        return( mem_avail )
    # else continue

    # estimate available memory using fallback procedure
    # this doesn't give exactly the same answer on my system, but it's close enough
    # these have all been in linux for longer so I don't expect any of them to be missing (all are mandatory)
    mem_free <- get_mem_by_type( data, 'MemFree' )
    mem_buff <- get_mem_by_type( data, 'Buffers' )
    mem_cache <- get_mem_by_type( data, 'Cached' )
    # done with fallback procedure, returning free+buff+cache
    return( mem_free + mem_buff + mem_cache )
}

# used by: parse_meminfo_linux
get_mem_by_type <- function( data, type, fatal = TRUE ) {
    index <- which( data$type == type )
    # if this fails, index is `integer(0)`
    if ( length( index ) == 0 ) {
        if ( fatal ) {
            stop( 'Unexpected: `', type, '` not found!' )
        } else
            return( NA ) # this means it wasn't found
    }
    if ( length( index ) != 1 )
        stop( 'Unexpected: more than one row match for: ', type )
    # get value, scaled correctly by unit
    mem_type <- scale_unit( data$value[ index ], data$unit[ index ] )
    return( mem_type )
}

# used by: parse_meminfo_linux
scale_unit <- function( value, unit ) {
    # scale by unit, convert to bytes
    if ( unit == 'kB' ) {
        value <- value * 1024
    } else
        warning( 'Unhandled unit (parsing /proc/meminfo): ', unit )
}
