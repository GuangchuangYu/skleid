##' print Info
##'
##' 
##' @title printInfo
##' @return NULL
##' @author ygc
##' @export
printInfo <- function() {
    cat('\n')
    cat('         /-S\n')
    cat('      /-|\n')
    cat('   /-|   \\-K\n')
    cat('  |  |\n')
    cat('  |   \\----L\n')
    cat('--|\n')
    cat('  |      /-E\n')
    cat('  |   /-|\n')
    cat('   \\-|   \\-I\n')
    cat('     |\n')
    cat('      \\----D\n')  
    
    cat("\n###################################################################\n")
    cat("##                                                               ##\n")
    cat("##                                                               ##\n")
    ## cat("##\t\t\t\t\t\t\t\t ##\n")
    cat("##  Author: Guangchuang Yu (gcyu@connect.hku.hk)                 ##\n")
    cat(paste("##  skleid package, version = ", as.character(packageVersion("skleid")), 
              sep = ""), "                             ##\n")
    cat("##  use help(package=\"skleid\") to view online manuals            ##\n")
    cat("##  This package is designed for internal use of SKLEID          ##\n")
    cat("##  If you got any problem, please contact me by email           ##\n")
    cat("##                                                               ##\n")
    cat("##  For first time user:                                         ##\n")
    cat("##      Please send your email to me,                            ##\n")
    cat("##      so that I can keep you updated.                          ##\n")    
    cat("##                                                               ##\n")
    cat("##                                                               ##\n")    
    cat("###################################################################\n\n")
}

says <- function() {
    n <- 4
    i <- sample.int(n, size=1)
    if (i == 1) {
        mooseSay()
    } else if (i == 2) {
        cowSay()
    } else if (i == 3) {
        catSay()
    } else if (i == 4) {
        daemonSay()
    }
}

mooseSay <- function() {
    cat("  \\\n")
    cat("   \\   \\_\\_    _/_/\n")
    cat("    \\      \\__/\n")
    cat("           (oo)\\_______\n")
    cat("           (__)\\       )\\/\\\n")
    cat("               ||----w |\n")
    cat("               ||     ||\n")
}

cheeseSay <- function() {
    cat("   \\\n")
    cat("    \\\n")
    cat("      _____   _________\n")
    cat("     /     \\_/         | \n")
    cat("    |                 ||\n")
    cat("    |                 ||\n")
    cat("   |    ###\\  /###   | |\n")
        cat("   |     0  \\/  0    | |\n")
    cat("  /|                 | |\n")
    cat(" / |        <        |\\ \\\n")
    cat("| /|                 | | |\n")
    cat("| |     \\_______/   |  | |\n")
    cat("| |                 | / /\n")
    cat("/||                 /|||\n")
    cat("   ----------------|\n")
    cat("        | |    | |\n")
    cat("        ***    ***\n")
    cat("       /___\\  /___\\\n")
    
}

cowSay <- function() {
    cat("       \\   ^__^\n")
    cat("        \\  (oo)\\_______\n")
    cat("           (__)\\       )\\/\\\n")
    cat("               ||----w |\n")
    cat("               ||     ||\n")
}

dragonSay <- function() {
    cat("      \\                    / \\  //\\\n")
    cat("       \\    |\\___/|      /   \\//  \\\\\n")
    cat("            /0  0  \\__  /    //  | \\ \\    \n")
    cat("           /     /  \\/_/    //   |  \\  \\\n")  
    cat("           @_\\^_@\'\\/   \\/_   //    |   \\   \\\n")
    cat("           //_^_/     \\/_ //     |    \\    \\\n")
    cat("        ( //) |        \\///      |     \\     \\\n")
    cat("      ( / /) _|_ /   )  //       |      \\     _\\\n")
    cat("    ( // /) '/,_ _ _/  ( ; -.    |    _ _\\.-~        .-~~~^-.\n")
    cat("  (( / / )) ,-{        _      `-.|.-~-.           .~         `.\n")
    cat(" (( // / ))  '/\      /                 ~-. _ .-~      .-~^-.  \\\n")
    cat(" (( /// ))      `.   {            }                   /      \\  \\\n")
    cat("  (( / ))     .----~-.\\        \\-'                 .~         \\  `. \\^-.\n")
    cat("             ///.----..>        \\             _ -~             `.  ^-`  ^-_\n")
    cat("               ///-._ _ _ _ _ _ _}^ - - - - ~                     ~-- ,.-~\n")
    cat("                                                                  /.-~\n")
}

daemonSay <- function() {
    cat("   \\         ,        ,\n")
    cat("    \\       /(        )`\n")
    cat("     \\      \\ \\___   / |\n")
    cat("            /- _  `-/  '\n")
    cat("           (/\\/ \\ \\   /\\\n")
    cat("           / /   | `    \\\n")
    cat("           O O   ) /    |\n")
    cat("           `-^--'`<     '\n")
    cat("          (_.)  _  )   /\n")
    cat("           `.___/`    /\n")
    cat("             `-----' /\n")
    cat("<----.     __ / __   \\\n")
    cat("<----|====O)))==) \\) /====\n")
    cat("<----'    `--' `.__,' \\\n")
    cat("             |        |\n")
    cat("              \\       /\n")
    cat("        ______( (_  / \\______\n")
    cat("      ,'  ,-----'   |        \\\n")
    cat("      `--{__________)        \\/\n")
}


catSay <- function() {
    cat("    \\\n")
    cat("      \\\n")
    cat("        \\\n")
    cat("            |\\___/|\n")
    cat("          ==) ^Y^ (==\n")
    cat("            \\  ^  /\n")
    cat("             )=*=(\n")
    cat("            /     \\\n")
    cat("            |     |\n")
    cat("           /| | | |\\\n")
    cat("           \\| | |_|/\\\n")
    cat("           //_// ___/\n")
    cat("               \\_)\n")
}
