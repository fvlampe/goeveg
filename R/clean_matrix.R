#' Clean species matrix
#' @description
#' The function cleans a species matrix by removing species without occurrence (frequency = 0) and samples without any species (species number = 0) in one step.
#'  
#' It will also subset the corresponding observations of environmental data (samples in rows) or species trait data (species in rows), if passed to the function.
#'
#' @param matrix Community data, a matrix-like object with samples in rows and species in columns. Missing values (NA) or empty character cell values will be transformed to 0.
#' @param env Optionally, a data frame of environmental variables, with samples in rows and variables in columns
#' @param traits Optionally, a data frame of species traits, with species in rows and trait variables in columns
#'
#' @return
#' If only a species matrix is provided, the return will be a cleaned species matrix. 
#' If environmental and/or trait data are also provided, the result will be a list of the cleaned and subsetted matrices/data frames. 
#'
#' @examples
#' # Clean species matrix
#' schedenveg.clean <- clean_matrix(schedenveg)
#' 
#' # Clean species matrix and subset environmental data
#' scheden.clean <- clean_matrix(schedenveg, schedenenv)
#' schedenveg.clean <- scheden.clean$matrix
#' schedenenv.clean <- scheden.clean$env
#' 
#' @author Friedemann von Lampe
#' @export
#' @import stats utils


clean_matrix <- function(matrix, env = NULL, traits = NULL) {
  
  matrix.clean <- NULL
  env.clean <- NULL
  traits.clean <- NULL
  
  # Check for NA values
  if (any(is.na(matrix))) {
    matrix[is.na(matrix)] <- 0
    print("NA values in matrix transformed into 0.")
  }
  
  # Check for empty character cells 
  if (any(matrix == "")) {
    matrix[matrix == ""] <- 0
    print("Empty cell values in matrix transformed into 0.")
  }
  
  # Calculate frequencies and richness
  freq <- apply(matrix != 0, 2, sum)
  freq
  freq.0 <- length(freq[freq == 0])  
  
  rich <- apply(matrix != 0, 1, sum)
  rich
  rich.0 <- length(rich[rich == 0])  
  
  # Remove empty species and samples
  matrix.clean <- matrix[rich > 0, freq > 0]
  
  # Select environmental data
  if(!is.null(env)) {
    env.clean <- env[rich > 0, ]
    env.clean <- droplevels(env.clean)   # Remove possible empty factor levels
    
  }
  
  # Select traits data
  if(!is.null(traits)) {
    traits.clean <- env[freq > 0, ]
  }
  
  print(paste0("Removed ", freq.0, " species without occurence and ", 
               rich.0, " samples without any species"))
  
  if(is.null(traits) && is.null(env)) {
    return(matrix.clean)
  } else {
    return(list(matrix = matrix.clean, env = env.clean, traits = traits.clean))
  }
}
